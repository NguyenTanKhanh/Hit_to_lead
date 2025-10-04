# -*- coding: utf-8 -*-
import os
import argparse
import pandas as pd
import numpy as np
from pathlib import Path

import torch
from chemprop import data, featurizers, models
from lightning import pytorch as pl

os.environ["PL_TRAINER_CLUSTER_ENV"] = "local"


def _build_dataset(smiles_list, num_workers=0):
    datapoints = [data.MoleculeDatapoint.from_smi(smi) for smi in smiles_list]
    featurizer = featurizers.SimpleMoleculeMolGraphFeaturizer()
    dataset = data.MoleculeDataset(datapoints, featurizer)
    loader = data.build_dataloader(dataset, shuffle=False, num_workers=num_workers)
    return dataset, loader


def _predict_single_ckpt(ckpt_path, dataset, num_workers=0):
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model: models.MPNN = models.MPNN.load_from_checkpoint(
        ckpt_path, map_location=device
    ).eval()
    loader = data.build_dataloader(dataset, shuffle=False, num_workers=num_workers)

    trainer = pl.Trainer(
        logger=False,
        enable_progress_bar=False,
        accelerator="gpu" if device == "cuda" else "cpu",
        devices=1,
    )
    preds_batches = trainer.predict(model, loader)
    preds = np.concatenate([b.cpu().numpy().flatten() for b in preds_batches])
    return preds


def run_three_stage_filter(
    smiles_list,
    cls_ckpt: str,
    reg1_ckpt: str = None,
    reg2_ckpt: str = None,
    cls_threshold: float = None,
    reg1_threshold: float = None,
    reg2_threshold: float = None,
    num_workers: int = 4,
    log_csv_path: str = None,
    out_smi_path: str = None,
):
    n = len(smiles_list)
    df_log = pd.DataFrame({"SMILES": smiles_list})

    if cls_ckpt is None or cls_threshold is None:
        raise ValueError("Classification checkpoint (--cls-ckpt) and threshold (--cls-threshold) are required.")

    ds_all, _ = _build_dataset(smiles_list, num_workers=num_workers)
    cls_preds = _predict_single_ckpt(cls_ckpt, ds_all, num_workers=num_workers).reshape(-1)
    df_log["cls_prob"] = cls_preds
    pass_cls_mask = cls_preds >= cls_threshold
    df_log["pass_cls"] = pass_cls_mask.astype(int)
    idx_cls = np.where(pass_cls_mask)[0]
    print(f"[OK] Classification pass: {len(idx_cls)} / {n} (threshold >= {cls_threshold})")

    if len(idx_cls) == 0:
        if log_csv_path:
            df_log.to_csv(log_csv_path, index=False)
        if out_smi_path:
            open(out_smi_path, "w").close()
        return []

    idx_final = idx_cls.copy()

    if reg1_ckpt is not None and reg1_threshold is not None:
        smiles_r1 = [smiles_list[i] for i in idx_final]
        ds_r1, _ = _build_dataset(smiles_r1, num_workers=num_workers)
        reg1_preds = _predict_single_ckpt(reg1_ckpt, ds_r1, num_workers=num_workers).reshape(-1)

        df_log["reg1_pred"] = np.nan
        df_log.loc[idx_final, "reg1_pred"] = reg1_preds

        pass_r1_local = reg1_preds > reg1_threshold
        pass_r1_global = np.zeros(n, dtype=int)
        pass_r1_global[idx_final[pass_r1_local]] = 1
        df_log["pass_reg1"] = pass_r1_global

        idx_final = idx_final[pass_r1_local]
        print(f"[OK] Regression#1 pass: {len(idx_final)} molecules (threshold > {reg1_threshold})")

    if reg2_ckpt is not None and reg2_threshold is not None and len(idx_final) > 0:
        smiles_r2 = [smiles_list[i] for i in idx_final]
        ds_r2, _ = _build_dataset(smiles_r2, num_workers=num_workers)
        reg2_preds = _predict_single_ckpt(reg2_ckpt, ds_r2, num_workers=num_workers).reshape(-1)

        df_log["reg2_pred"] = np.nan
        df_log.loc[idx_final, "reg2_pred"] = reg2_preds

        pass_r2_local = reg2_preds > reg2_threshold
        pass_r2_global = np.zeros(n, dtype=int)
        pass_r2_global[idx_final[pass_r2_local]] = 1
        df_log["pass_reg2"] = pass_r2_global

        idx_final = idx_final[pass_r2_local]
        print(f"[OK] Regression#2 pass: {len(idx_final)} molecules (threshold > {reg2_threshold})")

    df_log["final_pass"] = 0
    df_log.loc[idx_final, "final_pass"] = 1
    final_smiles = [smiles_list[i] for i in idx_final]
    print(f"[FINAL] Selected: {len(final_smiles)} / {n}")

    if log_csv_path:
        Path(os.path.dirname(log_csv_path) or ".").mkdir(parents=True, exist_ok=True)
        df_log.to_csv(log_csv_path, index=False)
        print(f"[SAVED] Log to {log_csv_path}")

    if out_smi_path:
        Path(os.path.dirname(out_smi_path) or ".").mkdir(parents=True, exist_ok=True)
        with open(out_smi_path, "w") as f:
            for smi in final_smiles:
                f.write(smi + "\n")
        print(f"[SAVED] {len(final_smiles)} SMILES to {out_smi_path}")

    return final_smiles


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Three-stage filtering (classification + optional regressions).")
    ap.add_argument("--in-smi", required=True, help="Input SMILES file (.smi)")
    ap.add_argument("--out-csv", required=True, help="Output CSV log file")
    ap.add_argument("--out-smi", required=True, help="Output SMILES file with passed molecules")
    ap.add_argument("--cls-ckpt", required=True, help="Classification checkpoint")
    ap.add_argument("--cls-threshold", type=float, required=True, help="Classification threshold (prob >= value)")
    ap.add_argument("--reg1-ckpt", default=None, help="Regression #1 checkpoint (or None)")
    ap.add_argument("--reg1-threshold", type=float, default=None, help="Regression #1 threshold (pred > value)")
    ap.add_argument("--reg2-ckpt", default=None, help="Regression #2 checkpoint (or None)")
    ap.add_argument("--reg2-threshold", type=float, default=None, help="Regression #2 threshold (pred > value)")
    ap.add_argument("--num-workers", type=int, default=4, help="Number of workers for DataLoader")

    args = ap.parse_args()
    smiles_list = [ln.strip() for ln in open(args.in_smi) if ln.strip()]

    run_three_stage_filter(
        smiles_list,
        cls_ckpt=None if args.cls_ckpt == "None" else args.cls_ckpt,
        reg1_ckpt=None if args.reg1_ckpt == "None" else args.reg1_ckpt,
        reg2_ckpt=None if args.reg2_ckpt == "None" else args.reg2_ckpt,
        cls_threshold=args.cls_threshold,
        reg1_threshold=args.reg1_threshold,
        reg2_threshold=args.reg2_threshold,
        num_workers=args.num_workers,
        log_csv_path=args.out_csv,
        out_smi_path=args.out_smi,
    )
