# -*- coding: utf-8 -*-
import os
import argparse
import pandas as pd
import numpy as np
from pathlib import Path

import torch
from chemprop import data, featurizers, models
from lightning import pytorch as pl


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
        strategy="auto",   # tránh detect MPI/Slurm
    )
    preds_batches = trainer.predict(model, loader)
    preds = np.concatenate([b.cpu().numpy().flatten() for b in preds_batches])
    return preds


def run_pka_filter(
    smiles_list,
    reg_ckpt: str,
    reg_min: float = None,
    reg_max: float = None,
    num_workers: int = 4,
    log_csv_path: str = None,
    out_smi_path: str = None,
):
    n = len(smiles_list)
    df_log = pd.DataFrame({"SMILES": smiles_list})

    # --- Regression prediction ---
    ds, _ = _build_dataset(smiles_list, num_workers=num_workers)
    preds = _predict_single_ckpt(reg_ckpt, ds, num_workers=num_workers).reshape(-1)

    df_log["pKa_pred"] = preds

    # --- Apply filters ---
    pass_mask = np.ones_like(preds, dtype=bool)
    if reg_min is not None:
        pass_mask &= preds >= reg_min
    if reg_max is not None:
        pass_mask &= preds <= reg_max

    df_log["pass"] = pass_mask.astype(int)

    idx_final = np.where(pass_mask)[0]
    final_smiles = [smiles_list[i] for i in idx_final]

    if reg_min is not None and reg_max is not None:
        print(f"[OK] pKa pass: {len(final_smiles)} / {n} (range [{reg_min}, {reg_max}])")
    elif reg_min is not None:
        print(f"[OK] pKa pass: {len(final_smiles)} / {n} (>= {reg_min})")
    elif reg_max is not None:
        print(f"[OK] pKa pass: {len(final_smiles)} / {n} (<= {reg_max})")
    else:
        print("[WARN] No threshold specified, all molecules pass.")

    # --- Xu?t file ---
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
    ap = argparse.ArgumentParser(description="Filter SMILES by regression model (pKa range).")
    ap.add_argument("--in-smi", required=True, help="Input SMILES file (.smi)")
    ap.add_argument("--out-csv", required=True, help="Output CSV log file")
    ap.add_argument("--out-smi", required=True, help="Output SMILES file with passed molecules")
    ap.add_argument("--reg-ckpt", required=True, help="Regression model checkpoint")
    ap.add_argument("--reg-min", type=float, default=None, help="Minimum cutoff (pred >= value)")
    ap.add_argument("--reg-max", type=float, default=None, help="Maximum cutoff (pred <= value)")
    ap.add_argument("--num-workers", type=int, default=4, help="Number of workers for DataLoader")

    args = ap.parse_args()
    smiles_list = [ln.strip() for ln in open(args.in_smi) if ln.strip()]

    run_pka_filter(
        smiles_list,
        reg_ckpt=args.reg_ckpt,
        reg_min=args.reg_min,
        reg_max=args.reg_max,
        num_workers=args.num_workers,
        log_csv_path=args.out_csv,
        out_smi_path=args.out_smi,
    )
