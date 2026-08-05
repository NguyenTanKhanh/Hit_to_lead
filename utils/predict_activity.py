# -*- coding: utf-8 -*-
import os

# Prevent Lightning from auto-detecting MPI/SLURM environments.
# These variables must be set before importing Lightning.
os.environ["PL_TRAINER_CLUSTER_ENV"] = "local"
os.environ["LIGHTNING_CLOUD_DISABLE"] = "1"

import argparse
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd
import torch
from rdkit import Chem
from chemprop import data, featurizers, models
from lightning import pytorch as pl
from lightning.fabric.plugins.environments import LightningEnvironment


def _validate_smiles(
    smiles_list: List[str],
) -> Tuple[List[str], np.ndarray, List[str]]:
    """
    Validate and canonicalize SMILES.

    Returns
    -------
    valid_smiles
        Canonical SMILES accepted by RDKit.
    valid_indices
        Original row indices corresponding to valid_smiles.
    errors
        One entry per input row. Empty string means valid.
    """
    valid_smiles = []
    valid_indices = []
    errors = [""] * len(smiles_list)

    for index, raw_smi in enumerate(smiles_list):
        smi = str(raw_smi).strip()

        if not smi:
            errors[index] = "Empty SMILES"
            continue

        mol = Chem.MolFromSmiles(smi)

        if mol is None:
            errors[index] = "Invalid SMILES"
            print(
                f"[WARNING] Invalid SMILES skipped at row {index + 1}: {smi}"
            )
            continue

        canonical_smi = Chem.MolToSmiles(
            mol,
            canonical=True,
            isomericSmiles=True,
        )

        valid_smiles.append(canonical_smi)
        valid_indices.append(index)

    return valid_smiles, np.asarray(valid_indices, dtype=int), errors


def _build_dataset(smiles_list, num_workers=0):
    """
    Build a Chemprop dataset from already validated SMILES.
    """
    if not smiles_list:
        raise ValueError("No valid SMILES were provided to Chemprop.")

    datapoints = [
        data.MoleculeDatapoint.from_smi(smi)
        for smi in smiles_list
    ]

    featurizer = featurizers.SimpleMoleculeMolGraphFeaturizer()
    dataset = data.MoleculeDataset(datapoints, featurizer)

    loader = data.build_dataloader(
        dataset,
        shuffle=False,
        num_workers=num_workers,
    )

    return dataset, loader


def _predict_single_ckpt(ckpt_path, dataset, num_workers=0):
    device = "cuda" if torch.cuda.is_available() else "cpu"

    model: models.MPNN = models.MPNN.load_from_checkpoint(
        ckpt_path,
        map_location=device,
    ).eval()

    loader = data.build_dataloader(
        dataset,
        shuffle=False,
        num_workers=num_workers,
    )

    trainer = pl.Trainer(
        logger=False,
        enable_progress_bar=False,
        accelerator="gpu" if device == "cuda" else "cpu",
        devices=1,
        plugins=[LightningEnvironment()],
    )

    preds_batches = trainer.predict(model, loader)

    if not preds_batches:
        return np.asarray([], dtype=float)

    preds = np.concatenate(
        [batch.detach().cpu().numpy().reshape(-1) for batch in preds_batches]
    )

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
    """
    Run classification and optional regression models.

    Invalid SMILES are preserved in the output CSV. Their prediction values
    remain NA, valid_smiles is set to 0, and smiles_error describes the issue.
    """
    original_smiles = [str(smi).strip() for smi in smiles_list]
    n = len(original_smiles)

    df_log = pd.DataFrame(
        {
            "SMILES": original_smiles,
            "canonical_SMILES": pd.Series([pd.NA] * n, dtype="string"),
            "valid_smiles": np.zeros(n, dtype=int),
            "smiles_error": pd.Series([""] * n, dtype="string"),
        }
    )

    if cls_ckpt is None or cls_threshold is None:
        raise ValueError(
            "Classification checkpoint (--cls-ckpt) and "
            "threshold (--cls-threshold) are required."
        )

    valid_smiles, valid_indices, errors = _validate_smiles(original_smiles)
    df_log["smiles_error"] = errors

    if len(valid_indices) > 0:
        df_log.loc[valid_indices, "canonical_SMILES"] = valid_smiles
        df_log.loc[valid_indices, "valid_smiles"] = 1

    # Initialize all output columns so invalid rows remain in the CSV with NA.
    df_log["cls_prob"] = np.nan
    df_log["pass_cls"] = pd.Series([pd.NA] * n, dtype="Int64")

    if reg1_ckpt is not None and reg1_threshold is not None:
        df_log["reg1_pred"] = np.nan
        df_log["pass_reg1"] = pd.Series([pd.NA] * n, dtype="Int64")

    if reg2_ckpt is not None and reg2_threshold is not None:
        df_log["reg2_pred"] = np.nan
        df_log["pass_reg2"] = pd.Series([pd.NA] * n, dtype="Int64")

    df_log["final_pass"] = pd.Series([pd.NA] * n, dtype="Int64")

    if len(valid_indices) == 0:
        print("[WARNING] No valid SMILES were found. No predictions were run.")

        if log_csv_path:
            Path(os.path.dirname(log_csv_path) or ".").mkdir(
                parents=True,
                exist_ok=True,
            )
            df_log.to_csv(log_csv_path, index=False)
            print(f"[SAVED] Log to {log_csv_path}")

        if out_smi_path:
            Path(os.path.dirname(out_smi_path) or ".").mkdir(
                parents=True,
                exist_ok=True,
            )
            Path(out_smi_path).write_text("", encoding="utf-8")
            print(f"[SAVED] 0 SMILES to {out_smi_path}")

        return []

    # Stage 1: classification
    ds_all, _ = _build_dataset(
        valid_smiles,
        num_workers=num_workers,
    )

    cls_preds = _predict_single_ckpt(
        cls_ckpt,
        ds_all,
        num_workers=num_workers,
    ).reshape(-1)

    if len(cls_preds) != len(valid_indices):
        raise RuntimeError(
            "Classification prediction count does not match the number "
            "of valid input SMILES."
        )

    df_log.loc[valid_indices, "cls_prob"] = cls_preds

    pass_cls_local = cls_preds >= cls_threshold
    df_log.loc[valid_indices, "pass_cls"] = pass_cls_local.astype(int)

    idx_final = valid_indices[pass_cls_local]

    print(
        f"[OK] Classification pass: {len(idx_final)} / {len(valid_indices)} "
        f"valid molecules (threshold >= {cls_threshold})"
    )

    # Stage 2: regression model 1
    if (
        reg1_ckpt is not None
        and reg1_threshold is not None
        and len(idx_final) > 0
    ):
        smiles_r1 = [
            str(df_log.loc[index, "canonical_SMILES"])
            for index in idx_final
        ]

        ds_r1, _ = _build_dataset(
            smiles_r1,
            num_workers=num_workers,
        )

        reg1_preds = _predict_single_ckpt(
            reg1_ckpt,
            ds_r1,
            num_workers=num_workers,
        ).reshape(-1)

        if len(reg1_preds) != len(idx_final):
            raise RuntimeError(
                "Regression #1 prediction count does not match its input."
            )

        df_log.loc[idx_final, "reg1_pred"] = reg1_preds

        pass_r1_local = reg1_preds > reg1_threshold
        df_log.loc[idx_final, "pass_reg1"] = pass_r1_local.astype(int)

        idx_final = idx_final[pass_r1_local]

        print(
            f"[OK] Regression #1 pass: {len(idx_final)} molecules "
            f"(threshold > {reg1_threshold})"
        )

    # Stage 3: regression model 2
    if (
        reg2_ckpt is not None
        and reg2_threshold is not None
        and len(idx_final) > 0
    ):
        smiles_r2 = [
            str(df_log.loc[index, "canonical_SMILES"])
            for index in idx_final
        ]

        ds_r2, _ = _build_dataset(
            smiles_r2,
            num_workers=num_workers,
        )

        reg2_preds = _predict_single_ckpt(
            reg2_ckpt,
            ds_r2,
            num_workers=num_workers,
        ).reshape(-1)

        if len(reg2_preds) != len(idx_final):
            raise RuntimeError(
                "Regression #2 prediction count does not match its input."
            )

        df_log.loc[idx_final, "reg2_pred"] = reg2_preds

        pass_r2_local = reg2_preds > reg2_threshold
        df_log.loc[idx_final, "pass_reg2"] = pass_r2_local.astype(int)

        idx_final = idx_final[pass_r2_local]

        print(
            f"[OK] Regression #2 pass: {len(idx_final)} molecules "
            f"(threshold > {reg2_threshold})"
        )

    # Valid molecules that reached the final stage receive 0 or 1.
    df_log.loc[valid_indices, "final_pass"] = 0
    df_log.loc[idx_final, "final_pass"] = 1

    final_smiles = [
        str(df_log.loc[index, "canonical_SMILES"])
        for index in idx_final
    ]

    print(f"[FINAL] Selected: {len(final_smiles)} / {n}")
    print(f"[INFO] Invalid SMILES retained with NA values: {n - len(valid_indices)}")

    if log_csv_path:
        Path(os.path.dirname(log_csv_path) or ".").mkdir(
            parents=True,
            exist_ok=True,
        )
        df_log.to_csv(log_csv_path, index=False)
        print(f"[SAVED] Log to {log_csv_path}")

    if out_smi_path:
        Path(os.path.dirname(out_smi_path) or ".").mkdir(
            parents=True,
            exist_ok=True,
        )

        with open(out_smi_path, "w", encoding="utf-8") as handle:
            for smi in final_smiles:
                handle.write(smi + "\n")

        print(f"[SAVED] {len(final_smiles)} SMILES to {out_smi_path}")

    return final_smiles


def _read_smiles_file(path: str) -> List[str]:
    """
    Read the first whitespace-separated field from each non-empty line.

    This supports:
        SMILES
        SMILES ID
        tab-separated SMILES files

    A header whose first field is SMILES is skipped.
    """
    smiles = []

    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()

            if not stripped or stripped.startswith("#"):
                continue

            first_field = stripped.split()[0].strip()

            if line_number == 1 and first_field.lower() in {
                "smiles",
                "smile",
                "canonical_smiles",
                "canonicalsmiles",
            }:
                continue

            smiles.append(first_field)

    return smiles


if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description=(
            "Three-stage Chemprop filtering with invalid-SMILES handling."
        )
    )

    ap.add_argument(
        "--in-smi",
        required=True,
        help="Input SMILES file (.smi)",
    )
    ap.add_argument(
        "--out-csv",
        required=True,
        help="Output CSV log file",
    )
    ap.add_argument(
        "--out-smi",
        required=True,
        help="Output SMILES file containing passed molecules",
    )
    ap.add_argument(
        "--cls-ckpt",
        required=True,
        help="Classification checkpoint",
    )
    ap.add_argument(
        "--cls-threshold",
        type=float,
        required=True,
        help="Classification threshold (probability >= value)",
    )
    ap.add_argument(
        "--reg1-ckpt",
        default=None,
        help="Regression #1 checkpoint or None",
    )
    ap.add_argument(
        "--reg1-threshold",
        type=float,
        default=None,
        help="Regression #1 threshold (prediction > value)",
    )
    ap.add_argument(
        "--reg2-ckpt",
        default=None,
        help="Regression #2 checkpoint or None",
    )
    ap.add_argument(
        "--reg2-threshold",
        type=float,
        default=None,
        help="Regression #2 threshold (prediction > value)",
    )
    ap.add_argument(
        "--num-workers",
        type=int,
        default=4,
        help="Number of DataLoader workers",
    )

    args = ap.parse_args()

    smiles_list = _read_smiles_file(args.in_smi)

    if not smiles_list:
        raise ValueError(f"No SMILES were found in input file: {args.in_smi}")

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