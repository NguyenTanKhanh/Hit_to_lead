#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compound properties filter — Open Babel (protonation) + RDKit (descriptors)

Functions:
  1) protonate_smiles(smiles, ph) -> str
  2) rdkit_features(smiles) -> dict
  3) filter_smiles(smiles_list, ph, ranges) -> list[str]

CLI:
  # Protonate toàn bộ file .smi ở pH nhất định
  python utils/compound_properties_filter.py --mode protonate --in input.smi --out pH6_3.smi --ph 6.3

  # Tính features (RDKit) cho file .smi (Không protonate trong mode này)
  python utils/compound_properties_filter.py --mode features --in some.smi --out feats.csv

  # Lọc theo tiêu chí tại pH (protonate -> compute -> filter)
  python utils/compound_properties_filter.py --mode filter --in input.smi --out pass.smi --ph 7.4 \
      --mw 300 480 --logp 2 4 --tpsa 50 90 --hbd 0 1 --hba 3 6 --rotb 0 7 --logs -6 0
"""

import argparse
import csv
import shutil
import subprocess
from typing import List, Dict, Optional, Tuple

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors


# =========================
# Function 1: Protonation (Open Babel only)
# =========================

def protonate_smiles(smiles: str, ph: float) -> str:
    """
    Protonate a SMILES at given pH using Open Babel.
    Returns ONE protomer SMILES (string). Raises RuntimeError if obabel is missing or fails.
    """
    if shutil.which("obabel") is None:
        raise RuntimeError("Open Babel (obabel) not found in PATH.")
    cmd = ["obabel", f"-:{smiles}", "-p", str(ph), "-h", "-osmi"]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"obabel failed: {proc.stderr.strip()}")
    lines = [ln.strip() for ln in proc.stdout.splitlines() if ln.strip()]
    if not lines:
        raise RuntimeError("obabel produced no output.")
    return lines[0].split()[0]


# =========================
# Function 2: Features (RDKit)
# =========================

def calc_logS_esol(m: Chem.Mol) -> float:
    """
    Predict aqueous solubility (logS, mol/L) using Delaney ESOL model.
    Ref: Delaney JS. J Chem Inf Comput Sci. 2004.
    """
    molwt = Descriptors.MolWt(m)
    logp = Crippen.MolLogP(m)
    rotb = rdMolDescriptors.CalcNumRotatableBonds(m)
    hbd = rdMolDescriptors.CalcNumHBD(m)
    # ESOL linear model
    logS = 0.16 - 0.63 * logp - 0.0062 * molwt + 0.066 * rotb - 0.74 * hbd
    return logS


def rdkit_features(smiles: str) -> Optional[Dict[str, float]]:
    """Compute MW, LogP, TPSA, HBD, HBA, RotBonds, LogS from a SMILES."""
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    return {
        "SMILES": smiles,
        "MW": Descriptors.MolWt(m),
        "LogP": Crippen.MolLogP(m),
        "TPSA": rdMolDescriptors.CalcTPSA(m),
        "HBA": rdMolDescriptors.CalcNumHBA(m),
        "HBD": rdMolDescriptors.CalcNumHBD(m),
        "RotBonds": rdMolDescriptors.CalcNumRotatableBonds(m),
        "LogS": calc_logS_esol(m),
    }


# =========================
# Function 3: Filtering
# =========================

Ranges = Dict[str, Tuple[float, float]]

def filter_smiles(smiles_list: List[str], ph: float, ranges: Ranges) -> List[str]:
    """
    Protonate each SMILES at pH, compute features, filter by ranges; return passed SMILES (protomer at that pH).
    """
    def in_range(v: float, lo: float, hi: float) -> bool:
        return lo <= float(v) <= hi

    passed: List[str] = []
    for smi in smiles_list:
        try:
            prot = protonate_smiles(smi, ph=ph)
        except Exception:
            continue
        feats = rdkit_features(prot)
        if not feats:
            continue
        if all([
            in_range(feats["MW"],   *ranges["MW"]),
            in_range(feats["LogP"], *ranges["LogP"]),
            in_range(feats["TPSA"], *ranges["TPSA"]),
            in_range(feats["HBD"],  *ranges["HBD"]),
            in_range(feats["HBA"],  *ranges["HBA"]),
            in_range(feats["RotBonds"], *ranges["RotBonds"]),
            in_range(feats["LogS"], *ranges["LogS"]),
        ]):
            passed.append(prot)

    # de-dup while preserving order
    seen = set()
    out: List[str] = []
    for s in passed:
        if s not in seen:
            seen.add(s)
            out.append(s)
    return out


# =========================
# File helpers + CLI
# =========================

def _read_smi(path: str) -> List[str]:
    arr = []
    with open(path, "r", encoding="utf-8") as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            arr.append(ln.split()[0])
    return arr


def _write_smi(path: str, smiles: List[str]) -> None:
    with open(path, "w", encoding="utf-8") as f:
        for s in smiles:
            f.write(s + "\n")


def _write_csv(path: str, rows: List[Dict[str, float]]) -> None:
    headers = ["SMILES", "MW", "LogP", "TPSA", "HBA", "HBD", "RotBonds", "LogS"]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=headers)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in headers})


def main():
    ap = argparse.ArgumentParser(description="OpenBabel protonation + RDKit features + filtering (no fallbacks).")
    ap.add_argument("--mode", choices=["protonate", "features", "filter"], required=True)
    ap.add_argument("--in", dest="infile", required=True, help="Input .smi (one SMILES per line)")
    ap.add_argument("--out", dest="outfile", required=True, help="Output path (.smi for protonate/filter; .csv for features)")
    ap.add_argument("--ph", type=float, default=7.4, help="pH for protonation (protonate/filter)")
    # Ranges for filter mode
    ap.add_argument("--mw",   nargs=2, type=float, metavar=("LO","HI"))
    ap.add_argument("--logp", nargs=2, type=float, metavar=("LO","HI"))
    ap.add_argument("--tpsa", nargs=2, type=float, metavar=("LO","HI"))
    ap.add_argument("--hbd",  nargs=2, type=float, metavar=("LO","HI"))
    ap.add_argument("--hba",  nargs=2, type=float, metavar=("LO","HI"))
    ap.add_argument("--rotb", nargs=2, type=float, metavar=("LO","HI"))
    ap.add_argument("--logs", nargs=2, type=float, metavar=("LO","HI"))
    args = ap.parse_args()

    smiles = _read_smi(args.infile)

    if args.mode == "protonate":
        out_smis = []
        for s in smiles:
            try:
                out_smis.append(protonate_smiles(s, ph=args.ph))
            except Exception:
                continue
        _write_smi(args.outfile, out_smis)

    elif args.mode == "features":
        rows = []
        for s in smiles:
            f = rdkit_features(s)
            if f:
                rows.append(f)
        _write_csv(args.outfile, rows)

    elif args.mode == "filter":
        if not (args.mw and args.logp and args.tpsa and args.hbd and args.hba and args.rotb and args.logs):
            raise SystemExit("For --mode filter you must provide --mw --logp --tpsa --hbd --hba --rotb --logs ranges.")
        ranges: Ranges = {
            "MW": tuple(args.mw),
            "LogP": tuple(args.logp),
            "TPSA": tuple(args.tpsa),
            "HBD": tuple(args.hbd),
            "HBA": tuple(args.hba),
            "RotBonds": tuple(args.rotb),
            "LogS": tuple(args.logs),
        }  # type: ignore
        passed = filter_smiles(smiles, ph=args.ph, ranges=ranges)
        _write_smi(args.outfile, passed)


if __name__ == "__main__":
    main()
