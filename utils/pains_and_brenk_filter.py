#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
from rdkit import Chem
from rdkit.Chem import FilterCatalog


def build_catalog():
    """Tạo catalog gồm PAINS_A/B/C + BRENK"""
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A)
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B)
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C)
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
    return FilterCatalog.FilterCatalog(params)


def filter_pains_brenk(df, smiles_col="SMILES"):
    """Lọc compound theo PAINS + BRENK, chỉ trả về list pass"""
    catalog = build_catalog()
    kept_rows, removed_count = [], 0

    for _, row in df.iterrows():
        smi = str(row[smiles_col]).strip()
        mol = Chem.MolFromSmiles(smi)
        if not mol:
            removed_count += 1
            continue

        entries = catalog.GetMatches(mol)
        if entries:
            removed_count += 1
        else:
            kept_rows.append(row.to_dict())

    df_kept = pd.DataFrame(kept_rows).reset_index(drop=True)
    return df_kept, removed_count


def main():
    parser = argparse.ArgumentParser(
        description="Filter compounds using PAINS + BRENK (RDKit)"
    )
    parser.add_argument("-i", "--input", required=True, help="Input .smi (tab-separated)")
    parser.add_argument("-o", "--output", required=True, help="Output .smi (pass compounds)")
    parser.add_argument("--smiles-col", default="SMILES", help="Column containing SMILES (default: SMILES)")
    args = parser.parse_args()

    # Load input file
    try:
        df = pd.read_csv(args.input, sep="\t")
        # Nếu file không có header thì thêm header thủ công
        if args.smiles_col not in df.columns:
            df = pd.read_csv(args.input, sep="\t", header=None, names=[args.smiles_col])
    except Exception:
        # fallback nếu file quá đơn giản (chỉ 1 cột không header)
        df = pd.read_csv(args.input, sep="\t", header=None, names=[args.smiles_col])

    print(f"📥 Loaded {len(df)} compounds from {args.input}")

    # Filter
    df_kept, removed_count = filter_pains_brenk(df, smiles_col=args.smiles_col)
    print(f"✅ Pass: {len(df_kept)}")
    print(f"❌ Removed: {removed_count}")

    # Save only pass
    df_kept.to_csv(args.output, sep="\t", index=False, header=True)
    print("💾 Saved pass compounds:", args.output)


if __name__ == "__main__":
    main()
