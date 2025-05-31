#!/usr/bin/env python3
"""
draw_smiles.py

Draw a SMILES string with heavy‐atom indices, save to PNG.
"""

import argparse
import io
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image
import sys

def show_and_save_smiles_with_indices(smiles, img_size, filename):
    """
    Draw a SMILES with each heavy atom’s 0-based index, and save the PNG to `filename`.
    """
    # 1. Parse & canonicalize
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles!r}")
    mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))

    # 2. Generate 2D coords
    AllChem.Compute2DCoords(mol)

    # 3. Set up RDKit’s Cairo drawer
    width, height = img_size
    drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
    opts = drawer.drawOptions()
    opts.addAtomIndices = False

    # 4. Add index notes to heavy atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() > 1:
            atom.SetProp("atomNote", str(atom.GetIdx()))

    # 5. Draw molecule
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    png_bytes = drawer.GetDrawingText()

    # 6. Ensure directory exists
    os.makedirs(os.path.dirname(filename), exist_ok=True)

    # 7. Save to file
    with open(filename, "wb") as f:
        f.write(png_bytes)
    print(f"✔ Image saved as {filename}")

def parse_args():
    p = argparse.ArgumentParser(
        description="Draw a SMILES string with heavy‐atom indices and save as PNG."
    )
    p.add_argument("-s", "--smiles", required=True, help="SMILES string to draw")
    p.add_argument("-o", "--output", default="molecule_indices.png", help="Output PNG filename")
    p.add_argument("-W", "--width", type=int, default=400, help="Image width in pixels")
    p.add_argument("-H", "--height", type=int, default=400, help="Image height in pixels")
    return p.parse_args()

def main():
    args = parse_args()
    try:
        show_and_save_smiles_with_indices(
            args.smiles,
            img_size=(args.width, args.height),
            filename=args.output
        )
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
