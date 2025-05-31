#!/usr/bin/env python3

import argparse
import io
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image

def attach_iodine(smiles, heavy_indices):
    """
    Attach iodine atom ([I]) to each specified heavy-atom index in a molecule.
    Returns the new SMILES and the modified RDKit Mol object.
    """
    base_mol = Chem.MolFromSmiles(smiles)
    if base_mol is None:
        raise ValueError("Invalid SMILES string")
    base_mol = Chem.MolFromSmiles(Chem.MolToSmiles(base_mol))  # canonicalize

    rw = Chem.RWMol(base_mol)
    max_idx = rw.GetNumAtoms() - 1

    for hi in heavy_indices:
        if not (0 <= hi <= max_idx):
            raise IndexError(f"Atom index {hi} out of range (0–{max_idx})")

        iodine_idx = rw.AddAtom(Chem.Atom(53))  # Atomic number of Iodine
        rw.AddBond(hi, iodine_idx, Chem.BondType.SINGLE)

    new_mol = rw.GetMol()
    Chem.SanitizeMol(new_mol)
    return Chem.MolToSmiles(new_mol), new_mol

def parse_args():
    p = argparse.ArgumentParser(
        description="Attach iodine ([I]) to specified atom indices, save SMILES and image."
    )
    p.add_argument("-s", "--smiles", required=True, help="Input SMILES string")
    p.add_argument("-i", "--indices", required=True,
                   help="Comma-separated list of 0-based heavy-atom indices, e.g. 5,6,21")
    p.add_argument("--output-smi", default="molecule_with_I.smi", help="Output .smi filename")
    p.add_argument("--output-img", default="molecule_with_I.png", help="Output PNG filename")
    p.add_argument("-W", "--width", type=int, default=400, help="Image width in pixels")
    p.add_argument("-H", "--height", type=int, default=400, help="Image height in pixels")
    return p.parse_args()

def main():
    args = parse_args()
    try:
        idxs = [int(x) for x in args.indices.split(",") if x.strip()]
    except ValueError:
        print("Error: --indices must be comma-separated integers", file=sys.stderr)
        sys.exit(1)

    try:
        new_smiles, new_mol = attach_iodine(args.smiles, idxs)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    # Write SMILES to file
    with open(args.output_smi, "w") as f:
        f.write(new_smiles + "\n")
    print(f"✔ New SMILES written to {args.output_smi}")

    # Draw molecule and save image
    AllChem.Compute2DCoords(new_mol)
    drawer = rdMolDraw2D.MolDraw2DCairo(args.width, args.height)
    drawer.DrawMolecule(new_mol)
    drawer.FinishDrawing()
    png_bytes = drawer.GetDrawingText()
    img = Image.open(io.BytesIO(png_bytes))
    img.save(args.output_img)
    print(f"✔ Image saved to {args.output_img}")

if __name__ == "__main__":
    main()


