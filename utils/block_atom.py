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
    Attach iodine atoms ([I]) to each specified heavy-atom index in a molecule.
    For each specified atom, all attached hydrogens are replaced by iodines.
    Example: CH2 -> CI2, CH3 -> CI3 at that position.
    Returns the new SMILES and the modified RDKit Mol object.
    """
    base_mol = Chem.MolFromSmiles(smiles)
    if base_mol is None:
        raise ValueError("Invalid SMILES string")

    # Keep canonicalization as in the original code
    base_mol = Chem.MolFromSmiles(Chem.MolToSmiles(base_mol))  # canonicalize

    # Number of heavy atoms in the original molecule (without explicit H)
    n_heavy = base_mol.GetNumAtoms()

    # Add explicit hydrogens so we can see and replace them
    mol_with_H = Chem.AddHs(base_mol)
    rw = Chem.RWMol(mol_with_H)

    # Maximum valid heavy-atom index
    max_idx = n_heavy - 1

    for hi in heavy_indices:
        if not (0 <= hi <= max_idx):
            raise IndexError(f"Atom index {hi} out of range (0-{max_idx})")

        # Heavy atom corresponding to index hi
        # After AddHs, hydrogens are appended at the end, so heavy-atom indices stay the same
        atom = rw.GetAtomWithIdx(hi)

        # Collect all neighboring hydrogens of this atom
        h_neighbors = [
            nbr.GetIdx()
            for nbr in atom.GetNeighbors()
            if nbr.GetAtomicNum() == 1  # Hydrogen
        ]

        # If there are no hydrogens attached, skip
        if not h_neighbors:
            continue

        # Add the same number of iodines as hydrogens
        for _ in h_neighbors:
            iodine_idx = rw.AddAtom(Chem.Atom(53))  # Atomic number of iodine
            rw.AddBond(hi, iodine_idx, Chem.BondType.SINGLE)

        # Remove all hydrogens that were replaced by iodines
        for h_idx in sorted(h_neighbors, reverse=True):
            rw.RemoveAtom(h_idx)

    new_mol = rw.GetMol()
    Chem.SanitizeMol(new_mol)

    # Remove explicit hydrogens to get a clean SMILES
    new_mol = Chem.RemoveHs(new_mol)

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


