import os
import sys
import argparse
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from pathlib import Path

def calculate_molecular_weight(mol):
    mol_with_h = Chem.AddHs(mol)
    return AllChem.CalcExactMolWt(mol_with_h)

def draw_molecule_with_scaffold(smiles, scaffold_smiles=None, mol_id=0, size=(400, 400)):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Could not parse SMILES: {smiles}")
        return None

    mw = calculate_molecular_weight(mol)
    legend = f"Mol {mol_id} | MW: {mw:.2f}"

    highlight_atoms = []
    highlight_bonds = []
    atom_colors = {}
    bond_colors = {}

    if scaffold_smiles:
        scaffold = Chem.MolFromSmiles(scaffold_smiles)
        if scaffold:
            matches = mol.GetSubstructMatches(scaffold)
            if matches:
                highlight_atoms = list(matches[0])
                for b in scaffold.GetBonds():
                    i1 = highlight_atoms[b.GetBeginAtomIdx()]
                    i2 = highlight_atoms[b.GetEndAtomIdx()]
                    bid = mol.GetBondBetweenAtoms(i1, i2).GetIdx()
                    highlight_bonds.append(bid)
                atom_colors = {i: (0.8, 0.8, 0.8) for i in highlight_atoms}
                bond_colors = {b: (0.8, 0.8, 0.8) for b in highlight_bonds}
            else:
                legend = f"Scaffold not found | {legend}"
        else:
            legend = f"Invalid scaffold | {legend}"

    img = Draw.MolToImage(
        mol,
        size=size,
        highlightAtoms=highlight_atoms or None,
        highlightBonds=highlight_bonds or None,
        highlightAtomColors=atom_colors or None,
        highlightBondColors=bond_colors or None,
        legend=legend
    )
    return img

def generate_molecule_images(smiles_file, output_dir="image_lead_molecule", scaffold_smiles=None):
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    with open(smiles_file) as f:
        smiles_list = [line.strip() for line in f if line.strip()]

    for idx, smi in enumerate(smiles_list):
        try:
            img = draw_molecule_with_scaffold(smi, scaffold_smiles, mol_id=idx)
            if img:
                out_path = os.path.join(output_dir, f"molecule_{idx}.png")
                img.save(out_path)
                print(f"Saved molecule {idx} to {out_path}")
        except Exception as e:
            print(f"Error processing molecule {idx}: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("smiles_file", help="Path to SMILES file (1 SMILES per line)")
    parser.add_argument("output_dir", help="Directory to save molecule images")
    parser.add_argument("scaffold_smiles", nargs="?", default=None, help="Optional scaffold SMILES to highlight")
    args = parser.parse_args()

    generate_molecule_images(args.smiles_file, args.output_dir, args.scaffold_smiles)

