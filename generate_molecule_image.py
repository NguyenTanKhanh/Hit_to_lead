import os
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from pathlib import Path

def calculate_molecular_weight(mol):
    """Calculate molecular weight including hydrogens (implicit H’s counted)"""
    mol_with_h = Chem.AddHs(mol)
    return AllChem.CalcExactMolWt(mol_with_h)

def draw_molecule_with_scaffold(smiles, scaffold_smiles=None, mol_id=0, size=(400, 400)):
    """Return a PIL Image of the molecule, optionally highlighting the scaffold"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Could not parse SMILES: {smiles}")
        return None

    # Molecular weight with implicit H’s
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
                # map scaffold bonds onto the parent molecule
                for b in scaffold.GetBonds():
                    i1 = highlight_atoms[b.GetBeginAtomIdx()]
                    i2 = highlight_atoms[b.GetEndAtomIdx()]
                    bid = mol.GetBondBetweenAtoms(i1, i2).GetIdx()
                    highlight_bonds.append(bid)
                # light grey
                atom_colors = {i: (0.8, 0.8, 0.8) for i in highlight_atoms}
                bond_colors = {b: (0.8, 0.8, 0.8) for b in highlight_bonds}
            else:
                legend = f"Scaffold not found | {legend}"
        else:
            legend = f"Invalid scaffold | {legend}"

    # Draw via RDKit's high-level MolToImage
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
    """Read SMILES file (one per line) and save PNGs with optional scaffold highlighting"""
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
    # Path to your SMILES (one per line)
    generate_molecule_images(
        smiles_file="/home/khanhnt/Desktop/Hit_to_lead/lead_generation_output/generated_smiles.smi",
        output_dir="image_lead_molecule",
        scaffold_smiles="c1cnccn1"  # or None
    )


