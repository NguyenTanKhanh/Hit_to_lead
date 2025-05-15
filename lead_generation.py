# -*- coding: utf-8 -*-
import os
import subprocess
import argparse
import pandas as pd
from pathlib import Path
from rdkit import Chem
from drugex.training.generators import GraphTransformer
from drugex.data.corpus.vocabulary import VocGraph
from drugex.training.scorers.properties import Property
from drugex.training.scorers.modifiers import ClippedScore
from drugex.training.environment import DrugExEnvironment
from drugex.training.rewards import WeightedSum

def generate_molecules(model_path, frags, num_samples, output_dir, gpu_id=0, min_sa=0.9):
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    use_gpus = (gpu_id,) if gpu_id >= 0 else None

    sascore = Property("SA")
    sa_clipped = ClippedScore(lower_x=7, upper_x=3)  # Normalize to 0–1
    qed = Property("QED")
    environment = DrugExEnvironment([sascore, qed], [0.5, 0.5], reward_scheme=WeightedSum())

    model = GraphTransformer(voc_trg=VocGraph(), use_gpus=use_gpus)
    model.loadStatesFromFile(model_path)

    print("\n?? Generating molecules...")
    results = model.generate(input_frags=frags, num_samples=num_samples, evaluator=environment)

    if isinstance(results, dict) and 'smiles' in results:
        smiles = results['smiles']
    elif isinstance(results, pd.DataFrame):
        smiles_col = next((col for col in results.columns if col.lower() == 'smiles'), None)
        smiles = results[smiles_col] if smiles_col else []
    else:
        smiles = [Chem.MolToSmiles(mol) for mol in results]

    valid_smiles = []
    log_lines = []
    for s in smiles:
        mol = Chem.MolFromSmiles(s)
        if mol:
            sa_value = sascore([mol])[0]         # float
            sa_score = sa_clipped(sa_value)      # normalized float (0–1)
            if sa_score > min_sa:
                canonical = Chem.MolToSmiles(mol, canonical=True)
                valid_smiles.append(canonical)
                log_lines.append(f"? Accepted molecule (SA={sa_score:.3f}): {canonical}")
            else:
                log_lines.append(f"? Skipped (SA={sa_score:.3f}): {s}")

    log_path = os.path.join(output_dir, "sa_filter_log.txt")
    with open(log_path, "w") as f:
        f.write("\n".join(log_lines))
    for line in log_lines:
        print(line)

    if not valid_smiles:
        print("? No valid SMILES passed the SA filter.")
        return ""

    smiles_path = os.path.join(output_dir, "generated_smiles.smi")
    pd.DataFrame({'SMILES': valid_smiles}).to_csv(smiles_path, index=False, header=False)
    print(f"\n? Saved {len(valid_smiles)} SMILES to: {smiles_path}")
    return smiles_path

def convert_smiles_to_pdbqt(smiles_file, output_dir):
    pdbqt_dir = os.path.join(output_dir, "pdbqt_files")
    temp_dir = os.path.join(output_dir, "temp_sdf")
    Path(pdbqt_dir).mkdir(parents=True, exist_ok=True)
    Path(temp_dir).mkdir(parents=True, exist_ok=True)

    with open(smiles_file) as f:
        smiles_list = [line.strip() for line in f if line.strip()]

    success_count = 0
    for i, smile in enumerate(smiles_list):
        sdf_path = os.path.join(temp_dir, f"mol_{i}.sdf")
        pdbqt_path = os.path.join(pdbqt_dir, f"mol_{i}.pdbqt")

        try:
            subprocess.run(['obabel', f'-:{smile}', '-O', sdf_path, '--gen3D'], check=True)
            subprocess.run(['obabel', sdf_path, '-O', pdbqt_path, '--minimize', '--ff', 'MMFF94', '--steps', '500'], check=True)

            if os.path.exists(pdbqt_path) and os.path.getsize(pdbqt_path) > 0:
                print(f"? Success: {smile[:25]}... -> {pdbqt_path}")
                success_count += 1
            else:
                print(f"?? Failed to create: {pdbqt_path}")
        except subprocess.CalledProcessError as e:
            print(f"? Open Babel error: {smile[:25]} - {e}")
        except Exception as e:
            print(f"? Unexpected error: {smile[:25]} - {e}")

    for f in Path(temp_dir).glob("*.sdf"):
        f.unlink()
    Path(temp_dir).rmdir()

    print(f"\n? PDBQT conversion complete. Success: {success_count}/{len(smiles_list)}")
    print(f"?? PDBQT files in: {os.path.abspath(pdbqt_dir)}")
    return pdbqt_dir

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="DrugEx Lead Generation Pipeline")
    parser.add_argument('--pkg', required=True, help='Path to trained .pkg model')
    parser.add_argument('--frags', nargs='+', required=True, help='List of scaffold fragments')
    parser.add_argument('--num_samples', type=int, default=100, help='Number of molecules to generate')
    parser.add_argument('--output_dir', default='lead_generation_output', help='Main output directory')
    parser.add_argument('--gpu', type=int, default=0, help='GPU ID (or -1 for CPU)')
    parser.add_argument('--min_sa', type=float, default=0.9, help='Minimum normalized SA score (0-1)')
    args = parser.parse_args()

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    smiles_file = generate_molecules(
        args.pkg,
        args.frags,
        args.num_samples,
        output_dir=args.output_dir,
        gpu_id=args.gpu,
        min_sa=args.min_sa
    )

    if smiles_file:
        convert_smiles_to_pdbqt(smiles_file, output_dir=args.output_dir)
