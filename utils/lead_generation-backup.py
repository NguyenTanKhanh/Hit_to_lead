import argparse
import os
import subprocess
import pandas as pd
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem
from drugex.training.generators import GraphTransformer
from drugex.data.corpus.vocabulary import VocGraph
from drugex.training.scorers.properties import Property
from drugex.training.scorers.modifiers import ClippedScore
from drugex.training.environment import DrugExEnvironment
from drugex.training.rewards import WeightedSum

def replace_iodine_with_hydrogen(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")
    rw = Chem.RWMol(Chem.MolFromSmiles(Chem.MolToSmiles(mol)))
    for atom in rw.GetAtoms():
        if atom.GetAtomicNum() == 53:
            atom.SetAtomicNum(1)
    new_mol = rw.GetMol()
    Chem.SanitizeMol(new_mol)
    return Chem.MolToSmiles(new_mol)

def generate_molecules(model_path, frags, num_samples, output_dir, gpu_id=0, min_sa=0.9):
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    use_gpus = (gpu_id,) if gpu_id >= 0 else None

    sa = Property("SA", modifier=ClippedScore(lower_x=7, upper_x=3))
    env = DrugExEnvironment([sa], [1.0], reward_scheme=WeightedSum())

    model = GraphTransformer(voc_trg=VocGraph(), use_gpus=use_gpus)
    model.loadStatesFromFile(model_path)

    print("\n\U0001F680 Generating molecules...")
    raw = model.generate(input_frags=frags, num_samples=num_samples)

    if isinstance(raw, dict) and 'smiles' in raw:
        smiles = raw['smiles']
    elif isinstance(raw, pd.DataFrame):
        col = next((c for c in raw.columns if c.lower() == 'smiles'), None)
        smiles = raw[col] if col else []
    else:
        smiles = [Chem.MolToSmiles(m) for m in raw]

    smiles = [replace_iodine_with_hydrogen(s) for s in smiles]

    valid = []
    logs = []
    for s in smiles:
        mol = Chem.MolFromSmiles(s)
        if not mol:
            logs.append(f"❌ Invalid SMILES: {s}")
            continue
        score = sa([mol])[0]
        if score > min_sa:
            cano = Chem.MolToSmiles(mol, canonical=True)
            valid.append(cano)
            logs.append(f"✅ {cano} (SA={score:.3f})")
        else:
            logs.append(f"⏭ {s} (SA={score:.3f})")

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    log_path = Path(output_dir) / "sa_filter_log.txt"
    log_path.write_text("\n".join(logs))
    print("\n".join(logs))

    if not valid:
        print("\u26A0 No valid SMILES passed the SA filter.")
        return ""

    smi_path = Path(output_dir) / "generated_smiles.smi"
    pd.DataFrame({'SMILES': valid}).to_csv(smi_path, index=False, header=False)
    print(f"\n\U0001F4BE Saved {len(valid)} SMILES to: {smi_path}")
    return str(smi_path)

def predict_ensemble_conservative(smiles_list, ckpt_dir, output_path):
    if ckpt_dir is None or not os.path.isdir(ckpt_dir):
        print("⚠ No checkpoint directory provided or found, skipping prediction.")
        return None

    import torch
    from chemprop import data, featurizers, models
    from lightning import pytorch as pl
    import numpy as np

    test_datapoints = [data.MoleculeDatapoint.from_smi(smi) for smi in smiles_list]
    featurizer = featurizers.SimpleMoleculeMolGraphFeaturizer()
    test_dataset = data.MoleculeDataset(test_datapoints, featurizer)
    test_loader = data.build_dataloader(test_dataset, shuffle=False, num_workers=4)

    ckpt_paths = [os.path.join(ckpt_dir, f) for f in os.listdir(ckpt_dir) if f.endswith(".ckpt")]
    ckpt_paths.sort()

    all_preds = []
    device = "cuda" if torch.cuda.is_available() else "cpu"

    for ckpt in ckpt_paths:
        print(f"\U0001F50D Loading model: {ckpt}")
        model: models.MPNN = models.MPNN.load_from_checkpoint(ckpt, map_location=device).eval()

        trainer = pl.Trainer(
            logger=False,
            enable_progress_bar=False,
            accelerator="gpu" if device == "cuda" else "cpu",
            devices=1
        )

        pred_batches = trainer.predict(model, test_loader)
        probs = np.concatenate([b.cpu().numpy().flatten() for b in pred_batches])
        preds = (probs >= 0.5).astype(int)
        all_preds.append(preds)

    votes = np.stack(all_preds, axis=1)
    ensemble_pred = (votes.sum(axis=1) == len(all_preds)).astype(int)

    df_results = pd.DataFrame({"smile": smiles_list})
    for i, preds in enumerate(all_preds, start=1):
        df_results[f"model_{i}"] = preds
    df_results["ensemble_pred"] = ensemble_pred

    df_results.to_csv(output_path, index=False)
    print(f"✅ Predictions saved to {output_path}")

    active_smis = df_results[df_results["ensemble_pred"] == 1]["smile"].tolist()
    return active_smis

def convert_smiles_to_pdbqt(smiles_file, output_dir):
    pdbqt_dir = Path(output_dir) / "pdbqt_files"
    tmp_dir = Path(output_dir) / "temp_sdf"
    pdbqt_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    with open(smiles_file) as f:
        smiles_list = [l.strip() for l in f if l.strip()]

    success = 0
    for i, smi in enumerate(smiles_list):
        sdf = tmp_dir / f"mol_{i}.sdf"
        pdbqt = pdbqt_dir / f"mol_{i}.pdbqt"
        try:
            subprocess.run(['obabel', f'-:{smi}', '-O', str(sdf), '--gen3D'], check=True)
            subprocess.run(['obabel', str(sdf), '-O', str(pdbqt), '--minimize', '--ff', 'MMFF94', '--steps', '500'], check=True)
            if pdbqt.exists() and pdbqt.stat().st_size > 0:
                print(f"✔ Success: {smi[:20]}... -> {pdbqt}")
                success += 1
            else:
                print(f"⚠ Failed to create {pdbqt}")
        except Exception as e:
            print(f"❌ Error for {smi}: {e}")

    for f in tmp_dir.glob('*.sdf'):
        f.unlink()
    tmp_dir.rmdir()

    print(f"\n\U0001F389 PDBQT conversion complete. Success: {success}/{len(smiles_list)}")
    print(f"\U0001F4C2 Files in: {pdbqt_dir}")
    return str(pdbqt_dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="DrugEx Lead Generation Pipeline")
    parser.add_argument('--generation_model', required=True, help='Path to trained .pkg model used for molecule generation')
    parser.add_argument('--prediction_model', required=False, default=None, help='Path to directory of .ckpt files (optional, used for prediction)')
    parser.add_argument('--frags', nargs='+', required=True, help='List of scaffold fragments')
    parser.add_argument('--num_samples', type=int, default=100, help='Number of molecules to generate')
    parser.add_argument('--output_dir', default='lead_generation_output', help='Main output directory')
    parser.add_argument('--gpu', type=int, default=0, help='GPU ID (or -1 for CPU)')
    parser.add_argument('--min_sa', type=float, default=0.9, help='Minimum normalized SA score (0-1)')
    args = parser.parse_args()

    smiles_file = generate_molecules(
        model_path=args.generation_model,
        frags=args.frags,
        num_samples=args.num_samples,
        output_dir=args.output_dir,
        gpu_id=args.gpu,
        min_sa=args.min_sa
    )

    if smiles_file and args.prediction_model:
        with open(smiles_file) as f:
            smiles = [line.strip() for line in f if line.strip()]

        pred_path = Path(args.output_dir) / "predicted_activity.csv"
        active_smiles = predict_ensemble_conservative(smiles, ckpt_dir=args.prediction_model, output_path=pred_path)

        if active_smiles:
            active_file = Path(args.output_dir) / "active_smiles.smi"
            pd.DataFrame({'SMILES': active_smiles}).to_csv(active_file, index=False, header=False)
            convert_smiles_to_pdbqt(active_file, output_dir=args.output_dir)
    elif smiles_file:
        convert_smiles_to_pdbqt(smiles_file, output_dir=args.output_dir)

