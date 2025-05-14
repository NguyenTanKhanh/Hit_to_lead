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

def generate_molecules(model_path, frags, num_samples, output_dir, gpu_id=0):
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    use_gpus = (gpu_id,) if gpu_id >= 0 else None

    # Set up scoring environment
    sascore = Property("SA", modifier=ClippedScore(lower_x=7, upper_x=3))
    qed = Property("QED", modifier=ClippedScore(lower_x=0.2, upper_x=0.8))
    environment = DrugExEnvironment([sascore, qed], [0.5, 0.5], reward_scheme=WeightedSum())

    model = GraphTransformer(voc_trg=VocGraph(), use_gpus=use_gpus)
    model.loadStatesFromFile(model_path)

    print("Generating molecules...")
    results = model.generate(input_frags=frags, num_samples=num_samples, evaluator=environment)
    
    # Handle different output formats
    if isinstance(results, dict) and 'smiles' in results:
        smiles = results['smiles']
    elif isinstance(results, pd.DataFrame):
        smiles_col = next((col for col in results.columns if col.lower() == 'smiles'), None)
        if smiles_col:
            smiles = results[smiles_col]
        else:
            raise KeyError("No 'smiles' column found in DataFrame")
    else:
        try:
            smiles = [Chem.MolToSmiles(mol) for mol in results]
        except Exception as e:
            raise ValueError(f"Couldn't convert results to SMILES: {str(e)}")

    # Save SMILES to output directory
    smiles_path = os.path.join(output_dir, "generated_smiles.smi")
    df = pd.DataFrame({'SMILES': smiles})
    df.to_csv(smiles_path, index=False, header=False)
    print(f"Generated {len(df)} molecules saved to '{smiles_path}'")
    return smiles_path

def convert_smiles_to_pdbqt(smiles_file, output_dir):
    """
    Convert SMILES to PDBQT using OpenBabel
    """
    pdbqt_dir = os.path.join(output_dir, "pdbqt_files")
    temp_dir = os.path.join(output_dir, "temp_obabel")
    
    # Create directories
    Path(pdbqt_dir).mkdir(parents=True, exist_ok=True)
    Path(temp_dir).mkdir(exist_ok=True)
    
    # Read SMILES
    with open(smiles_file) as f:
        smiles_list = [line.strip() for line in f if line.strip()]
    
    success_count = 0
    
    for i, smile in enumerate(smiles_list):
        try:
            sdf_path = os.path.join(temp_dir, f"mol_{i}.sdf")
            pdbqt_path = os.path.join(pdbqt_dir, f"mol_{i}.pdbqt")
            
            # Convert SMILES to 3D SDF
            obabel_cmd = f'obabel -:"{smile}" -O "{sdf_path}" --gen3D 2>/dev/null'
            subprocess.run(obabel_cmd, shell=True, check=True)
            
            # Convert SDF to PDBQT with minimization
            obabel_cmd = f'obabel "{sdf_path}" -O "{pdbqt_path}" --minimize --ff MMFF94 --steps 500 2>/dev/null'
            subprocess.run(obabel_cmd, shell=True, check=True)
            
            if os.path.exists(pdbqt_path) and os.path.getsize(pdbqt_path) > 0:
                success_count += 1
                print(f"✅ Success: {smile[:20]}... -> pdbqt_files/mol_{i}.pdbqt")
                
        except subprocess.CalledProcessError:
            print(f"❌ Failed: {smile[:20]}... (invalid SMILES or conversion error)")
        except Exception as e:
            print(f"❌ Unexpected error with {smile[:20]}...: {str(e)}")

    # Cleanup
    for f in Path(temp_dir).glob("*"):
        f.unlink()
    Path(temp_dir).rmdir()
    
    print(f"\n🎯 Conversion complete. Success rate: {success_count}/{len(smiles_list)}")
    print(f"All files saved in: {os.path.abspath(output_dir)}")
    return pdbqt_dir

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="DrugEx Lead Generation Pipeline")
    parser.add_argument('--pkg', required=True, help='Path to trained .pkg model')
    parser.add_argument('--frags', nargs='+', required=True, help='List of scaffold fragments')
    parser.add_argument('--num_samples', type=int, default=100, help='Number of molecules to generate')
    parser.add_argument('--output_dir', default='lead_generation_output', help='Main output directory')
    parser.add_argument('--gpu', type=int, default=0, help='GPU ID (or -1 for CPU)')
    args = parser.parse_args()

    # Create main output directory
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    
    # Run pipeline
    smiles_file = generate_molecules(
        args.pkg, 
        args.frags, 
        args.num_samples, 
        output_dir=args.output_dir,
        gpu_id=args.gpu
    )
    
    pdbqt_dir = convert_smiles_to_pdbqt(
        smiles_file,
        output_dir=args.output_dir
    )