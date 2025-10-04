import subprocess
from pathlib import Path

def convert_smiles_to_pdbqt(smiles_file, output_dir):
    """
    Convert a list of SMILES into PDBQT files using Open Babel.
    - smiles_file: đường dẫn file .smi (mỗi dòng = 1 SMILES)
    - output_dir: thư mục output
    """
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
            # Tạo cấu trúc 3D
            subprocess.run(['obabel', f'-:{smi}', '-O', str(sdf), '--gen3D'], check=True)
            # Tối ưu hóa và chuyển sang PDBQT
            subprocess.run(['obabel', str(sdf), '-O', str(pdbqt), '--minimize', '--ff', 'MMFF94', '--steps', '500'], check=True)
            if pdbqt.exists() and pdbqt.stat().st_size > 0:
                print(f"✔ Success: {smi[:30]}... -> {pdbqt}")
                success += 1
            else:
                print(f"⚠ Failed to create {pdbqt}")
        except Exception as e:
            print(f"❌ Error for {smi}: {e}")

    # Xóa file trung gian
    for f in tmp_dir.glob('*.sdf'):
        f.unlink()
    tmp_dir.rmdir()

    print(f"\n🎉 PDBQT conversion complete. Success: {success}/{len(smiles_list)}")
    print(f"📂 Files saved in: {pdbqt_dir}")
    return str(pdbqt_dir)
