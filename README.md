# Hit-to-Lead Optimization Pipeline

This project integrates **DrugEx** and **AutoDock Vina** into a streamlined pipeline for hit-to-lead optimization. It supports scaffold-based molecule generation with fragment attachment, and virtual screening via molecular docking.

---

## 🔧 Installation

Make sure you have [Conda](https://docs.conda.io/en/latest/) installed. Then follow these steps:

```bash
git clone https://github.com/NguyenTanKhanh/Hit_to_lead.git
cd Hit_to_lead

# Create and activate a new Conda environment
conda create -n Hit_to_lead python=3.11
conda activate Hit_to_lead

# Install DrugEx and ChemProp
pip install git+https://github.com/CDDLeiden/DrugEx.git@master
pip install chemprop

# Install cheminformatics and scientific libraries
conda install -c conda-forge rdkit pandas numpy matplotlib scikit-learn
```

---

## 🚀 Pipeline Overview

### 1. Fragment Preparation

Generate dummy molecules from fragment combinations using DrugEx.

```bash
bash run_training.sh
```

---

### 2. Scaffold Preparation

Prepare a scaffold in SMILES format, visualize indexed atoms, and define attachment constraints.

```bash
bash run_scaffold_molecule_index.sh
```

After identifying atoms to block, apply the modification:

```bash
bash run_block_atom.sh
```


---

### 3. Molecule Generation from Scaffold

Use the modified scaffold to grow lead-like molecules:

```bash
bash run_lead_generation.sh
```


---

### 4. Molecular Docking with AutoDock Vina

Dock generated molecules against your target protein:

```bash
bash run_vinascreen.sh
```

---

### 🧪 Optional: Molecule Visualization

To visualize generated molecules as 2D images:

```bash
bash run_generate_molecule_image.sh
```

---

## 📚 References

- [DrugEx – CDD Leiden](https://github.com/CDDLeiden/DrugEx)
- [VinaScreen – Yassir Boulaamane](https://github.com/yboulaamane/VinaScreen/tree/main)

---

## ⚖️ License

**This project has no formal license. It is provided freely for educational, academic, and non-commercial use.**  
Feel free to use, modify, and share this repository at your own discretion.

> ⚠️ Note: If you build upon this for publication or product development, please consider citing the original tools (DrugEx, Vina) and this repository if relevant.


