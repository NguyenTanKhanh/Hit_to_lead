# Hit-to-Lead Optimization Pipeline

Hello, I am Khanh, a very handsome boy :)) 
<img src="https://flagcdn.com/w20/vn.png" width="20"> <img src="https://flagcdn.com/w20/vn.png" width="20"> <img src="https://flagcdn.com/w20/vn.png" width="20">



This project aims to integrate DrugEx with a new feature that preserves key hydrogen atoms, and to combine it with Chemprop and AutoDock Vina into a unified, streamlined pipeline for hit-to-lead optimization.
It supports **scaffold-based molecule generation** and **virtual screening via molecular docking**.

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
conda install -c conda-forge rdkit pandas numpy matplotlib scikit-learn openbabel
```

---

## 🚀 Pipeline Overview

### 1. Fragment Preparation

Generate dummy molecules from fragment combinations using **DrugEx**.

```bash
bash run_training.sh
```

---

### 2. Scaffold Preparation

Prepare a scaffold in **SMILES** format, visualize atom indices, and define attachment constraints.

```bash
bash run_scaffold_molecule_index.sh
```

After identifying atoms to block, apply modifications:

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

### 4. PAINS + BRENK filtering *(optional)*

Filter molecules like PAINS and BRENK using Rdkit

```bash
bash run_pains_brenk.sh
```
---
### 5. Property-Based Filtering *(optional)*

Filter molecules by predicted properties (e.g., **cell membrane and Golgi permeability**) to increase the chance of successful intracellular transport:

```bash
bash run_filter_pH_7.4.sh
bash run_filter_pH_6.3.sh
```

---

### 6. pKa-Based Filtering *(optional)*

Filter molecules by predicted properties (e.g., **cell membrane and Golgi permeability**) to increase the chance of successful intracellular transport:

```bash
bash run_pKa.sh
```

---

### 7. Activity-Based Filtering *(optional)*

Predict compound activity against the target protein:

```bash
bash run_predict_activity.sh
```

---

### 8. Molecular Docking with AutoDock Vina *(optional)*

Dock generated molecules against your target protein:

```bash
bash run_vinascreen.sh
```

---

### 9. Molecule Visualization *(optional)*

Visualize generated molecules as 2D images:

```bash
bash run_generate_molecule_image.sh
```

---

## 📚 References

- [DrugEx – CDD Leiden](https://github.com/CDDLeiden/DrugEx)  
- [VinaScreen – Yassir Boulaamane](https://github.com/yboulaamane/VinaScreen/tree/main)  
- [Chemprop-Esther Heid](https://github.com/chemprop/chemprop)
---

## ⚖️ License

This project has **no formal license** and is provided freely for **educational, academic, and non-commercial use**.  
Feel free to use, modify, and share this repository.

> ⚠️ Note: If you use this workflow for publication or product development, please consider citing the original papers ....
