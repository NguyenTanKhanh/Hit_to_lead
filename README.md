# Hit-to-Lead Optimization Pipeline

This tool integrates **DrugEx** and **AutoDock Vina** into a streamlined pipeline for hit-to-lead optimization. It supports scaffold-based molecule generation with fragment attachment, and virtual screening via molecular docking.

---
![github_hit2lead](https://github.com/user-attachments/assets/52eb5392-fc83-438f-86c7-725a56b0f635)



## 🔧 Installation

Make sure you have [Conda](https://docs.conda.io/en/latest/) installed. Then follow these steps:

```bash
git clone https://github.com/NguyenTanKhanh/Hit_to_lead.git
cd Hit_to_lead

# Create and activate a new Conda environment
conda create -n Hit_to_lead python=3.11
conda activate Hit_to_lead


# Install DrugEx (latest from GitHub master)
pip install git+https://github.com/CDDLeiden/DrugEx.git@master

# Install Chemprop version 2.1.0 (requires Python >=3.11)
pip install chemprop==2.1.0

# Install openbabel
conda install -c conda-forge openbabel


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

###  Limitation 

This method uses an iodine (I) atom to temporarily block specific atoms in the scaffold during molecule growth. After generation, the iodine atoms are replaced back with hydrogen.
If your original scaffold already contains iodine atoms or if you intend to grow molecules that include iodine, this method will not work correctly.

---

## 📚 References

- [DrugEx – CDD Leiden](https://github.com/CDDLeiden/DrugEx)
- [VinaScreen – Yassir Boulaamane](https://github.com/yboulaamane/VinaScreen/tree/main)

---

## ⚖️ License

**This project has no formal license. It is provided freely for educational, academic, and non-commercial use.**  
Feel free to use, modify, and share this repository at your own discretion.

> ⚠️ Note: If you build upon this for publication or product development, please consider citing the original tools (DrugEx, Vina) and this repository if relevant.


