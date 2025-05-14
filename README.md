# Hit2Lead_K_test
 
This is process to optimize hit to lead by intergated DrugEx and Vina

Step 1: Run 1.run_training.sh to train the generation model

Step 2: Run 2.run_lead_generation.sh to create new molecules containing scaffold

Step 2.1: Run python generate_molecule_image.py to see image of generation_compound

Step 3: Run 3.run_vinascreen.sh to dock list compounds to target protein


## Installation


```bash
git clone https://github.com/NguyenTanKhanh/Hit_to_lead.git
cd Hit_to_lead
conda env create -f env.yml
conda activate hit2lead
```

### Hardware Requirements

I build based the DrugEx and Vina for optimze hit2lead process. 

Ref: 
1. https://github.com/CDDLeiden/DrugEx
2. https://github.com/yboulaamane/VinaScreen/tree/main



