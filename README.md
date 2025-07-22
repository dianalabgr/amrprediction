## Machine Learning Predicts Antimicrobial Resistance from Genomic Data across ESKAPEE Pathogens

### 📖 Description of the Study

Antimicrobial resistance (AMR) is a mounting global crisis, fueled by the rapid emergence of multidrug-resistant bacteria. Among the most concerning culprits are the ESKAPEE bacteria — *Enterococcus faecium, Staphylococcus aureus, Klebsiella pneumoniae, Acinetobacter baumannii, Pseudomonas aeruginosa, Enterobacter spp.*, and *Escherichia coli* — which are leading causes of hospital-acquired infections worldwide.

In this study, we developed and validated machine learning models for predicting antimicrobial resistance phenotypes directly from genomic data. We assembled a robust dataset of **18,916 ESKAPEE genome assemblies**, each paired with its corresponding antibiogram, covering susceptibility results for **40 different antibiotics**.

Using this data, we trained **Random Forest** and **Extreme Gradient Boosting (XGBoost)** models for each antibiotic separately, consistently achieving over **90% recall and F1 score** for most pathogen–antibiotic combinations.

To maximize utility and accessibility:
- We developed an **interactive web platform**:  
  👉 https://dianalab.e-ce.uth.gr/amrpredictor/  
  allowing users to explore predictions and SHAP-derived feature importance.

- We validated our pipeline on **36 blood culture-positive ESKAPEE samples**, showing strong concordance with conventional phenotypic testing.

- We provide all data and models in Zenodo:  
  👉 https://zenodo.org/records/16213507

Our work underscores the transformative potential of integrating genomics and machine learning to deliver accurate, interpretable, and clinically actionable predictions against AMR.

---

### 📦 Repository Contents

This repository includes:

<ul>
  <li> Scripts for data acquisition and preprocessing</li>
  <li> Annotation tools for AMR genes, promoters, and rRNA features</li>
  <li> K-mer encoding pipelines</li>
  <li> Model training and evaluation scripts (Random Forest, XGBoost, DRO models)</li>
  <li> Clinical pipelines for metagenomic sample prediction</li>
</ul>

---

### 📂 Folder Overview (Workflow Order)

| Folder                   | Description                                                                                                                                     |
|-------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------|
| `data_acquisition/`     | Scripts for downloading, cleaning, and quality-controlling genome assemblies and matched antibiograms from public repositories.                 |
| `annotation/`           | Scripts for annotating genomes with AMR genes, upstream promoters, core genes, and rRNA regions; includes redundancy filtering and clustering.  |
| `maachine_learning/`    | Machine learning pipeline: k-mer encoding, model training (Random Forest, XGBoost), cross-validation, and evaluation across 40 antibiotics.     |
| `blood_cultures_predictions/` | Clinical pipeline to process metagenomic blood culture samples: taxonomic identification, assembly, k-mer extraction, and AMR prediction.    |
| `scripts/`              | Helper scripts for data formatting, statistics, figure creation, and reproducibility tasks.                                                     |
| `DRO_models/`          | Models trained using a Distributionally Robust Optimization (DRO)-inspired strategy to improve predictions for harder antibiotics.              |
| `annotate_newAssemblies/` | Tools to apply annotation steps on new or external assemblies, matching the training set format.                                              |

---

### 🔗 Summary

This repository provides a **start-to-finish workflow** for genome-driven AMR prediction:  
from data collection → annotation → machine learning → clinical metagenomic validation.

---

### 📜 License

This project is under the **MIT License**.

---

### 📫 Contact

For questions, collaborations, or licensing inquiries:  
📧 skulakis@gmail.com
