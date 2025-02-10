# 🧬 PRIME

[![DOI](https://zenodo.org/badge/DOI/your-doi-here.svg)](https://doi.org/your-doi-here)

> **Paper Citation:** Koskela et al. (2025). Integrative systems immunology analysis develops a new anti-PD1 prognostication model based on functional priming and reveals new targets involved in antigen presentation in melanoma. *Journal Name*. DOI: [link]

## 📋 Overview

This repository contains the analysis code and documentation for our paper investigating [brief description of your research]. The analysis primarily focuses on [main aspects of your analysis, e.g., "single-cell RNA sequencing analysis of immune cell populations in melanoma samples"].

## 🗂️ Repository Structure

```
.
├── data/           # Data files
├── R/              # Analysis scripts in R
└── Jupyter/        # Analysis scripts in Python as Jupyter notebbok files
   

## 🔧 Requirements

Python scripts use mostly standard libraries such as

- pandas
- numpy
- matplotlib
= scikit-learn
- scanpy
- celltypist

## 📊 Analysis Workflow

1. **Data Preprocessing** (`scripts/preprocessing/`)
   - `01_quality_control.py`: Initial QC and filtering
   - `02_normalization.py`: Data normalization and scaling

2. **Main Analysis** (`scripts/analysis/`)
   - `01_clustering.py`: Cell clustering analysis
   - `02_differential_expression.py`: DE analysis
   
3. **Visualization** (`scripts/visualization/`)
   - `01_umap_plots.py`: UMAP visualization scripts
   - `02_violin_plots.py`: Gene expression visualization

## 🔍 Data Availability

SAARA: Raw data is available at [repository/database name] under accession number [XXX]. Due to size constraints, this repository contains only processed data and analysis scripts.


## 📜 License

This project is licensed under the [appropriate license] - see the [LICENSE](LICENSE) file for details. 

## ✍️ Citation

If you use this code or data, please cite:

```bibtex
@article{Koskela2025,
  title={Integrative systems immunology analysis develops a new anti-PD1 prognostication model based on functional priming and reveals new targets involved in antigen presentation in melanoma},
  author={Koskela, Saara and Pulkkinen, Otto and },
  journal={Journal Name},
  year={2025},
  doi={your-doi}
}
```

## 🤝 Contributing

While this repository is primarily for archival purposes, if you find any issues or have suggestions, please open an issue in the GitLab repository.

## 📫 Contact

* **Principal Investigator:** Carlos Rogerio de Figueiredo (crdefi@utu.fi)
* **Lead Author:** Saara Koskela (saankos@utu.fi)

---
*Repository maintained by Arfa Mahmood and Otto Pulkkinen*
