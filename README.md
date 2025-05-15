# 🧬 PRIME

[![DOI](https://zenodo.org/badge/DOI/your-doi-here.svg)](https://doi.org/your-doi-here)

> **Paper Citation:** Koskela et al. (2025). Novel Anti-PD1 Predictive Signature and Functional Dendritic-Cell Biomarkers in Melanoma Identified with Systems Immunology, *Journal Name*. DOI: [link]

## 📋 Overview

This repository contains scripts and preprocessed data (within size limitations, see below) related to the development and testing of the Prime model.

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

1. **R Codes** (`R/`)
   - Contain analyses related to circaplot visualization  
   - A script for integration and DEA of Gide, Riaz and Hugo count data

2. **Jupyter Notebooks** (`Jupyter/`)
   - Preprocessing of TPM-normalized data
   - Prime model development and testing the predictive performance
   - DEA and other analyses of Li cohort scRNA-seq data

## 🔍 Data Availability

Raw data is available at [repository/database name] under accession number [XXX]. Due to size constraints, this repository contains only processed data and analysis scripts.


## 📜 License

This project is licensed under the [appropriate license] - see the [LICENSE](LICENSE) file for details. 

## ✍️ Citation

If you use this code or data, please cite:

```bibtex
@article{Koskela2025,
  title={Novel Anti-PD1 Predictive Signature and Functional Dendritic-Cell Biomarkers in Melanoma Identified with Systems Immunology},
  author={Koskela, Saara and Pulkkinen, Otto and },
  journal={Journal Name},
  year={2025},
  doi={your-doi}
}
```

## 📫 Contact

* **Principal Investigator:** Carlos Rogerio de Figueiredo (crdefi@utu.fi)
* **Lead Author:** Saara Koskela (saankos@utu.fi)

---
*Repository maintained by Arfa Mahmood and Otto Pulkkinen*
