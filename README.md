# 🧬 PRIME

[![DOI](https://zenodo.org/badge/DOI/your-doi-here.svg)](https://doi.org/your-doi-here)

> **Paper Citation:** Koskela et al. (2025). Novel Anti-PD1 Predictive Signature and Functional Dendritic-Cell Biomarkers in Melanoma Identified with Systems Immunology, *Journal Name*. DOI: [link]

## 📋 Overview

This repository contains scripts related to the development and testing of the Prime model, and codes to produce figures in the manuscript. 

## 🗂️ Repository Structure

```
.
├── data/           # Data files
├── R/              # Analysis scripts in R
└── Jupyter/        # Analysis scripts in Python as Jupyter notebbok files
   

## 🔧 Requirements

Python codes use mostly standard libraries such as

- pandas
- numpy
- seaborn
- matplotlib
- scikit-learn
- scanpy
- celltypist

Jupyter notebooks were tested with Python version 3.10.12.

## 📊 Analysis Workflow

1. **R Codes** (`R/`)
   - Contain analyses related to circos plot visualization  
   - A script for integration and DEA of Gide, Riaz and Hugo count data

2. **Jupyter Notebooks** (`Jupyter/`)
   - Preprocessing of TPM-normalized data
   - Prime model development and testing its predictive performance (response, OS, and PFS)
   - DEA and other analyses of Li cohort scRNA-seq data
   - Anti-MAA analysis

## 🔍 Data Availability

Data to run the codes are available to reviewers only.


## 📜 License

This project is licensed under the XXX - see the [LICENSE](LICENSE) file for details. 

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
