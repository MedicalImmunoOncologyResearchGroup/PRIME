# 🧬 PRIME

[![DOI](https://zenodo.org/badge/DOI/your-doi-here.svg)](https://doi.org/your-doi-here)

> **Paper Citation:** Koskela et al. (2025). Novel Anti-PD1 Predictive Signature and Functional Dendritic-Cell Biomarkers in Melanoma Identified with Systems Immunology, *Journal Name*. DOI: [link]

## 📋 Overview

Abstract
Despite the success of anti-PD1 immunotherapy in melanoma, identifying biomarkers for response prediction and rational treatment combinations remains a major challenge. Effective dendritic cell (DC):CD8+ T cell crosstalk in the tumor microenvironment enhances anti-PD1 responses; however, no predictive signatures reflect this interaction. To address this, we developed a systems immunology approach based on DC:CD8+ T cell crosstalk reference biomarkers CD74 and CD8A, leading to the discovery of a 15-gene immune signature called PRIME, incorporated into a logistic regression framework for anti-PD1 response prediction in melanoma. PRIME outperforms previous signatures, especially in tumors with high CD8+ T cells, reflecting a functional immune response. Integrating clinical, molecular, spatial, and single-cell profiling identified SLAMF7 and TYMP as new functional biomarker candidates in DCs. Mechanistic studies revealed that these biomarkers regulate antigen processing and presentation in DCs. Our study demonstrates a new strategy for predicting anti-PD1 responses in melanoma, revealing biomarkers with functional roles that may inform future therapeutic development.


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
  title={ Novel Anti-PD1 Predictive Signature and Functional Dendritic-Cell Biomarkers in Melanoma Identified with Systems Immunology},
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
