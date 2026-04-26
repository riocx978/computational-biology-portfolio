# CompBioCodes

> A working collection of R and Python scripts spanning my experience in computational biology, statistical genetics, and machine learning — from data cleaning through to publication-level analyses.

This repo is intentionally broad. It reflects the range of methods I've applied across research projects at the University of South Florida, including EWAS, meQTL mapping, differential expression, and causal inference.

---

## What's here

### 📊 Differential Gene Expression
| File | Description |
|------|-------------|
| `celegans_deg_analysis.R` | DEG analysis pipeline using DESeq2 — normalization, statistical testing, volcano plots | Extended DEG workflow with additional visualization and filtering steps |

---

### 🤖 Machine Learning & Statistical Modeling
| File | Description |
|------|-------------|
| `heart_disease_prediction.R` | Benchmarks multiple ML models (random forest, SVM, logistic regression) with cross-validation |
| `covid_viral_classification.R` | Cross-validation framework for model selection and performance evaluation |
| `viral_status_classification.R` | Logistic regression and LDA applied to biological classification problems |

---

### 🐍 Python: Data Processing
| File | Description |
|------|-------------|
| `protein_disorder_pipeline.py` | Protein disorder analysis and separation metrics for high-dimensional biological data |

---

## Skills demonstrated

`R` `Python` `DESeq2` `TwoSampleMR` `scikit-learn` `ggplot2` `tidyverse` `Mendelian Randomization` `differential expression` `machine learning` `cross-validation` `statistical modeling` `data visualization`

---

## Related repositories

| Repo | What it contains |
|------|-----------------|
| [Mendelian-Randomization](https://github.com/riocx978/Mendelian-Randomization) | Full MR analysis — published Master's thesis on periodontal disease |
| [MeQTL-Mapping](https://github.com/riocx978/MeQTL-Mapping) | HPC pipeline for extracting top meQTL signals per CpG site |

---

## Requirements

**R ≥ 4.0**
```r
install.packages(c("tidyverse", "DESeq2", "glmnet", "TwoSampleMR", "MendelianRandomization", "caret"))
```

**Python ≥ 3.8**
```bash
pip install numpy pandas scikit-learn matplotlib seaborn
```

---

## Contact

**Rhea Charles** · [riocx1997@gmail.com](mailto:riocx1997@gmail.com) · [LinkedIn](https://www.linkedin.com/in/rhea-charles/) · [GitHub](https://github.com/riocx978)
