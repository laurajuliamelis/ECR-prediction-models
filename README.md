# EOR-prediction-FEP

**Early Onset Response Prediction in First Episode Psychosis using Clinical and Genetic Data**

This repository contains the full analysis pipeline used to develop predictive models for early onset response (EOR) to antipsychotic treatment in patients with a first episode of psychosis (FEP). The study integrates clinical and polygenic data and applies both traditional statistical modelling and machine learning techniques.

---

## 📚 Project Overview

The goal of this project is to identify early predictors of treatment response by developing a clinically interpretable model. The model includes baseline clinical variables, cognitive assessments, and polygenic risk scores (PRS) for obstetric, psychiatric and biological traits.

- **Cohort**: PEPs (n = 236)
- **Target variable**: Early Onset Response (EOR) at 2 months
- **Predictors**: Baseline clinical variables, neurocognition, and 7 PRSs

---

## 🧪 Analysis Pipeline

The analysis is structured across 6 scripts, executed in sequential order:

| Script | Description |
|--------|-------------|
| `1_data_preprocessing_and_descriptives.R` | Initial data cleaning, exclusion criteria, and exploratory analyses |
| `2_data_imputation.R` | Multiple imputation of missing values using MICE |
| `3_variable_selection_LASSO.R` | Feature selection via logistic LASSO with 10-fold cross-validation |
| `4_model_fitting.R` | Logistic regression and ML models fitted on imputed data (clinical, genetic, combined) |
| `5_model_evaluation.R` | Model calibration, discrimination (AUC), bootstrap-based comparisons and variable contribution|
| `6_decision_curve_analysis.R` | Evaluation of clinical utility using decision curve analysis (DCA) |

---

## 📦 Required Packages

Key R packages used include:

- `openxlsx`, `mice`, `VIM`, `naniar`
- `glmnet`, `dplyr`
- `rms`, `caret`, `pROC`, `boot`, `gbm`
- `ggplot2`, `ggsignif`, `fastshap`, `shapviz`
- `dcurves`, `gtsummary` 

You can install missing packages using:

```r
packages <- c("openxlsx", "mice", "VIM", "naniar", "glmnet", "dplyr", "rms", "caret", "pROC", "boot", "gbm", "ggplot2", "ggsignif", "fastshap", "shapviz", "dcurves", "gtsummary")
install.packages(setdiff(packages, rownames(installed.packages())))
```

---

## ⚙️ Reproducibility

To reproduce the full analysis:

1. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/EOR-prediction-FEP.git
   ```

2. Open R or RStudio and set the working directory to the repo folder.

3. Run the scripts in the order specified above.

*Note: Raw data are not included in this repository due to privacy restrictions.*

---

## 📄 License

This code is released under the [MIT License](LICENSE).

---

## 📬 Contact

For questions, suggestions, or collaboration opportunities, please contact:

**Laura Julià Melis**  
Bioestatistician · Fundació Clínic-IDIBAPS  
PhD Candidate · University of Barcelona  
📧 laurajulia@ub.edu

and 

**Sergi Mas Herrero**  
Accredited Researcher · Fundació Clínic-IDIBAPS

Associate Professor · University of Barcelona  
📧 sergimash@ub.edu
