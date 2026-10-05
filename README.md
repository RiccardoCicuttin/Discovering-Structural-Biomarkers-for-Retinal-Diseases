# Discovering Structural Biomarkers for Retinal Diseases

In this project retinal layer thickness profiles extracted from Optical Coherence Tomography (OCT) scans are used to identify the structural signatures that separate four retinal pathologies (AMD, CSR, DR, MH) from healthy controls.


## Repository structure
```text
Biomarker-discovery-for-retinal-diseases/
│
├── data/
│   └── PatientsData.xlsx              # raw dataset; not uploaded
│
├── R/
│   ├── packages.R                     # libraries used across the project
│   ├── feature_utils.R                # feature engineering on the raw profiles
│   ├── fda_utils.R                    # smoothing and functional plotting helpers
│   ├── model_defs.R                   # nested-CV fit function, one per model
│   └── eval_utils.R                   # evaluation metrics and result tables
│
├── scripts/
│   ├── 01_preprocessing.R             # train/test split and feature datasets
│   ├── 02_eda_multivariate.R          # exploratory analysis of the layer means
│   ├── 03_covariance_analysis.R       # inter-layer covariance analysis
│   ├── 04_baseline_models.R           # model selection by nested cross-validation
│   ├── 05_eda_functional.R            # smoothing, MFPCA, functional separation
│   ├── 06_fda_plots.R                 # figures for the functional analysis
│   └── 07_restricted_models.R         # restricted feature sets, test-set evaluation
│
├── datasets/                          # processed .rds datasets
│
├── results/                           # model outputs and .rds results
│
├── figures/                           # generated figures
│
├── .gitignore
├── Functional-Data-Analysis-for-Retinal-Disease-Classification.Rproj
└── README.md
```

The raw dataset is not public. `datasets/`, `results/` and `figures/` are generated locally and are not committed.

Full results and analysis are detailed in the [project report](report/Biomarker-Discovery-for-Retinal-Diseases.pdf).
