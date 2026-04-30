# COMP 7745 Machine Learning

## Project

**Predicting Cancer Drug Sensitivity from Gene Expression using Machine Learning**

This project builds a regression pipeline to predict Erlotinib drug response (IC50 / LN_IC50) from CCLE gene-expression features using classical ML models.

## Repository Structure

- `notebooks/01_data_ingestion_and_alignment.ipynb` - load and align CCLE + GDSC2 data
- `notebooks/02_preprocessing_and_feature_reduction.ipynb` - split, variance filter, scale, PCA
- `notebooks/03_baseline_and_linear_models.ipynb` - Dummy, Ridge, Lasso, ElasticNet, SVR
- `notebooks/04_tree_boosted_models.ipynb` - Random Forest, XGBoost
- `notebooks/05_evaluation_and_diagnostics.ipynb` - test metrics and diagnostic plots
- `notebooks/06_feature_importance_and_biological_interpretation.ipynb` - feature importance back-projection and pathway checks
- `notebooks/utils/` - helper utilities for data and model metrics
- `reports/tables/` - model comparison and interpretation outputs
- `reports/figures/` - diagnostic plots

## Dataset Download Links

### 1) GDSC2 Drug Sensitivity Data (CSV)

Use the Cell Model Passports downloads page:  
[https://cellmodelpassports.sanger.ac.uk/downloads](https://cellmodelpassports.sanger.ac.uk/downloads)

Steps:
1. Open **Drug Sensitivity Data**
2. Click **Download GDSC2 IC50 data**
3. Click **View all versions**
4. Select and download your target CSV (the file you choose for this project should be saved as `GDSC2_fitted_dose_response_27Oct23.csv` in the `Data/` folder)

### 2) CCLE Expression Data (GCT/RES source)

Use DepMap data page:  
[https://depmap.org/portal/data_page/?releasename=mRNA+expression&filename=CCLE_Expression_2012-09-29.res&tab=allData](https://depmap.org/portal/data_page/?releasename=mRNA+expression&filename=CCLE_Expression_2012-09-29.res&tab=allData)

Steps:
1. Locate **CCLE_Expression_2012-09-29.res**
2. Download the expression dataset in **GCT** format
3. Save it as `CCLE_Expression_Entrez_2012-09-29.gct` in the `Data/` folder

General DepMap downloads landing page:  
[https://depmap.org/portal/data_page/?tab=allData](https://depmap.org/portal/data_page/?tab=allData)

## Project Setup

### 1) Create and activate environment (recommended)

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -r requirements.txt
```

### 2) Place raw datasets

Create a `Data/` folder at project root and place:

- `Data/CCLE_Expression_Entrez_2012-09-29.gct`
- `Data/GDSC2_fitted_dose_response_27Oct23.csv`

## How to Run

Run notebooks in order:

1. `notebooks/01_data_ingestion_and_alignment.ipynb`
2. `notebooks/02_preprocessing_and_feature_reduction.ipynb`
3. `notebooks/03_baseline_and_linear_models.ipynb`
4. `notebooks/04_tree_boosted_models.ipynb`
5. `notebooks/05_evaluation_and_diagnostics.ipynb`
6. `notebooks/06_feature_importance_and_biological_interpretation.ipynb`

You can run from JupyterLab/Notebook or execute via VS Code/Cursor notebook runner.

## Expected Outputs

After execution, outputs are written to:

- `data/processed/` - aligned features/targets, train-val-test splits, PCA features
- `models/` - fitted preprocessing objects and trained models (`.joblib`)
- `reports/tables/` - CSV result tables (model comparison, top genes, pathway ranks)
- `reports/figures/` - residual and predicted-vs-actual plots

## Notes

- `notebooks/utils/data_utils.py` supports both `.csv` and Excel (`.xlsx`/`.xls`) for GDSC input.
- Drug filtering is set to `Erlotinib` in Notebook 01.
- Keep notebook execution order unchanged to ensure downstream artifacts exist.
