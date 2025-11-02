## Plasma proteome mediates the associations between air pollution exposure and disease risk

![graphical overview](./src/overall.png)

## Description
This repository provides a reproducible R codebase and figure scripts for a study on the air pollution–proteome–disease axis. It covers data derivation, proteome-wide association study (PWAS), time-to-event modeling for both air pollution and proteome predictors, causal mediation analysis identifying protein mediators, construction of an air-pollution–mediated proteomic risk score (APPRS), and predictive modeling via LASSO. The goal is to quantify the effects of ambient air pollution on disease risk and to dissect the proteomic pathways that mediate these effects.

![methods schematic](./src/workflow.png)

## Methods
We implement PWAS, Cox models for air pollution ~ disease outcomes and proteome ~ disease outcomes, mediation to identify potential protein mediators, APPRS construction and airpollution interaction analyses, and LASSO-based prediction in R.



## Assets
All data used in this study were obtained from the [UK Biobank](https://www.ukbiobank.ac.uk/).
- `DataDerive/`
  - `Outcome_BasicVariables.R`: derive outcomes and baseline variables.
  - `Other_Covirates.R`: additional covariates and preprocessing helpers.
- `StatisticalAnalysis/`
  - `PWAS_Analysis.R`: PWAS workflow (limma); batch across exposures.
  - `Proteome_DiseaseOutcome_Cox.R`: proteome-disease Cox.
  - `AP_DiseaseOutcome_Cox.R`: air pollution-disease Cox.
  - `MediationAnalysis.R`: mediation analysis.
  - `APPRS_LASSO.R`: build APPRS and compare three LASSO models.
  - `APPRS_AP_Diseases_Interaction.R`: Exposure×APPRS groups and Cox interaction.
- `Figures/`: scripts to main figures
- Dependencies (R 4.1): `data.table`, `dplyr`, `tidyr`, `survival`, `limma`, `mediation`, `glmnet`, `caret`, `pROC`, `ggplot2`, `ggpubr`, `GGally`, `ggradar`, `ggrepel`, `scales`, `openxlsx`.



## Citation
If you use this repository, please cite the associated manuscript once available. Provisional citation:
```
Li W, Li K, Zhou P, et al. Plasma proteome mediates the associations between air pollution exposure and disease risk. Year. Nat Commun. 
```
Manuscript under review
