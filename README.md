## AirPollution-Proteome-Diseases-Mediation

![graphical overview](./src/overall.png)

## Description
This repository provides a reproducible R codebase and figure scripts for a study on the air pollution–proteome–disease axis. It covers data derivation, proteome-wide association (PWAS), time-to-event modeling for both air pollution and proteome predictors, causal mediation analysis identifying protein mediators, construction of an air-pollution–mediated proteomic risk score (ApPRS), and predictive modeling via LASSO. The goal is to quantify the effects of ambient air pollution on disease risk and to dissect the proteomic pathways that mediate these effects.

![methods schematic](./src/workflow.png)

## Methods
We implement PWAS, Cox models for air pollution and proteomic predictors, causal mediation to identify protein mediators, ApPRS construction and interaction analyses, and LASSO-based prediction in R.



## Assets
- `DataDerive/`
  - `Outcome_BasicVariables.R`: derive outcomes and baseline variables.
  - `Other_Covirates.R`: additional covariates and preprocessing helpers.
- `StatisticalAnalysis/`
  - `PWAS_Analysis.R`: PWAS workflow (limma); batch across exposures.
  - `Proteome_DiseaseOutcome_Cox.R`: proteome→disease Cox (scalable public API).
  - `AP_DiseaseOutcome_Cox.R`: air pollution→disease Cox (IQR scaling; start-date truncation).
  - `MediationAnalysis.R`: single/batch mediation and overlap-based mediator filtering.
  - `APPRS_LASSO.R`: build ApPRS and compare three LASSO models.
  - `APPRS_AP_Diseases_Interaction.R`: exposure×ApPRS groups and Cox interaction.
- `Figures/`: scripts to reproduce main figures (`Figure1.R`–`Figure5.R`).
- Dependencies (R): `data.table`, `dplyr`, `tidyr`, `survival`, `limma`, `mediation`, `glmnet`, `caret`, `pROC`, `ggplot2`, `ggpubr`, `GGally`, `ggradar`, `ggrepel`, `scales`, `openxlsx`.



## Citation
If you use this repository, please cite the associated manuscript once available. Provisional citation:
```
[Authors]. Air pollution, proteomic profiles, and disease mediation: a reproducible analysis repository. Year. (Manuscript under review.)
```
Replace with the final bibliographic details upon publication.
