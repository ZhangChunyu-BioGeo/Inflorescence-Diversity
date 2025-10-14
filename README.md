# README: Code for "Ecological and Evolutionary Drivers of Inflorescence Diversity in Chinese Angiosperms"

This repository contains the R code and supporting data files for the analyses presented in the manuscript:

> Shao, M., Zhang, C., Wang, R., Wang, H., & Zheng, P. (submitted to *Journal of Ecology*). Ecological and Evolutionary Drivers of Inflorescence Diversity in Chinese Angiosperms. 
>
> *Mingxue Shao and Chunyu Zhang contributed equally to this study. Hui Wang and Peiming Zheng contributed equally to this study.*

The code in this repository, written by Chunyu Zhang, implements the complete analytical pipeline for the manuscript. This includes data processing, null model testing, machine learning model comparisons (RF, GBM, etc.), phylogenetic comparative methods (ancestral state reconstruction), plant-pollinator co-occurrence analysis, and structural equation modeling (SEM).

All analyses were performed in R (version 4.4.2 or higher). The code is structured as a sequential workflow to ensure full reproducibility of the findings on the spatial patterns and drivers of inflorescence diversity across 25,778 Chinese angiosperm species.

## License

This project is licensed under the Apache License, Version 2.0. See the `LICENSE` file for full details.

## Repository Link

For peer review: `https://anonymous.4open.science/XXXXX`
Upon publication, the repository will be made publicly available on GitHub: `https://github.com/Chunyu-Zhang-BioGeo/Inflorescence-Diversity` (link to be finalized).

## System Requirements

* **Operating Systems**: Tested on Windows 10/11, macOS, and Linux.
* **Software Dependencies**: R (version 4.4.2 or higher).
* **Required R Packages** (with tested versions): `circlize` (0.4.15), `CoordinateCleaner` (Zizka et al. 2019), `corHMM` (2.8), `DALEX` (2.4.3), `dplyr` (1.1.2), `e1071` (1.7-14), `fastshap` (0.1.0), `gbm` (2.1.8), `ggplot2` (3.4.2), `ggparty` (1.0.0), `ggpubr` (0.6.0), `glmnet` (4.1-8), `iml` (0.11.1), `kknn` (1.3.1), `lavaan` (0.6-15), `nnet` (7.3-19), `partykit` (1.2-20), `phytools` (2.0), `randomForest` (4.7-1.1), `rgbif` (3.7.8), `rinat` (0.1.8), `tidyr` (1.3.0), `V.PhyloMaker2` (0.1.0).
* **Hardware**: A standard desktop computer with a multi-core CPU and at least 16 GB of RAM is recommended for the more computationally intensive steps (e.g., null models, SCM simulations).

## Installation Guide

1.  Install R from [https://cran.r-project.org/](https://cran.r-project.org/) and optionally RStudio from [https://posit.co/download/rstudio-desktop/](https://posit.co/download/rstudio-desktop/).
2.  Clone this repository or download and unzip the ZIP file.
3.  Open R/RStudio and set the working directory to the repository's root folder.
4.  Install all required R packages by running the following command in the R console:
    ```r
    install.packages(c("dplyr", "tidyr", "ggplot2", "ggpubr", "circlize", "lavaan",
                       "corHMM", "randomForest", "gbm", "e1071", "nnet", "glmnet",
                       "kknn", "iml", "DALEX", "fastshap", "partykit", "ggparty",
                       "V.PhyloMaker2"))
    ```
5.  Typical installation time for all packages is 5-15 minutes on a standard desktop computer withan internet connection.

## Reproduction Instructions

To fully reproduce all results and figures presented in the manuscript, follow these steps:

1.  **Data**:
All processed data files required to run the analysis scripts are provided in the `/data/` folder. Raw data were originally sourced from public databases as described in the manuscript's Methods section (Species 2000 China Node, Flora of China, GBIF, WorldClim, etc.).

2.  **Execute Analysis Workflow**: Run the R scripts located in the `/code/` folder in their numerical order (from `01` to `18`).
Each block of scripts corresponds to a major analytical section and figure in the paper.

* **Figure 1: Spatial Patterns & Environmental Filtering**
* `01-Calculate_Adjusted_Indeterminate_Inflorescence_and_Diversity.R`
* `02-LOO_Validation_for_Adjusted_Inflorescence_Diversity.R`
* `03-Null_Model_for_Inflorescence_Proportions_and_Diversity.R`

* **Figure 2: Machine Learning for Environmental Drivers**
* `04-ML_Models_P1_Generate_Simulated_Data.R`
* `05-ML_Models_P2_Process_Simulated_Data.R`
* `06-ML_Models_P3_Process_Observed_Data.R`
* `07-ML_Models_P4_Calculate_Means_and_P_Values.R
* `08.1-ML_Models_P5_Visualize_Results_Indeterminate_adj.R`
* `08.2-ML_Models_P6_Visualize_Results_H_Rank.R`

* **Figure 3:Evolutionary History & Paleoclimate**
* `09-Evolutionary_History_P1_Fit_Compare_Models.R`
* `10-Evolutionary_History_P2_SCM_Simulation_HR2_ARD.R`
* `11-Evolutionary_History_P3_Summarize_SCM_Results.R`
* `12-Evolutionary_History_P4_Extract_Paleoclimate_Data.R`
* `13-Evolutionary_History_P5_Visualize_SCM_Results.R`
* `14-Evolutionary_History_P6_Regression_Analysis_GMST.R

* **Figure 4: Co-occurrence with Potential Pollinators**
* `15-Pollinator_Cooccurrence_Strength_Observed.R`
* `16-Pollinator_Cooccurrence_Strength_Null_Model.R`
* `17-Pollinator_Cooccurrence_Strength_Visualization.R`

* **Figure 5: Synthesis via Structural Equation Modeling**
* `18-Structural_Equation_Modeling.R`

3.  **Expected Total Run Time**: The full reproduction may take several hours to complete on a standard multi-core desktop computer. The most time-intensive steps are the null model simulations (Script 16) and the stochastic character mapping simulations (Script 10), which involve thousands of iterations. Please be aware that fitting all 10 models in the 09-Evolutionary_History_P1_Fit_Compare_Models.R script, including the HR2_ARD model, is computationally intensive and may take several days to complete.

4.  **Outputs**: The execution of the workflow will generate figures identical to those in the manuscript (Figures 1-5) and print statistical summaries (e.g., model fit indices, p-values) to the Rconsole. Figures can be saved from the R plotting device or by adding `ggsave()` commands to the scripts.
