# **Structured heterogeneity in residential retrofit adoption: profile-specific evidence from Spain**
This repository contains the R code for the paper of the same title, **currently under review**.

## **Overview**
The analysis identifies latent household profiles in Spain and estimates how the determinants of residential retrofit adoption differ across them. The pipeline combines Factor Analysis of Mixed Data (FAMD) for dimensionality reduction, Gaussian Mixture Modelling (GMM) to recover four household profiles, and profile-specific weighted logit models with bootstrapped Average Marginal Effects (AMEs). Data come from the Spanish Living Conditions Survey (ECV / EU-SILC) 2023.

## **System and Software Information**
R version: developed under R ≥ 4.3.0.
Package management: all required packages are installed and loaded automatically by `_setup.R` and the preamble of `01_build_data.R`. Users do not need to install or manage dependencies manually.

## **Data Access**
Automatic download: `01_build_data.R` downloads the 2023 ECV microdata directly from the Spanish Statistical Office (Instituto Nacional de Estadística, INE). An active internet connection is required for this stage.
Required local file: place `HDD_CDD.xlsx` (heating and cooling degree days) in the repository directory before running. It is merged with the INE microdata in stage 01.
Manual download (optional): the ECV microdata are also available from the [INE Living Conditions Survey page](https://www.ine.es/dyngs/INEbase/es/operacion.htm?c=Estadistica_C&cid=1254736176807&menu=resultados&idp=1254735976608); select the 2023 file. Paths are handled dynamically through the `here` package.
Stage 01 writes `data.RData`; if it is already present, stage 01 can be skipped.

## **Structure of the Repository**
The code is a staged pipeline. Each stage saves its intermediate objects to disk, so stages can be run individually or end to end.

1. `00_run_all.R` - master runner; sources stages 01 to 04 in order.
2. `_setup.R` - shared preamble for stages 02 to 04 (packages, working directory, publication ggplot2 theme, colour palettes, save helpers). Not a pipeline stage on its own.
3. `01_build_data.R` - downloads and cleans the 2023 ECV microdata, merges climatic variables (HDD/CDD), and builds the analysis dataset (`data.RData`).
4. `02_descriptives_logit.R` - descriptive statistics and the pooled weighted logit adoption model, with bootstrapped AMEs (Figure 1 and supporting tables).
5. `03_famd_clustering.R` - FAMD and GMM clustering recovering four household profiles (scree plot, profile-characterisation table).
6. `04_profiles_policy.R` - profile-specific analysis: adoption-rate heterogeneity, the effect-heterogeneity likelihood-ratio test, profile-specific bootstrapped AMEs, and the profile-faceted AME figure (Figure 3).

## **How to Run the Code**
**Clone this repository:**
git clone https://github.com/merceamvi/rehabSpain.git

Place `HDD_CDD.xlsx` in the repository directory. Open the project in RStudio or your preferred R environment and run `00_run_all.R` from start to finish. The pipeline installs/updates packages, downloads the data, builds the dataset, and generates all outputs (tables, figures, models) in the working directory.

To run a single stage instead, source `_setup.R` first (required by stages 02 to 04), since each stage loads the intermediate objects saved by the previous one.

## **Reproducibility and Support**
This pipeline has been designed to facilitate reproducibility and transparency. If you encounter any issues or have questions, please open an issue in this repository or contact the corresponding author of the article.

## **License and Citation**

This code is licensed under the [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).
You are free to copy, modify, and distribute the work as long as you provide appropriate credit to the original authors.

### Citation

If you use or adapt this code for your work, please cite:
*Amich, M., Arto, I., Mariel, P. (2026). Structured heterogeneity in residential retrofit adoption: profile-specific evidence from Spain. Manuscript under review.*

- **Mercè Amich**, Basque Centre for Climate Change (BC3), Scientific Campus of the University of the Basque Country (UPV/EHU), Building 1, 1st floor, Sarriena s/n, E48940 Leioa, Spain; and University of the Basque Country (EHU), Avda. Lehendakari Aguirre, 83, E48015 Bilbao, Spain. E-mail: [merce.amich@bc3research.org](mailto:merce.amich@bc3research.org)
- **Iñaki Arto**, Basque Centre for Climate Change (BC3), Scientific Campus of the University of the Basque Country, Building 1, 1st floor, Sarriena s/n, E48940 Leioa, Spain. E-mail: [inaki.arto@bc3research.org](mailto:inaki.arto@bc3research.org)
- **Petr Mariel**, Department of Quantitative Methods, University of the Basque Country (EHU), Avda. Lehendakari Aguirre, 83, E48015 Bilbao, Spain. E-mail: [petr.mariel@ehu.es](mailto:petr.mariel@ehu.es), Tel: +34.94.601.3848, Fax: +34.94.601.3754

Corresponding author for code queries: merce.amich@bc3research.org

Please ensure that proper attribution is given when citing or using the code in your research or projects.
