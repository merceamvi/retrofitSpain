# *Beyond adoption rates: What lies behind uneven residential retrofit uptake? Evidence from Spain*
This repository contains the R code for the paper of the same title, **currently under review**.

## **Overview**
The analysis examines how residential retrofit adoption differs across household profiles in Spain and what accounts for these differences. The pipeline combines Factor Analysis of Mixed Data (FAMD) and Gaussian Mixture Models (GMM-VEV) to identify four multidimensional household profiles, a survey-weighted logit with Wald tests of common slopes across profiles, and a Fairlie decomposition of the adoption differences relative to the highest-adopting profile, with a Shapley attribution of the compositional component by block of characteristics and stratified bootstrap confidence intervals. Data come from the Spanish Survey on Income and Living Conditions (Encuesta de Condiciones de Vida, ECV, the Spanish component of EU-SILC) 2023 and the JRC AGRI4CAST Data Portal.

## **System and Software Information**
R version: developed under R ≥ 4.3.0.
Package management: all required packages are installed and loaded automatically by `_setup.R` and the preamble of `01_build_data.R`. Users do not need to install or manage dependencies manually.

## **Data Access**
Automatic download: `01_build_data.R` downloads the 2023 ECV microdata directly from the Spanish Statistical Office (Instituto Nacional de Estadística, INE). An active internet connection is required for this stage.
Climatic data: `HDD_CDD.xlsx` (included in this repository) contains regional (NUTS 2) heating and cooling degree days, computed as monthly averages over 2017–2023 from the [JRC AGRI4CAST Data Portal](https://agri4cast.jrc.ec.europa.eu/DataPortal/). It is merged with the INE microdata in stage 01.
Manual download (optional): the ECV microdata are also available from the [INE Living Conditions Survey page](https://www.ine.es/dyngs/INEbase/es/operacion.htm?c=Estadistica_C&cid=1254736176807&menu=resultados&idp=1254735976608); select the 2023 file. Paths are handled dynamically through the `here` package.
Stage 01 writes `data.RData`; if it is already present, stage 01 can be skipped.

## **Structure of the Repository**
The code is a staged pipeline. Each stage saves its intermediate objects to disk, so stages can be run individually or end to end. Tables and figures of the manuscript and appendix are written to `paper/`, named by their number in the paper; diagnostics and figures cited in the text without a table are written to `output/`.

1. `00_run_all.R` - master runner; sources stages 01 to 05 in order.
2. `_setup.R` - shared preamble for stages 02 to 05 (packages, working directory, publication ggplot2 theme, colour palettes, save helpers). Not a pipeline stage on its own.
3. `01_build_data.R` - downloads and cleans the 2023 ECV microdata, merges the climatic variables (HDD/CDD) and builds the analysis dataset (`data.RData`).
4. `02_descriptives_logit.R` - outcome distribution, descriptive statistics and the baseline survey-weighted logit with HC1-robust estimates and bootstrapped Average Marginal Effects (Tables 1, A.1, B.1, B.2 and Figure B.1).
5. `03_famd_clustering.R` - FAMD, multi-criteria assessment of GMM-VEV solutions across candidate dimensions (multi-start EM and bootstrap cluster-wise Jaccard stability), final partition and profile characterisation. Profiles are numbered by weighted adoption, from the highest (P1) to the lowest (P4) (Tables 3, C.1, C.2, C.3 and Figure C.1).
6. `04_profiles_ames.R` - profile-specific logit coefficients and bootstrapped AMEs, and observed versus predicted adoption by profile (Tables D.1, D.2 and Figure C.2).
7. `05_decomposition.R` - Wald tests of common slopes, aggregate decomposition of each adoption difference relative to P1, Shapley attribution of the compositional component and robustness to alternative coefficient specifications (Tables 4, 5, 6 and D.3).

Approximate run time: about 6 h in total, almost all of it in stage 03 (multi-start GMM and bootstrap stability). Stage 05 takes about 16 min with the default 500 bootstrap replications (`N_BOOT_DECOMP`).

## **How to Run the Code**
**Clone this repository:**
git clone https://github.com/merceamvi/retrofitSpain.git

Open the project in RStudio or your preferred R environment and run `00_run_all.R` from start to finish. The pipeline installs/updates packages, downloads the data, builds the dataset and generates all outputs (tables, figures, models).

To run a single stage instead, source `_setup.R` first (required by stages 02 to 05), since each stage loads the intermediate objects saved by the previous one.

## **Reproducibility and Support**
This pipeline has been designed to facilitate reproducibility and transparency. If you encounter any issues or have questions, please open an issue in this repository or contact the corresponding author of the article.

## **License and Citation**

This code is licensed under the [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).
You are free to copy, modify, and distribute the work as long as you provide appropriate credit to the original authors.

### Citation

If you use or adapt this code for your work, please cite:
*Amich, M., Arto, I., Mariel, P. (2026). Beyond adoption rates: What lies behind uneven residential retrofit uptake? Evidence from Spain. Manuscript under review.*

- **Mercè Amich**, University of the Basque Country (EHU), Avda. Lehendakari Aguirre, 83, E48015 Bilbao, Spain; and Basque Centre for Climate Change (BC3), Scientific Campus of the University of the Basque Country, Building 1, 1st floor, Sarriena s/n, E48940 Leioa, Spain. E-mail: [merce.amich@ehu.eus] (mailto:merce.amich@ehu.eus)
- **Iñaki Arto**, Basque Centre for Climate Change (BC3), Scientific Campus of the University of the Basque Country, Building 1, 1st floor, Sarriena s/n, E48940 Leioa, Spain. E-mail: [inaki.arto@bc3research.org](mailto:inaki.arto@bc3research.org)
- **Petr Mariel**, Department of Quantitative Methods, University of the Basque Country (EHU), Avda. Lehendakari Aguirre, 83, E48015 Bilbao, Spain. E-mail: [petr.mariel@ehu.es](mailto:petr.mariel@ehu.es), Tel: +34.94.601.3848, Fax: +34.94.601.3754

Corresponding author for code queries: merce.amich@ehu.eus (alternatively, merce.amich@bc3research.org)

Please ensure that proper attribution is given when citing or using the code in your research or projects.
