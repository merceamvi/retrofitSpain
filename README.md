# **Structured heterogeneity in residential retrofit adoption: profile-specific evidence from Spain**
This repository contains the R script developed for a paper currently Under Review.

## **System and Software Information**
R Version: The script was built using R version 4.3.1 (2023-06-16 ucrt).  
Package Management: All required R packages are automatically installed and updated as part of the script’s "Preliminaries" section.  
Users do not need to manually install or manage dependencies.  
The script ensures that the R packages are installed and loaded for all users.  

## **Data download**
The script automates the download of data from the public database of the Spanish Statistical Office (Instituto Nacional de Estadística). 
Users do not need to manually locate or download data, as the script handles the process entirely.

## **Structure of the Script**
The script is organized into six sections for clarity:  

1. [Objective](https://github.com/merceamvi/cleanvehiclesSpain/blob/8788fb15821c4d2f59032ee9acd0ac6eaf38d6de/Script_CleanCars.R#L8): An overview of the analysis objectives and research questions addressed in the article  
2. [Preliminaries](https://github.com/merceamvi/cleanvehiclesSpain/blob/8788fb15821c4d2f59032ee9acd0ac6eaf38d6de/Script_CleanCars.R#L15): Automatic installation or updating of required packages, loading libraries, and setting up the working environment  
3. [Download Data and Preprocessing](): Automated data download, data cleaning, and preparation for analysis  
4. [Descriptive Statistics and Frequency Tables](https://github.com/merceamvi/cleanvehiclesSpain/blob/8788fb15821c4d2f59032ee9acd0ac6eaf38d6de/Script_CleanCars.R#L527): Summary statistics and frequency tables for the variables of interest, providing an overview of the dataset  
5. [Model Estimation](https://github.com/merceamvi/cleanvehiclesSpain/blob/8788fb15821c4d2f59032ee9acd0ac6eaf38d6de/Script_CleanCars.R#L731): Estimation of an Error-Component Mixed Logit model using the Apollo package  
6. [Post-Estimation](https://github.com/merceamvi/cleanvehiclesSpain/blob/8788fb15821c4d2f59032ee9acd0ac6eaf38d6de/Script_CleanCars.R#L1276): AME computation + CIs 95% via 1k parametric MC simulation and final plot.  

## **Data Access**
Download Instructions: The script automatically downloads the required datasets as part of the preprocessing step. Ensure an active internet connection when running the script.  
Manual Download: If needed, datasets can also be downloaded manually from the following [URL](https://www.ine.es/dyngs/INEbase/es/operacion.htm?c=Estadistica_C&cid=1254736177092&menu=ultiDatos&idp=1254735572981)  
Manually downloaded files should be placed in the appropriate directory as specified in the script. The here package is used to manage paths dynamically.  

## **How to Run the Script**
This script is designed for one-click execution. To run:  

**Clone this repository:**
git clone https://github.com/merceamvi/cleanvehiclesSpain.git  
Open the R script ([Script_Cleancars].R) in RStudio or your preferred R environment.  

Execute the script from start to finish.   
The script automatically installs/updates packages, downloads data, and generates all outputs (tables, figures, and models) in the designated output directories.  

## **Reproducibility and Support**
This script has been carefully designed to facilitate reproducibility and transparency.   
If you encounter any issues or have questions, please open an issue in this repository or contact the corresponding author of the article. 

## **License and Citation**

This script is licensed under the [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).  
You are free to copy, modify, and distribute the work as long as you provide appropriate credit to the original authors.

### Citation

If you use or adapt this script for your work, please cite:
*Amich, M., Tomás, M., Mariel, P., & Arto, I. (2026). Clean vehicle ownership: Implications for effective policy interventions. Energy Economics, In Press. https://doi.org/10.1016/j.eneco.2026.109172*

- **Mercè Amich**, Basque Centre for Climate Change (BC3), Scientific Campus of the University of the Basque Country (UPV/EHU), Building 1, 1st floor, Sarriena s/n, E48940 Leioa, Spain; and University of the Basque Country (UPV/EHU), Avda. Lehendakari Aguirre, 83, E48015 Bilbao, Spain. E-mail: [merce.amich@bc3research.eus](mailto:merce.amich@bc3research@ehu.eus)
- **Manuel Tomás**, Basque Centre for Climate Change (BC3), Scientific Campus of the University of the Basque Country, Building 1, 1st floor, Sarriena s/n, E48940 Leioa, Spain. E-mail: [manuel.tomas@bc3research.org](mailto:manuel.tomas@bc3research.org)
- **Iñaki Arto**, Basque Centre for Climate Change (BC3), Scientific Campus of the University of the Basque Country, Building 1, 1st floor, Sarriena s/n, E48940 Leioa, Spain. E-mail: [inaki.arto@bc3research.org](mailto:inaki.arto@bc3research.org)
- **Petr Mariel**, Department of Quantitative Methods, University of the Basque Country (UPV/EHU), Avda. Lehendakari Aguirre, 83, E48015 Bilbao, Spain. E-mail: [petr.mariel@ehu.es](mailto:petr.mariel@ehu.es), Tel: +34.94.601.3848, Fax: +34.94.601.3754

Corresponding author for code queries: merce.amich@bc3research.org

Please ensure that proper attribution is given when citing or using the script in your research or projects.



# cleanvehiclesSpain

