###############################################################################
# Project       : [rehabSpain] [1] DOWNLOAD DATA & MERGE FUSED AND CLIMATIC
# Creation date : 02/12/2024
# Last update   : 26/06/2026
# Author        : Mercè Amich (merce.amich@ehu.eus)
# Institution   : UPV/EHU, BC3
# Last run time : 2.06 min.

# Script Overview:
#   This script processes SILC-23 data from INE, merges additional data 
#   (HDD/CDD), performs preliminary transformations, and prepares the data
#   for binary logit + FAMD + GMM

# Requirements:
#   - This script must be in the same directory as:
#       b) HDD_CDD.xlsx

# Output  : data.RData

###############################################################################

# =============================== PRELIMINARIES ==============================

rm(list = ls(all = TRUE)) # Clear workspace
Sys.setenv(LANG = "en")   # Set system language
start.time <- Sys.time()  # Track start time

# ---- Load Required Packages ----
packages.needed <- c(
  "here", "rvest", "httr", "dplyr", "ggplot2", "pscl", "AER", "scales",
  "survey", "gridExtra", "readxl", "effects", "foreign", "aod", "lmtest",
  "DescTools", "stargazer", "Hmisc", "margins", "plotly", "lme4", "ggeffects",
  "patchwork", "marginaleffects", "forcats", "rnaturalearth", "sf", "mgcv",
  "gratia", "glmx", "ResourceSelection", "pROC"
)

# Install & load packages
packages.loaded <- installed.packages()
for (p in packages.needed) {
  if (!p %in% row.names(packages.loaded)) install.packages(p)
  library(p, character.only = TRUE)
}

# Define working directory
path <- here()
setwd(path)

# ============================ 1. SILC-23 DATA ==============================

# ---- a) Download and unzip ----
# Explanation: Downloads microdata for 2023 from INE, unzips and cleans files

page_url <- "https://www.ine.es/dyngs/INEbase/es/operacion.htm?c=Estadistica_C&cid=1254736176807&menu=resultados&idp=1254735976608#_tabs-1254736195153"
page <- read_html(page_url)

option_selector <- "//select[@id='ir']//option[contains(text(),'2023')]"
option_element  <- page %>% html_nodes(xpath = option_selector)
file_path       <- option_element %>% html_attr("value")

base_url <- "https://www.ine.es"
file_url <- paste0(base_url, file_path[1])  # Take only first link

output_dir <- "DATA" # Create output directory
if (!dir.exists(output_dir)) dir.create(output_dir)

# Download & unzip
file_name <- paste0(output_dir, "/microdatos_2023.zip")
download.file(file_url, destfile = file_name, mode = "wb")
unzip(file_name, exdir = output_dir)
unlink(file_name)

# Unzip any remaining zip files in the directory
zip_files <- list.files(output_dir, pattern = "\\.zip$", full.names = TRUE)
for (zip_file in zip_files) {
  unzip(zip_file, exdir = output_dir)
  unlink(zip_file)
}

# ---- b) Load RData files ----
# Explanation: Load individual hh & personal datasets, rename, and clean envir.

load_and_clean_data <- function(file_path, new_name) { # create function
  load(file_path)
  assign(new_name, Microdatos, envir = .GlobalEnv)
  rm(Microdatos, Metadatos)
}

datasets <- list( 
  list(file = paste0(path, "/DATA/R/ECV_Td_2023.RData"), name = "Td"),
  list(file = paste0(path, "/DATA/R/ECV_Th_2023.RData"), name = "Th"),
  list(file = paste0(path, "/DATA/R/ECV_Tp_2023.RData"), name = "Tp"),
  list(file = paste0(path, "/DATA/R/ECV_Tr_2023.RData"), name = "Tr")
)

for (dataset in datasets) load_and_clean_data(dataset$file, dataset$name)

# Remove unnecessary folders
folders_to_remove <- c("CSV", "SAS", "SPSS", "STATA")
for (folder in folders_to_remove) {
  folder_path <- file.path(output_dir, folder)
  if (dir.exists(folder_path)) unlink(folder_path, recursive = TRUE)
}

# Keep only relevant objects
rm(list = setdiff(ls(), c("Td", "Th", "Tp", "Tr", "path", "start.time")))

# ---- c) Merge Household & Individual data objects by shared key ----
household  <- left_join(Td, Th, by = c("DB030" = "HB030"))
individual <- left_join(Tr, Tp, by = c("RB030" = "PB030"))

# Keep only merged data frames, path and timing
rm(list = setdiff(ls(), c("household", "individual", "path", "start.time")))

# ---- d) Select and summarise variables ----
# Explanation: Select key hh & individual-level variables, create dummies
# and collapse individual-level variables to household level.

# Household:
household <- household %>%
  dplyr::select(
    DB030      # Household identifier
    , DB090    # Weights
    , DB040    # CCAA
    , DB100    # Urbanisation grade
    , HB120    # Number of household members
    , vhRentaa # Income 
    , HH010    # Dwelling type
    , HH021    # Dwelling tenure regime
    , HH030    # Number of rooms
    , HH040    # Leaks, dampness... in walls, floors...
    , HX060    # Household type
    , HC020    # Square meters dwelling
    , HS160    # Insufficient natural light
    , HS170    # Noises (industries, street...)
    , HS180    # Environmental pollution in the zone
    , HC001    # Heating type
    
    , HC060    # Enough temperature in winter
    , HC070    # Enough temperature in summer
    , HS022    # Social bonus
    
    , HH050    # Affordability enough temp. winter
    , vhMATDEP # Material deprivation
    , HC080    # Dwelling satisfaction
    
    , HH060    # Rental costs for main dwelling    
    , HH070    # Housing expenses  
    
    , HC002    # Main energy source dwelling
    , HC003    # Retrofit-renovation measures count
  ) %>%
  
  # Group by household identifier
  group_by(DB030) %>%  
  
  # Indications on how to collapse variables
  dplyr::summarize(
    weight         = as.numeric(DB090),
    
    NUTS2 = DB040,
    galicia        = ifelse(any(DB040 == "ES11"), 1, 0),
    asturias       = ifelse(any(DB040 == "ES12"), 1, 0),
    cantabria      = ifelse(any(DB040 == "ES13"), 1, 0),
    paisvasco      = ifelse(any(DB040 == "ES21"), 1, 0),
    navarra        = ifelse(any(DB040 == "ES22"), 1, 0),
    larioja        = ifelse(any(DB040 == "ES23"), 1, 0),
    aragon         = ifelse(any(DB040 == "ES24"), 1, 0),
    madrid         = ifelse(any(DB040 == "ES30"), 1, 0),
    castillaleon   = ifelse(any(DB040 == "ES41"), 1, 0),
    castillamancha = ifelse(any(DB040 == "ES42"), 1, 0),
    extremadura    = ifelse(any(DB040 == "ES43"), 1, 0),
    catalunya      = ifelse(any(DB040 == "ES51"), 1, 0),
    comvalenciana  = ifelse(any(DB040 == "ES52"), 1, 0),
    illesbalears   = ifelse(any(DB040 == "ES53"), 1, 0),
    andalucia      = ifelse(any(DB040 == "ES61"), 1, 0),
    murcia         = ifelse(any(DB040 == "ES62"), 1, 0),
    ceuta          = ifelse(any(DB040 == "ES63"), 1, 0),
    melilla        = ifelse(any(DB040 == "ES64"), 1, 0),
    canarias       = ifelse(any(DB040 == "ES70"), 1, 0),
    extraregio     = ifelse(any(DB040 == "ESZZ"), 1, 0),
    
    bonosocial     = ifelse(any(HS022 == "1"), 1, 0),
    naturalight    = ifelse(any(HS160 == "1"), 1, 0),
    noises         = ifelse(any(HS170 == "1"), 1, 0),
    pollution      = ifelse(any(HS180 == "1"), 1, 0),
    
    pagar.enoughtemp.wint = ifelse(any(HH050 == "1"), 1, 0),
    material.depriv       = ifelse(any(vhMATDEP == "1"), 1, 0),
    
    satisfac.vivienda     = as.numeric(HC080),
    satisfac.vivienda_n   = ifelse(any(HC080 %in% c("1","2")), 0,
                                   ifelse(any(HC080 %in% c("3","4")), 1, NA_real_)),
    
    enough.temp.wint = ifelse(any(HC060 == "1"), 1, 0),
    enough.temp.summ = ifelse(any(HC070 == "1"), 1, 0),
    
    
    highlyurbanised  = ifelse(any(DB100 == 1), 1, 0),
    mediumpopulation = ifelse(any(DB100 == 2), 1, 0),
    lowpopulation    = ifelse(any(DB100 == 3), 1, 0),
    
    hhmembers        = as.numeric(HB120),
    
    detachedindep    = ifelse(any(HH010 == 1), 1, 0), 
    semidetached     = ifelse(any(HH010 == 2), 1, 0),
    flatless10       = ifelse(any(HH010 == 3), 1, 0),
    flatmore10       = ifelse(any(HH010 == 4), 1, 0),
    
    detached         = ifelse(any(HH010 %in% c(1,2)), 1, 0),
    flatorapt        = ifelse(any(HH010 %in% c(3,4)), 1, 0),
    
    ownership        = ifelse(any(HH021 == 1 | HH021 == 2), 1, 0),
    rental           = ifelse(any(HH021 == 3 | HH021 == 4), 1, 0),
    otherregimes     = ifelse(any(HH021 == 5), 1, 0),
        
    rooms            = as.numeric(HH030),
    
    oneadult         = ifelse(any(HX060 %in% c(1, 2, 3, 4, 5, 6, 10)), 1, 0),
    twoadults        = ifelse(any(HX060 %in% c(7, 8, 11, 12, 13)), 1, 0),
    otherhhtypes     = ifelse(any(HX060 %in% c(9, 14)), 1, 0),
    kids             = ifelse(any(HX060 %in% c(10, 11, 12, 13, 14)), 1, 0),
        
    sqmeters         = as.numeric(HC020),
    sqmeteroom       = as.numeric(HC020) / as.numeric(HH030),
     
    moisturedamage   = ifelse(any(HH040 == "1"), 1, 0),
     
    districtheating  = ifelse(any(HC001 == "1"), 1, 0),
    centralheating   = ifelse(any(HC001 == "2"), 1, 0),
    indivheating     = ifelse(any(HC001 == "3"), 1, 0),
    portableheating  = ifelse(any(HC001 == "4"), 1, 0),
    noheatingsystem  = ifelse(any(HC001 == "5"), 1, 0),
    
    electricity_heat = ifelse(any(HC002 == "1"), 1, 0),
    gas_heat         = ifelse(any(HC002 == "2"), 1, 0),
    gasoil_heat      = ifelse(any(HC002 == "3"), 1, 0),
    biomass_heat     = ifelse(any(HC002 == "4"), 1, 0),
    wood_heat        = ifelse(any(HC002 == "5"), 1, 0),
    coal_heat        = ifelse(any(HC002 == "6"), 1, 0),
    renewable_heat   = ifelse(any(HC002 == "7"), 1, 0),
    other_heat       = ifelse(any(HC002 == "8"), 1, 0),
    
    vhRentaa         = as.numeric(vhRentaa),
    
    actualrent       = as.numeric(HH060),
    housingcosts     = as.numeric(HH070),
    
    # Recode HC003 values to get the count variable "count"
    count        = case_when(
      HC003 == 1  ~ 3,
      HC003 == 2  ~ 2,
      HC003 == 3  ~ 1,
      HC003 == 4  ~ 0,
      HC003 == 99 ~ NA_real_,   
      TRUE        ~ NA_real_
    ), 
    
    # Recode HC003 values to get the binary "retrofit"
    retrofit = ifelse(HC003 %in% c(1, 2, 3), 1, 
                      ifelse(HC003 == 99, NA_real_, 0))
  )

individual <- individual %>%
  dplyr::select(
    RB030,  # Individual identifier (DB030 + two individual digits)
    RB050,  # Weights
    RB082,  # Age at the time of the interview
    RB090,  # Sex
    RB211,  # Situation in the activity in the actuality
    PL141,  # Type of contract (labor)
    PE041   # Level of education completed
  ) %>%
  
  # Create household identifier from "RB030" by removing the last two digits
  dplyr::mutate(
    DB030 = substr(RB030, 1, nchar(RB030) - 2),  # Extract household identifier
    DB030 = as.numeric(DB030) # Convert to numeric to remove leading zeros
  ) %>%
  
  # Group by household identifier
  group_by(DB030) %>%
  
  # Create household-level variables inside summarize()
  dplyr::summarize(
    # Weight: 
    weight          = first(as.numeric(RB050)),
    
    # Mean age of adults (age >= 18)
    meanageadults = ifelse(sum(RB082 >= 18, na.rm = TRUE) == 0, NA_real_,
                           mean(RB082[RB082 >= 18], na.rm = TRUE)),
    
    # Proportion of male adults (age >= 18)
    propmale        = ifelse(sum(RB082 >= 18, na.rm = TRUE) == 0, NA_real_,
                             sum(RB090 == 1 & RB082 >= 18, na.rm = TRUE) / sum(RB082 >= 18, na.rm = TRUE)),
    
    
    # Income stability: 1 if any member has RB211 == 3 or 4, or PL141 == 21 or 22
    incomestability  = ifelse(any(RB211 %in% c(3, 4), na.rm = TRUE) |
                                any(PL141 %in% c(21, 22), na.rm = TRUE), 1, 0),
    
    # All members in the household have high education (PE041 == 500)
    allhigheducation = ifelse(all(PE041 == 500, na.rm = TRUE), 1, 0)
  ) %>%
  
  # Ungroup after summarizing to get a tidy dataframe
  ungroup()

## e) Final merge ----
data <- left_join(household, individual, by = c("DB030", "weight"))

# Store original sample size BEFORE any cleaning
n_original_unweighted <- nrow(data)
n_original_weighted   <- sum(data$weight, na.rm = TRUE)

# ================== 2. LOAD & MERGE ADDITIONAL DATA ===========================

# Merge HDD & CDD by NUTS2 from AGRI4CAST JCR EU Portal ----

hdd_cdd <- "HDD_CDD.xlsx" # Load file downloaded from Agri4Cast 

HDD <- read_excel(hdd_cdd, sheet = 2) # Create HDD object (sheet 2 .xlsx)
names(HDD) # Check names
HDD <- HDD[, !names(HDD) %in% "CCAA"] # Remove unnecessary info ("CCAA")

CDD <- read_excel(hdd_cdd, sheet = 3) # Create CDD object (sheet 2 .xlsx)
names(CDD) # Check names
CDD <- CDD[, !names(CDD) %in% "CCAA"] # Remove unnecessary info ("CCAA")

names(HDD) # Check removal
names(CDD) # Check removal

# Merge "data" and "HDD" & "CDD" objects
data <- left_join(data, HDD, by = "NUTS2")
data <- left_join(data, CDD, by = "NUTS2")

# Keep variables of interest
data$HDD_1723 <- data$HDD_media_17_23
data$CDD_1723 <- data$CDD_media_17_23

remove("hdd_cdd") # Clean environment & workspace
gc() # Keep memory usage low

# ================ 3. PRELIMINARY TRANSFORMATIONS ==============================

## a. Create dummies & factor(region) ----

data$region <- NA
data$region[data$galicia        == 1] <- "galicia"
data$region[data$asturias       == 1] <- "asturias"
data$region[data$cantabria      == 1] <- "cantabria"
data$region[data$paisvasco      == 1] <- "paisvasco"
data$region[data$navarra        == 1] <- "navarra"
data$region[data$larioja        == 1] <- "larioja"
data$region[data$aragon         == 1] <- "aragon"
data$region[data$madrid         == 1] <- "madrid"
data$region[data$castillaleon   == 1] <- "castillaleon"
data$region[data$castillamancha == 1] <- "castillamancha"
data$region[data$extremadura    == 1] <- "extremadura"
data$region[data$catalunya      == 1] <- "catalunya"
data$region[data$comvalenciana  == 1] <- "comvalenciana"
data$region[data$illesbalears   == 1] <- "illesbalears"
data$region[data$murcia         == 1] <- "murcia"
data$region[data$ceuta          == 1] <- "ceuta"
data$region[data$melilla        == 1] <- "melilla"
data$region[data$andalucia      == 1] <- "andalucia"
data$region[data$canarias       == 1] <- "canarias"
data$region <- factor(data$region)

## b. Identify & remove NA's (retrofit) ----

# Count NA's in retrofit before removal
n_na_retrofit              <- sum(is.na(data$retrofit))
weight_na_retrofit         <- sum(data$weight[is.na(data$retrofit)], na.rm = TRUE)
pct_na_retrofit_unweighted <- 100 * n_na_retrofit / nrow(data)
pct_na_retrofit_weighted   <- 100 * weight_na_retrofit / sum(data$weight, na.rm = TRUE)

# Control accumulated % of removed observations (real and in the dataset)
removed_dataset  <- pct_na_retrofit_unweighted
removed_weighted <- pct_na_retrofit_weighted

data <- data[!is.na(data$retrofit),] # Remove NA's

## c. Clean "income" ----

boxplot(data$vhRentaa)
hist(data$vhRentaa)
summary(data$vhRentaa)

# First, we identify zero/negative income values --> Exclude if they are low % 
# Limit fixed at 100€ and not in 0€ because there are 22 observations between
# 0 and 100 which we assume are not realistic

# Count observations to be removed
n_lowincome_partial <- length(data$vhRentaa[data$vhRentaa <= 100 & data$vhRentaa > 0])
n_lowincome_total   <- length(data$vhRentaa[data$vhRentaa <= 100])
n_negative_income   <- length(data$vhRentaa[data$vhRentaa < 0])

pct_lowincome_unweighted <- 100 * n_lowincome_total / nrow(data)
weight_lowincome         <- sum(data$weight[data$vhRentaa <= 100], na.rm = TRUE)
pct_lowincome_weighted   <- 100 * weight_lowincome / sum(data$weight, na.rm = TRUE)

# Control accumulated % of removed observations (real and in the dataset)
removed_dataset  <- removed_dataset  + pct_lowincome_unweighted
removed_weighted <- removed_weighted + pct_lowincome_weighted

# Exclude the observations (to then log-transform & scale)
data <- data[data$vhRentaa > 100,]

# Control removed observations (real and in the dataset)
nrow(data)
removed_dataset  
removed_weighted 


# ================== 4. TRANSFORM VARIABLES ==================================

## a. Log & scale "income" ----

# Take logs vhRentaa to normalize it (right-skewed) & reduce asymmetry
data$ln_income <- log(data$vhRentaa)

# Check distribution of log(income)
summary(data$ln_income)
hist(data$ln_income)
boxplot(data$ln_income)

# ══════════════════════════════════════════════════════════════════════════════
# 5. TRACKING: SAMPLE SELECTION & REMOVAL SUMMARY
# ══════════════════════════════════════════════════════════════════════════════

# Store sample selection information for next script
sample_info <- list(
  n_original_unweighted           = n_original_unweighted,
  n_original_weighted             = n_original_weighted,
  
  # Step 1: Remove NA in retrofit
  n_removed_retrofit_unweighted   = n_na_retrofit,
  n_removed_retrofit_weighted     = weight_na_retrofit,
  pct_removed_retrofit_unweighted = pct_na_retrofit_unweighted,
  pct_removed_retrofit_weighted   = pct_na_retrofit_weighted,
  
  # Step 2: Remove low/negative income (<=100)
  n_removed_income_unweighted      = n_lowincome_total,
  n_removed_income_weighted_before = weight_lowincome,
  pct_removed_income_unweighted    = pct_lowincome_unweighted,
  pct_removed_income_weighted      = pct_lowincome_weighted,
  
  # Total removed after Step 1 + 2
  n_after_step2_unweighted          = nrow(data),
  n_after_step2_weighted            = sum(data$weight, na.rm = TRUE),
  pct_removed_cumulative_unweighted = removed_dataset,
  pct_removed_cumulative_weighted   = removed_weighted,
  
  # Final sample (will be updated in next script after listwise deletion)
  n_final_unweighted = nrow(data),
  n_final_weighted   = sum(data$weight, na.rm = TRUE)
)

# Save data + sample_info
save(data, sample_info, file = "data.RData")

# Print summary
cat(sprintf("Original sample (unweighted):        n = %s\n", 
            format(sample_info$n_original_unweighted, big.mark = ",")))
cat(sprintf("Original sample (weighted):          N = %s\n", 
            format(round(sample_info$n_original_weighted), big.mark = ",")))
cat("\n")
cat(sprintf("Removed (retrofit NA, unweighted):   n = %d (%.2f%%)\n",
            sample_info$n_removed_retrofit_unweighted,
            sample_info$pct_removed_retrofit_unweighted))
cat(sprintf("Removed (retrofit NA, weighted):     N = %.0f (%.2f%%)\n",
            sample_info$n_removed_retrofit_weighted,
            sample_info$pct_removed_retrofit_weighted))
cat("\n")
cat(sprintf("Removed (income ≤100, unweighted):   n = %d (%.2f%%)\n",
            sample_info$n_removed_income_unweighted,
            sample_info$pct_removed_income_unweighted))
cat(sprintf("Removed (income ≤100, weighted):     N = %.0f (%.2f%%)\n",
            sample_info$n_removed_income_weighted_before,
            sample_info$pct_removed_income_weighted))
cat("\n")
cat(sprintf("After cleaning (unweighted):         n = %s (%.2f%% retained)\n",
            format(sample_info$n_after_step2_unweighted, big.mark = ","),
            100 - sample_info$pct_removed_cumulative_unweighted))
cat(sprintf("After cleaning (weighted):           N = %s (%.2f%% retained)\n",
            format(round(sample_info$n_after_step2_weighted), big.mark = ","),
            100 - sample_info$pct_removed_cumulative_weighted))


# Track and print total time of execution
end.time <- Sys.time()
cat("Total time of execution:", round(end.time - start.time, 2), "min\n")
