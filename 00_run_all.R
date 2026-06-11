###############################################################################
# Project       : [rehabSpain] [00] MASTER RUNNER
# Creation date : 09/06/2026
# Last update   : 09/06/2026
# Author        : Mercè Amich (merce.amich@ehu.eus)
# Institution   : UPV/EHU, BC3
# Last run time : —
# Script Overview:
#   Runs the full pipeline in order. Each stage saves its intermediate
#   objects to disk, so stages can also be run individually.
# Requirements:
#   - All stage scripts and _setup.R in the same directory
#   - HDD_CDD.xlsx (for stage 01)
#   - Internet access for stage 01 (INE download); skip 01 if data.RData
#     is already present in the directory.
# Output  : data.RData + intermediate .RData + publication tables/figures
###############################################################################

rm(list = ls(all = TRUE))
Sys.setenv(LANG = "en")

# Stage 01 downloads and builds data.RData from INE microdata (needs internet).
# Comment it out if data.RData is already in the directory.
source("01_build_data.R")

# Analysis stages (each loads the previous stage's intermediate objects)
source("02_descriptives_logit.R")
source("03_famd_clustering.R")
source("04_profiles_policy.R")

cat("\n=== PIPELINE COMPLETE ===\n")
