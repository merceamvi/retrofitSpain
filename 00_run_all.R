###############################################################################
# Project       : [rehabSpain] [00] MASTER RUNNER
# Creation date : 09/06/2026
# Last update   : 23/11/2026
# Author        : Merce Amich (merce.amich@ehu.eus)
# Institution   : UPV/EHU, BC3

# Approximate run time : ~6 h in total
#   - 01_build_data.R          ~2 min  (INE download; needs internet)
#   - 02_descriptives_logit.R  ~11 min
#   - 03_famd_clustering.R     ~6 h    (multi-start GMM + bootstrap stability)
#   - 04_profiles_ames.R       ~5 min
#   - 05_decomposition.R       ~16 min (N_BOOT_DECOMP, default 500)

# Script Overview:
#   Runs the full pipeline in order. Each stage saves its intermediate objects
#   to disk, so stages can also be run individually. Stage 03 numbers the
#   profiles by weighted adoption (P1 highest), and P1 is the reference profile
#   of the decomposition in stage 05.
#   Tables and figures of the manuscript are written to paper/, named by their
#   number in the paper; diagnostics and the numbers cited in the text without
#   a table are written to output/.

# Requirements:
#   - All stage scripts and _setup.R in the same directory
#   - HDD_CDD.xlsx (for stage 01)
#   - Internet access for stage 01 (INE download); comment it out if
#     data.RData is already present in the directory.

# Output  : data.RData + intermediate .RData + paper/ + output/
###############################################################################

rm(list = ls(all = TRUE))
Sys.setenv(LANG = "en")

source("01_build_data.R")
source("02_descriptives_logit.R")
source("03_famd_clustering.R")
source("04_profiles_ames.R")
source("05_decomposition.R")

cat("PIPELINE COMPLETE\n")
