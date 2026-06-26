###############################################################################
# Project       : [rehabSpain] [0] SHARED SETUP (packages, theme, helpers)
# Creation date : 02/12/2024
# Last update   : 26/06/2026
# Author        : Mercè Amich (merce.amich@ehu.eus)
# Institution   : EHU, BC3
# Script Overview:
#   Common preamble sourced by scripts 02-04: package loading, working
#   directory, publication ggplot2 theme, color palettes and save helpers.
#   Not a pipeline stage on its own; it is sourced at the top of each
#   analysis script so they share an identical environment.
# Requirements:
#   - source("_setup.R") at the top of scripts 02, 03 and 04
# Output  : none (defines objects in the calling environment)
###############################################################################

Sys.setenv(LANG = "en")

# Install and load all required packages
packages_needed <- c(
  "here",                                          # wd handling
  "dplyr", "tidyr", "tibble", "scales", "Hmisc",   # Data wrangling
  "boot", "lmtest", "car", "AER", "broom",         # Modelling & inference
  "FactoMineR", "factoextra", "cluster", "mclust", # FAMD & clustering
  "ggplot2", "ggrepel", "patchwork", "forcats",    # Visualisation
  "kableExtra", "stringr",                         # Publication 
  "weights",                                       # Weighted correlations 
  "parallel", "foreach", "doParallel"              # Bootstrap parallel processing
)

for (p in packages_needed) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

# Working directory (tables and figures saved here) 
path <- here()
setwd(path)

# Set ggplot2 theme for publication-quality figures
theme_paper <- function(base_size = 11) {
  theme_minimal(base_size = base_size, base_family = "Times New Roman") +
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 1,
                                      hjust = 0, margin = margin(b = 5)),
      plot.subtitle    = element_text(size = base_size - 1, color = "grey40",
                                      hjust = 0, margin = margin(b = 8)),
      plot.caption     = element_text(size = base_size - 2, color = "grey40",
                                      hjust = 0, margin = margin(t = 6)),
      axis.title       = element_text(size = base_size),
      axis.text        = element_text(size = base_size - 1, color = "black"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
      strip.text       = element_text(face = "bold", size = base_size),
      legend.title     = element_text(face = "bold", size = base_size - 1),
      legend.text      = element_text(size = base_size - 1),
      legend.position  = "bottom"
    )
}

# Set color palettes
pal_profiles <- c("#e74c3c", "#3498db", "#2ecc71", "#f39c12", "#9b59b6",
                  "#d6a96a", "#b5703a", "#7b4b21")

# Helper functions for saving outputs
save_plot <- function(p, fname, w = 20, h = 8) {
  fp <- file.path(getwd(), fname)
  ggsave(fp, plot = p, device = "tiff", 
         width = w, height = h, units = "cm", dpi = 300)
  cat(sprintf("  Saved: %s\n", fname))
}

save_table <- function(kbl, fname) {
  fp <- file.path(getwd(), fname)
  kableExtra::save_kable(kbl, file = fp)
  cat(sprintf("  Saved: %s\n", fname))
}
