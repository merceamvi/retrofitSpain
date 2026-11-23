###############################################################################
# Project       : [rehabSpain] [0] SHARED SETUP (packages, theme, helpers)
# Creation date : 02/12/2024
# Last update   : 23/11/2026 (stage numbering 02-05)
# Author        : Mercè Amich (merce.amich@ehu.eus)
# Institution   : UPV/EHU, BC3
# Script Overview:
#   Common preamble sourced by scripts 02-05: package loading, working
#   directory, output folders, publication ggplot2 theme, color palettes and
#   save helpers.
#   Folders:
#     paper/   every table and figure that appears in the manuscript, named by
#              its number in the paper (Table_3.html, Figure_B1.pdf, ...)
#     output/  everything else: diagnostics, checkpoints and the numbers that
#              the text cites without a table (CSV / RDS / TXT)
#   Not a pipeline stage on its own; it is sourced at the top of each
#   analysis script so they share an identical environment.
# Requirements:
#   - source("_setup.R") at the top of scripts 02 to 05
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

# Random number generation: R's defaults since 3.6.0, set explicitly so that
# every seeded step (multi-start EM, bootstrap resamples, matching draws) gives
# the same result on any R version or configuration.
RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

# Working directory
path <- here()
setwd(path)

# Output folders. Their paths are stored as options so the save helpers can
# reach them from any stage.
dir_paper  <- file.path(path, "paper")
dir_output <- file.path(path, "output")
for (d in c(dir_paper, dir_output)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)
options(rehab.dir_paper = dir_paper, rehab.dir_output = dir_output)

# Times New Roman on the Windows screen and TIFF devices. Without this, every
# plot with base_family = "Times New Roman" raises "font family not found".
# On other systems, map the family to the PDF device's Times so that plots
# printed to a PDF screen device (Rscript) do not fail.
if (.Platform$OS.type == "windows") {
  grDevices::windowsFonts("Times New Roman" = grDevices::windowsFont("Times New Roman"))
} else if (!"Times New Roman" %in% names(grDevices::pdfFonts())) {
  grDevices::pdfFonts("Times New Roman" = grDevices::pdfFonts()$Times)
}

# Record the R and package versions of the run: exact reproduction of the
# seeded results requires the same versions (notably mclust and FactoMineR).
writeLines(c(format(Sys.time()), capture.output(sessionInfo()), "",
             paste("RNGkind:", paste(RNGkind(), collapse = ", "))),
           file.path(dir_output, "sessionInfo.txt"))

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
# Paper objects (tables and figures of the manuscript) go to paper/
paper_path  <- function(fname) file.path(getOption("rehab.dir_paper"),  fname)
# Everything else goes to output/
output_path <- function(fname) file.path(getOption("rehab.dir_output"), fname)

save_table <- function(kbl, fname) {
  kableExtra::save_kable(kbl, file = paper_path(fname))
  cat(sprintf("  Saved: paper/%s\n", fname))
}

# Paper figure in PDF (cairo) and TIFF, same file stem
save_figure <- function(p, stem, w_pdf, h_pdf, w_tiff, h_tiff, dpi = 900) {
  ggsave(paper_path(paste0(stem, ".pdf")), p,
         width = w_pdf, height = h_pdf, device = cairo_pdf)
  ggsave(paper_path(paste0(stem, ".tiff")), plot = p, device = "tiff",
         width = w_tiff, height = h_tiff, units = "cm", dpi = dpi)
  cat(sprintf("  Saved: paper/%s.pdf and .tiff\n", stem))
}

# Non-paper tables (diagnostics, numbers cited in the text)
save_output <- function(df, fname) {
  write.csv(df, output_path(fname), row.names = FALSE)
  cat(sprintf("  Saved: output/%s\n", fname))
}

# Canonic order of variables in all of the results tables
# Canonical variable order, shared by Tables 3, 4, A.1, B.1, B.2, D.1 and D.2.
# Rows are grouped in six blocks (climate as one block, with heating and cooling
# as two rows); the Shapley attribution in script 05 splits climate into two
# blocks, heating and cooling, and therefore uses seven.
var_order <- c("ownership",                                   # A. Tenure
               "ln_income",                                   # B. Income
               "allhigheducation",                            # C. Education
               "HDD_1723", "CDD_1723",                        # D. Climate
               "rooms", "detached", "moisturedamage", "pollution",  # E. Dwelling
               "meanageadults", "kids", "highlyurbanised")    # F. Household

# Block membership of each variable, used to group the rows of the results
# tables. Labels match the six blocks of Section 3.
block_map <- c(
  "ownership"        = "A. Tenure",
  "ln_income"        = "B. Income",
  "allhigheducation" = "C. Education",
  "HDD_1723"         = "D. Climatic exposure",
  "CDD_1723"         = "D. Climatic exposure",
  "rooms"            = "E. Dwelling characteristics",
  "detached"         = "E. Dwelling characteristics",
  "moisturedamage"   = "E. Dwelling characteristics",
  "pollution"        = "E. Dwelling characteristics",
  "meanageadults"    = "F. Household composition and location",
  "kids"             = "F. Household composition and location",
  "highlyurbanised"  = "F. Household composition and location"
)

block_order <- c("A. Tenure",
                 "B. Income",
                 "C. Education",
                 "D. Climatic exposure",
                 "E. Dwelling characteristics",
                 "F. Household composition and location")

# Short block labels, for figure facets where the full names do not fit
block_map_short <- c(
  "ownership"        = "Tenure",
  "ln_income"        = "Income",
  "allhigheducation" = "Education",
  "HDD_1723"         = "Climate",
  "CDD_1723"         = "Climate",
  "rooms"            = "Dwelling",
  "detached"         = "Dwelling",
  "moisturedamage"   = "Dwelling",
  "pollution"        = "Dwelling",
  "meanageadults"    = "Household",
  "kids"             = "Household",
  "highlyurbanised"  = "Household"
)

block_order_short <- c("Tenure", "Income", "Education",
                       "Climate", "Dwelling", "Household")
