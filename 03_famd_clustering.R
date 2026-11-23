###############################################################################
# Project       : [rehabSpain] [3] FAMD & GMM CLUSTERING
# Creation date : 02/12/2024
# Last update   : 23/11/2026 (profiles numbered by weighted adoption, P1 highest;
#                 short profile labels saved for stage 04; dimension retention by
#                 the equal-share benchmark; Tables C.1, C.2 and 3 as in the paper;
#                 English object names; selection by highest Jaccard, as in 4.2)
# Author        : Mercè Amich (merce.amich@ehu.eus)
# Institution   : UPV/EHU, BC3
# Last run time : ~6 h (bootstrap of 3 x 4 candidates)

# Script Overview:
#   Factor Analysis of Mixed Data (FAMD), multi-criteria assessment of GMM-VEV
#   solutions over candidate dimensions (multi-start EM reference fits,
#   bootstrap cluster-wise Jaccard with the FAMD re-estimated in each resample),
#   selection rule, diagnostics of the five-dimensional space, final GMM and
#   household-profile characterisation. Profiles are numbered by weighted
#   adoption, from the highest (P1) to the lowest (P4).

# Requirements:
#   - This script must be in the same directory as:
#       a) _setup.R
#       b) 02_descriptives_logit.RData   (output of 02_descriptives_logit.R)

# Output  : paper/  Table_C1, Table_C2, Table_C3, Table_3 (HTML), Figure_C1
#           output/ assessment of all candidates, cluster-wise Jaccard,
#                   dimension-5 diagnostics
#           + 03_famd_clustering.RData (objects for the next stage)

###############################################################################

# ══════════════════════════════════════════════════════════════════════════════
# 0. SETUP
# ══════════════════════════════════════════════════════════════════════════════

rm(list = ls(all = TRUE))
Sys.setenv(LANG = "en")
start_time <- Sys.time()

# Shared packages, working directory, theme, palettes and save helpers
source("_setup.R")

# Objects from the previous stage:
#   data_clean, logit_baseline, formula_baseline, ame_df,
#   active_vars, var_labels, outcome_var
load("02_descriptives_logit.RData")

# ══════════════════════════════════════════════════════════════════════════════
# FACTOR ANALYSIS OF MIXED DATA (FAMD)
# ══════════════════════════════════════════════════════════════════════════════

# Prepare data: factorise binary variables for FAMD
data_famd <- data_clean %>%
  select(all_of(c(active_vars, outcome_var))) %>%
  mutate(
    across(all_of(active_vars),
           ~ if (is.numeric(.) && all(na.omit(.) %in% 0:1))
             factor(., levels = 0:1, labels = c("No", "Yes")) else .),
    across(all_of(outcome_var),
           ~ factor(., levels = 0:1, labels = c("No", "Yes")))
  )

# Run FAMD with retrofit as supplementary variable
idx_sup <- match(outcome_var, names(data_famd))

set.seed(123) # Set seed for reproducibility
famd_result <- FactoMineR::FAMD(data_famd,
                                ncp     = 12, 
                                sup.var = idx_sup,
                                row.w   = data_clean$weight_final,
                                graph   = FALSE)

# Extract eigenvalues
eig_df           <- as.data.frame(famd_result$eig)
colnames(eig_df) <- c("Eigenvalue", "Variance_%", "Cumulative_%")
eig_df$Dimension <- seq_len(nrow(eig_df))

# Dimensions above the equal-share benchmark (100 / 12 = 8.33% of variance)
n_factors <- sum(eig_df$`Variance_%` > 100 / length(active_vars))

cat(sprintf("  Equal-share benchmark: %d dimensions above it (%.1f%% cumulative variance)\n",
            n_factors, eig_df$`Cumulative_%`[n_factors]))

# Candidate factor spaces carried forward to the clustering stage
dims <- c(3, 4, 5)

# TABLE C1: Eigenvalues (candidate dimensionalities in bold)
kblC1 <- eig_df %>%
  filter(Dimension <= 7) %>%
  mutate(across(c(Eigenvalue, `Variance_%`, `Cumulative_%`), ~ round(., 2))) %>%
  select(Dimension, Eigenvalue, `Variance_%`, `Cumulative_%`) %>%
  kbl(format = "html",
      caption = "Table C.1. FAMD dimension retention",
      col.names = c("Dim.", "Eigenvalue", "Var. (%)", "Cumul. (%)"),
      align = c("c", "r", "r", "r"),
      row.names = FALSE) %>%
  kable_classic(full_width = FALSE, html_font = "Times New Roman") %>%
  kable_styling(font_size = 11) %>%
  row_spec(0, bold = TRUE) %>%
  row_spec(dims, bold = TRUE) %>%
  footnote(general = paste("The first five dimensions exceed the 8.33% equal-share",
                           "benchmark; the candidate factor spaces retain three, four",
                           "or five."),
           general_title = "Notes:", footnote_as_chunk = TRUE)

print(kblC1)
save_table(kblC1, "Table_C1.html")

# FIGURE C1: Scree plot
figC1 <- ggplot(eig_df %>% filter(Dimension <= 8),
                aes(x = Dimension, y = `Variance_%`)) +
  geom_col(aes(fill = `Variance_%` > 100 / length(active_vars)),
           width = 0.65, color = "gray20", linewidth = 0.3) +
  geom_line(linewidth = 0.6, color = "gray40") +
  geom_point(size = 2.5, color = "black", fill = "white", shape = 21, stroke = 0.5) +
  geom_hline(yintercept = 100 / length(active_vars),
             linetype = "dashed", color = "gray40", linewidth = 0.4) +
  scale_fill_manual(values = c("TRUE" = "gray50", "FALSE" = "gray85"),
                    labels = c("TRUE" = "Carried forward", "FALSE" = "Excluded"),
                    name = NULL) +
  scale_x_continuous(breaks = 1:8) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title    = NULL,
       x        = "Dimension", 
       y        = "Variance explained (%)") +
  theme_minimal(base_size = 10) +
  theme(
    plot.caption          = element_text(size = 8, color = "gray30",
                                         hjust = 0, margin = margin(t = 8)),
    plot.caption.position = "plot",
    axis.text             = element_text(size = 9, color = "black"),
    axis.title.x          = element_text(size = 10, margin = margin(t = 6)),
    axis.title.y          = element_text(size = 10, margin = margin(r = 6)),
    panel.background      = element_blank(),
    panel.grid.major.x    = element_blank(),
    panel.grid.major.y    = element_line(color = "gray85", linewidth = 0.3),
    panel.grid.minor      = element_blank(),
    axis.ticks            = element_line(color = "gray40", linewidth = 0.3),
    axis.ticks.length     = unit(2, "pt"),
    legend.position       = c(0.95, 0.95),
    legend.justification  = c("right", "top"),
    legend.key.size       = unit(0.4, "cm"),
    legend.text           = element_text(size = 8),
    legend.background     = element_rect(fill = "white", color = "gray70", linewidth = 0.3),
    legend.margin         = margin(3, 4, 3, 4),
    plot.margin           = margin(t = 5, r = 10, b = 5, l = 5)
  )

print(figC1)
save_figure(figC1, "Figure_C1",
            w_pdf = 6.5, h_pdf = 3.5, w_tiff = 14, h_tiff = 9)

# Extract FAMD components
var_contrib   <- famd_result$var$contrib    [, 1:n_factors, drop = FALSE]
var_coord     <- famd_result$var$coord      [, 1:n_factors, drop = FALSE]
sup_coord     <- famd_result$quali.sup$coord[, 1:n_factors, drop = FALSE]
factor_scores <- famd_result$ind$coord      [, 1:n_factors, drop = FALSE]

# Extract retrofit (== 1) coordinates
ret_row <- grepl("^Yes$", rownames(sup_coord))
ret_yes <- as.numeric(sup_coord[ret_row, ])

# TABLE C.2: Dimension interpretation
# Variable labels as in Table A.1 for the climate variables, var_labels otherwise
labels_c2 <- var_labels
labels_c2["HDD_1723"] <- "HDD (mean 2017-2023)"
labels_c2["CDD_1723"] <- "CDD (mean 2017-2023)"

dim2 <- data.frame(
  Dim = 1:n_factors,
  `Var (%)` = round(eig_df$`Variance_%`[1:n_factors], 2),
  `Top 3`   = sapply(1:n_factors, function(d) {
    cc <- var_contrib[, d, drop = TRUE]
    top <- names(sort(cc[!is.na(cc)], decreasing = TRUE)[1:3])
    paste(labels_c2[top], collapse = ", ")
  }),
  `|Gs(·)|>0.4` = sapply(1:n_factors, function(d) {
    cl <- var_coord[, d, drop = TRUE]
    cl <- cl[!is.na(cl)]
    s  <- names(cl[abs(cl) > 0.4])
    if (length(s) == 0) "-" else paste(labels_c2[s], collapse = ", ")
  }),
  `Retrofit` = round(ret_yes[1:n_factors], 2),
  check.names = FALSE
)
kblC2 <- kbl(dim2, format = "html",
             caption = "Table C.2. FAMD dimension interpretation",
             align = c("c", "r", "l", "l", "r")) %>%
  kable_classic(full_width = FALSE, html_font = "Times New Roman") %>%
  kable_styling(font_size = 11) %>%
  row_spec(0, bold = TRUE) %>%
  footnote(general = "Retrofit = G_s(·) of 'Adopt = Yes'.",
           general_title = "Notes:", footnote_as_chunk = TRUE)

print(kblC2)
save_table(kblC2, "Table_C2.html")

# Clean up
remove(eig_df, var_contrib, var_coord, sup_coord, ret_row, ret_yes,
       dim2, kblC2, figC1, idx_sup, labels_c2)
gc()


# ══════════════════════════════════════════════════════════════════════════════
# MULTI-CRITERIA ASSESSMENT OF GMM-VEV SOLUTIONS (Section 4.2, Table C.3)
# ══════════════════════════════════════════════════════════════════════════════
# Design (revised 09/2026):
#   1. ONE reference fit per (FAMD dimensions, G), by MULTI-START EM (n_starts
#      random initialisations, highest log-likelihood kept; Biernacki et al.,
#      2003). With n > 2000, mclust initialises EM on a random 2000-row
#      subsample, so single fits depend on the seed. The seeds (1000 + s) are
#      the ones that produced the published partition and must not change.
#   2. Every diagnostic of a row is computed from that reference fit, and the
#      final model is that same object, so Table 3 and Table C.3 agree.
#   3. Bootstrap stability: cluster-wise Jaccard (Hennig, 2007) over B = 200
#      resamples. In each resample the FAMD is RE-ESTIMATED with survey weights
#      and the mixture re-fitted, so the coefficient reflects the whole
#      procedure. Each reference component is compared with the most similar
#      bootstrap component: J_g = max_h |C_g n C_h| / (|C_g| + |C_h| - |C_g n C_h|).
#      (The previous version used max(cont[, h]) instead of the column sum in
#      the denominator, which inflated the coefficient; corrected here.)
#   4. The same resample indices are used for every G within a dimension
#      (paired comparisons).
#   5. Ward clustering is no longer part of the procedure.
#   6. The mixture is fitted to the factor scores of the WEIGHTED FAMD without
#      further weighting (mclust takes no weights); shares and adoption rates
#      are weighted when the profiles are characterised.

# Candidate retained FAMD dimensions: dims, defined before Table C.1
G_values <- 3:6          # Candidate numbers of mixture components
n_boot   <- 200          # Bootstrap resamples for the Jaccard coefficient
n_starts <- 20           # Multi-start EM initialisations per reference fit

sup_idx  <- match(outcome_var, names(data_famd))

# ── Helper: multi-start Mclust fit, best log-likelihood kept ─────────────────
fit_gmm_multistart <- function(X, G, modelNames = "VEV",
                               n_starts = 20, seed_base = 1000) {
  ll   <- rep(NA_real_, n_starts)
  best <- NULL
  for (s in seq_len(n_starts)) {
    set.seed(seed_base + s)   # seeds the random init subsample of this start
    fit <- tryCatch(
      Mclust(X, G = G, modelNames = modelNames, verbose = FALSE),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    ll[s] <- fit$loglik
    if (is.null(best) || fit$loglik > best$loglik) best <- fit
  }
  if (is.null(best)) return(NULL)
  ll_ok <- ll[!is.na(ll)]
  # Tolerance for "same optimum": EM stops on a relative loglik change of
  # ~1e-5, so differences below that scale are convergence noise.
  tol   <- max(1, 1e-5 * abs(max(ll_ok)))
  list(fit        = best,
       loglik_all = ll,
       n_optima   = sum(diff(sort(ll_ok)) > tol) + 1,
       share_best = mean(ll_ok >= max(ll_ok) - tol))
}

# ── Helper: cluster-wise Jaccard (Hennig, 2007) ──────────────────────────────
# For each reference component g: max over bootstrap components h of
# n_gh / (rowSum_g + colSum_h - n_gh). Returns one value per reference label.
jaccard_by_cluster <- function(ref_labels, boot_labels, ref_levels) {
  cont <- table(factor(ref_labels, levels = ref_levels), boot_labels)
  rs   <- rowSums(cont)
  cs   <- colSums(cont)
  vapply(seq_len(nrow(cont)), function(g) {
    if (rs[g] == 0) return(NA_real_)
    n_gh <- cont[g, ]
    max(n_gh / (rs[g] + cs - n_gh))
  }, numeric(1))
}

# ── Helper: degeneracy diagnostics of a fitted GMM ───────────────────────────
# Near-deterministic posteriors or near-singular component covariances signal
# components collapsing onto thin bands of the factor space, not separation.
degeneracy_flags <- function(fit) {
  sigma <- fit$parameters$variance$sigma
  cond  <- vapply(seq_len(dim(sigma)[3]), function(g) {
    ev <- eigen(sigma[, , g], symmetric = TRUE, only.values = TRUE)$values
    max(ev) / min(ev)
  }, numeric(1))
  data.frame(
    pro_min  = min(fit$parameters$pro),
    size_min = min(table(fit$classification)),
    cond_max = max(cond),
    near_det = as.integer(mean(1 - apply(fit$z, 1, max)) < 0.005)
  )
}

gmm_store   <- list()   # reference fits, "d<dim>_G<G>"
assess_list <- list()   # one row per (dim, G)
clusterwise <- list()   # cluster-wise Jaccard, one row per component

for (dim in dims) {
  
  cat(sprintf("\n=== FAMD with %d retained dimensions ===\n", dim)); flush.console()
  factor_scores <- famd_result$ind$coord[, 1:dim, drop = FALSE]
  dist_mat      <- dist(factor_scores, method = "euclidean")
  
  # ── Reference fits (multi-start) ───────────────────────────────────────────
  ref <- list()
  for (g in G_values) {
    ms <- fit_gmm_multistart(factor_scores, G = g, modelNames = "VEV",
                             n_starts = n_starts, seed_base = 1000)
    if (is.null(ms)) { cat(sprintf("  G = %d: all starts failed\n", g)); next }
    ref[[as.character(g)]] <- ms
    gmm_store[[sprintf("d%d_G%d", dim, g)]] <- ms$fit
    cat(sprintf("  G = %d: loglik = %.2f | distinct optima = %d | share at best = %.2f\n",
                g, ms$fit$loglik, ms$n_optima, ms$share_best))
  }
  
  # ── Bootstrap: FAMD re-estimated in each resample ──────────────────────────
  set.seed(123)
  idx_mat <- replicate(n_boot, sample(nrow(factor_scores), replace = TRUE))
  J <- lapply(G_values, function(g) matrix(NA_real_, n_boot, g))
  names(J) <- as.character(G_values)
  
  t0 <- Sys.time()
  for (b in seq_len(n_boot)) {
    idx_b <- idx_mat[, b]
    fb <- tryCatch(
      FactoMineR::FAMD(data_famd[idx_b, ],
                       ncp     = dim,
                       sup.var = sup_idx,
                       row.w   = data_clean$weight_final[idx_b],
                       graph   = FALSE),
      error = function(e) NULL)
    if (is.null(fb)) next
    Xb <- fb$ind$coord[, 1:dim, drop = FALSE]
    
    set.seed(20000 + b)   # single-start replicate fits, reproducible
    for (g in G_values) {
      gg <- as.character(g)
      if (is.null(ref[[gg]])) next
      gb <- tryCatch(Mclust(Xb, G = g, modelNames = "VEV", verbose = FALSE),
                     error = function(e) NULL)
      if (is.null(gb)) next
      J[[gg]][b, ] <- jaccard_by_cluster(ref[[gg]]$fit$classification[idx_b],
                                         gb$classification, seq_len(g))
    }
    rm(fb, Xb)
    if (b %% 25 == 0) {
      el <- as.numeric(Sys.time() - t0, units = "mins")
      cat(sprintf("    bootstrap %d/%d  (%.1f min, about %.1f min left)\n",
                  b, n_boot, el, el / b * (n_boot - b)))
      flush.console(); gc()
    }
  }
  
  # ── Diagnostics of each reference fit ──────────────────────────────────────
  for (g in G_values) {
    gg <- as.character(g)
    if (is.null(ref[[gg]])) next
    fit <- ref[[gg]]$fit
    z   <- fit$z; z[z < 1e-10] <- 1e-10
    Jm  <- J[[gg]]
    j_rep <- rowMeans(Jm, na.rm = TRUE)          # mean over components, per resample
    
    assess_list[[length(assess_list) + 1]] <- cbind(
      data.frame(
        Dimension   = dim,
        G           = g,
        Jaccard     = mean(j_rep, na.rm = TRUE),
        MCse        = sd(j_rep, na.rm = TRUE) / sqrt(sum(!is.na(j_rep))),
        Entropy     = mean(-rowSums(z * log(z))) / log(g),
        Uncertainty = mean(1 - apply(fit$z, 1, max)),
        Silhouette  = mean(silhouette(fit$classification, dist_mat)[, "sil_width"]),
        loglik      = fit$loglik,
        BIC         = fit$bic,
        ICL         = as.numeric(mclust::icl(fit)),
        n_optima    = ref[[gg]]$n_optima,
        share_best  = ref[[gg]]$share_best),
      degeneracy_flags(fit))
    
    clusterwise[[length(clusterwise) + 1]] <- data.frame(
      Dimension     = dim,
      G             = g,
      Component     = seq_len(g),
      Jaccard       = colMeans(Jm, na.rm = TRUE),
      Dissolved_pct = 100 * colMeans(Jm <  0.50, na.rm = TRUE),
      Recovered_pct = 100 * colMeans(Jm >= 0.75, na.rm = TRUE))
  }
  
  save(gmm_store, assess_list, clusterwise,
       file = output_path(sprintf("checkpoint_03_dim%d.RData", dim)))
  rm(dist_mat, idx_mat, J, ref); gc()
}

assessment <- do.call(rbind, assess_list)
clusterwise <- do.call(rbind, clusterwise)
rownames(assessment) <- NULL

cat("\n-- Assessment of every candidate --\n")
print(assessment, digits = 3)
save_output(assessment,  "C3_assessment_all_candidates.csv")
save_output(clusterwise, "C3_clusterwise_jaccard.csv")

# ── SELECTION RULE (Section 4.2) ──────────────────────────────────────────────
# Admissible: mean Jaccard >= 0.60 and no degenerate component (near-
# deterministic posteriors or a component covariance with condition number
# above 100). Among admissible candidates, the most stable (highest mean
# Jaccard) is retained; entropy and uncertainty are reported alongside.
assessment$Degenerate <- assessment$near_det == 1 | assessment$cond_max > 100
assessment$Admissible <- assessment$Jaccard >= 0.60 & !assessment$Degenerate

cand <- assessment[assessment$Admissible, ]
if (nrow(cand) == 0) stop("No admissible GMM solution: review the assessment.")
top  <- cand[which.max(cand$Jaccard), ]

sel_dim <- top$Dimension
sel_G   <- top$G
cat(sprintf("\n  Selected solution: FAMD = %d dimensions, G = %d\n", sel_dim, sel_G))
cat("  Admissible candidates:\n"); print(cand[order(-cand$Jaccard), c("Dimension", "G", "Jaccard", "MCse", "Entropy", "Uncertainty")])

# The manuscript is built on FAMD-4, G-4. A different selection means the text,
# Table 3 and every later stage must be revised before continuing.
if (!(sel_dim == 4 && sel_G == 4))
  stop("The selection rule no longer picks FAMD-4, G-4. Review before continuing.")

# ── TABLE C.3 ─────────────────────────────────────────────────────────────────
# Three- and four-dimensional spaces; the five-dimensional space is excluded in
# the text (no admissible solution, degenerate components).
tC3 <- assessment[assessment$Dimension %in% c(3, 4), ]
tC3_df <- data.frame(
  G           = tC3$G,
  Jaccard     = sprintf("%.3f", tC3$Jaccard),
  Entropy     = sprintf("%.3f", tC3$Entropy),
  Uncertainty = sprintf("%.3f", tC3$Uncertainty),
  stringsAsFactors = FALSE)

row_sel <- which(tC3$Dimension == sel_dim & tC3$G == sel_G)
rows_d3 <- which(tC3$Dimension == 3)
rows_d4 <- which(tC3$Dimension == 4)
lcss    <- "background-color: #e8e8e8; font-weight: bold;"

kblC3 <- kbl(
  tC3_df,
  format    = "html",
  caption   = "Table C.3. GMM-VEV multi-criteria assessment across candidate FAMD dimensions",
  col.names = c("G", "Jaccard (≥ 0.60)", "Entropy (↓)", "Mean uncertainty (↓)"),
  align     = c("c", "r", "r", "r"),
  row.names = FALSE) %>%
  kable_classic(full_width = FALSE, html_font = "Times New Roman") %>%
  kable_styling(font_size = 11) %>%
  row_spec(0, bold = TRUE) %>%
  row_spec(row_sel, bold = TRUE) %>%
  pack_rows("FAMD, 3 retained dimensions", min(rows_d3), max(rows_d3), label_row_css = lcss) %>%
  pack_rows("FAMD, 4 retained dimensions", min(rows_d4), max(rows_d4), label_row_css = lcss) %>%
  footnote(
    general = paste("Bold marks the retained solution. The five-dimension space",
                    "yielded no admissible solution and is not shown."),
    general_title     = "Notes:",
    footnote_as_chunk = TRUE)

print(kblC3)
save_table(kblC3, "Table_C3.html")

remove(kblC3, tC3, tC3_df, row_sel, rows_d3, rows_d4, lcss, cand, top)
gc()

# ── Five-dimensional space: degeneracy diagnostics ───────────────────────────
# Condition numbers of the component covariances, purity of each component on
# the binary active variables, and the spread of dimension-5 scores within the
# cells defined by those variables. Written to output/, not tabulated.
bin_vars <- active_vars[vapply(data_famd[active_vars], is.factor, logical(1))]
fits_d5  <- c(sprintf("d4_G%d", sel_G), sprintf("d5_G%d", G_values))
fits_d5  <- intersect(fits_d5, names(gmm_store))

eig_diag <- do.call(rbind, lapply(fits_d5, function(nm) {
  fit   <- gmm_store[[nm]]
  sigma <- fit$parameters$variance$sigma
  do.call(rbind, lapply(seq_len(dim(sigma)[3]), function(g) {
    ev <- eigen(sigma[, , g], symmetric = TRUE, only.values = TRUE)$values
    data.frame(Fit = nm, Component = g,
               pro = fit$parameters$pro[g],
               size = sum(fit$classification == g),
               eig_min = min(ev), eig_max = max(ev),
               cond_num = max(ev) / min(ev))
  }))
}))
purity <- do.call(rbind, lapply(fits_d5, function(nm) {
  cl <- gmm_store[[nm]]$classification
  do.call(rbind, lapply(sort(unique(cl)), function(g) {
    m <- cl == g
    data.frame(Fit = nm, Component = g,
               t(vapply(bin_vars, function(v) max(prop.table(table(data_famd[[v]][m]))),
                        numeric(1))), check.names = FALSE)
  }))
}))
coord5  <- famd_result$ind$coord[, 5]
cell_id <- interaction(data_famd[bin_vars], drop = TRUE)
band    <- data.frame(Cell = levels(cell_id),
                      n = as.integer(tapply(coord5, cell_id, length)),
                      sd_ratio = as.numeric(tapply(coord5, cell_id, sd)) / sd(coord5))
cat(sprintf("\n  Dimension-5 scores: within-cell SD / overall SD (n-weighted) = %.3f\n",
            weighted.mean(band$sd_ratio, band$n, na.rm = TRUE)))
save_output(eig_diag, "C4_dim5_condition_numbers.csv")
save_output(purity,   "C4_dim5_binary_purity.csv")
save_output(band,     "C4_dim5_band_spread.csv")
remove(eig_diag, purity, coord5, cell_id, band, fits_d5)

# ── Final GMM: REUSE the stored reference fit (no refit) ─────────────────────
# The final model is the exact object whose diagnostics appear in Table C.3,
# so Table 3 uncertainty matches the C.3 row by construction.
gmm_model     <- gmm_store[[sprintf("d%d_G%d", sel_dim, sel_G)]]
factor_scores <- famd_result$ind$coord[, 1:sel_dim, drop = FALSE]
selected_G    <- sel_G

optimal_k                      <- selected_G
post                           <- gmm_model$z
data_clean$profile_uncertainty <- 1 - apply(post, 1, max)

# Profile numbering. The mixture labels its components arbitrarily; profiles
# are numbered by weighted adoption, from the highest (P1) to the lowest. The
# outcome plays no role in the classification: it only orders the labels.
gmm_class      <- gmm_model$classification
adopt_by_class <- tapply(data_clean[[outcome_var]] * data_clean$weight_final, gmm_class, sum) /
  tapply(data_clean$weight_final, gmm_class, sum)

# Mixture components from the highest to the lowest weighted adoption
adopt_order <- as.integer(names(sort(adopt_by_class, decreasing = TRUE)))

# Profile number of each household = position of its component in that order
data_clean$profile_gmm <- match(gmm_class, adopt_order)

cat(sprintf("  Profile sizes: %s\n",
            paste(table(data_clean$profile_gmm), collapse = ", ")))
cat(sprintf("  Mean uncertainty: %.4f | log-likelihood: %.2f\n",
            mean(data_clean$profile_uncertainty), gmm_model$loglik))

# Cluster-wise stability of the retained solution, by profile
cw_sel <- clusterwise[clusterwise$Dimension == sel_dim & clusterwise$G == sel_G, ]
cw_sel$Profile <- match(cw_sel$Component, adopt_order)
cw_sel <- cw_sel[order(cw_sel$Profile), ]
cat("\n-- Cluster-wise stability of the retained solution --\n"); print(cw_sel, digits = 3)

# Check: the retained fit reproduces the published log-likelihood
# (-166239.64 on the 2023 data)
if (abs(gmm_model$loglik - (-166239.64)) > 0.05)
  warning("The retained fit does not reproduce the published log-likelihood.")

# profile_labels and the characterisation below assume the 4-profile solution.
remove(post, gmm_model, gmm_store, gmm_class, adopt_by_class, adopt_order, cw_sel)
gc()

# ══════════════════════════════════════════════════════════════════════════════
# 6. PROFILE CHARACTERIZATION
# ══════════════════════════════════════════════════════════════════════════════

# Step 1: compute profile characteristics
# All means are WEIGHTED. The FAMD is fitted with row.w = weight_final, so the
# partition itself is population-representative; describing it with unweighted
# means would be internally inconsistent. "n" is the raw sample count,
# not a weighted quantity, and "profile_uncertainty" is a model diagnostic
# for the estimation sample, so both stay unweighted on purpose.
gmm_char <- data_clean %>%
  group_by(profile_gmm) %>%
  summarise(
    n            = n(),
    weight_total = sum(weight_final, na.rm = TRUE),
    adopt        = weighted.mean(retrofit,         w = weight_final, na.rm = TRUE),
    inc          = weighted.mean(ln_income,        w = weight_final, na.rm = TRUE),
    pct_own      = weighted.mean(ownership,        w = weight_final, na.rm = TRUE),
    pct_edu      = weighted.mean(allhigheducation, w = weight_final, na.rm = TRUE),
    pct_det      = weighted.mean(detached,         w = weight_final, na.rm = TRUE),
    pct_moisture = weighted.mean(moisturedamage,   w = weight_final, na.rm = TRUE),
    rooms        = weighted.mean(rooms,            w = weight_final, na.rm = TRUE),
    pct_poll     = weighted.mean(pollution,        w = weight_final, na.rm = TRUE),
    HDD          = weighted.mean(HDD_1723,         w = weight_final, na.rm = TRUE),
    CDD          = weighted.mean(CDD_1723,         w = weight_final, na.rm = TRUE),
    age          = weighted.mean(meanageadults,    w = weight_final, na.rm = TRUE),
    pct_kids     = weighted.mean(kids,             w = weight_final, na.rm = TRUE),
    pct_highpop  = weighted.mean(highlyurbanised,  w = weight_final, na.rm = TRUE),
    uncert       = mean(profile_uncertainty, na.rm = TRUE),
    .groups = "drop"
  )

# Step 2: weighted population share of each profile. This column is used in the
# table header row and was missing from the previous version of the script.
gmm_char <- gmm_char %>%
  mutate(weighted_pct = 100 * weight_total / sum(weight_total))

# Check: weighted shares must add up to 100
cat(sprintf("  Weighted profile shares: %s (sum = %.1f%%)\n",
            paste(sprintf("%.1f%%", gmm_char$weighted_pct), collapse = ", "),
            sum(gmm_char$weighted_pct)))

# Check: weighted adoption rates by profile, needed for the manuscript text
cat(sprintf("  Weighted adoption rates: %s\n",
            paste(sprintf("P%s = %.2f%%", gmm_char$profile_gmm,
                          100 * gmm_char$adopt), collapse = ", ")))

# Check: profiles are numbered in decreasing order of weighted adoption
stopifnot(!is.unsorted(rev(gmm_char$adopt[order(gmm_char$profile_gmm)])))

# Check: weighted population adoption rate, for comparison
cat(sprintf("  Weighted population adoption rate: %.2f%%\n",
            100 * weighted.mean(data_clean[[outcome_var]],
                                w = data_clean$weight_final, na.rm = TRUE)))

# Profile labels, in profile order (P1 highest adoption, P4 lowest)
profile_labels <- c(
  "1" = "Educated young families in owner-occupied urban homes",
  "2" = "Older urban homeowners in compact, well-maintained apartments",
  "3" = "Cold-climate mixed-tenure households in damaged dwellings",
  "4" = "Warm-climate homeowners in detached homes"
)

# Short versions for the axis of Figure C.2 (stage 04)
profile_labels_short <- c(
  "1" = "Educated young\nfamilies, urban\nowner-occupiers (P1)",
  "2" = "Older urban\nhomeowners in\napartments (P2)",
  "3" = "Cold-climate\nmixed-tenure, damaged\ndwellings (P3)",
  "4" = "Warm-climate\nhomeowners in\ndetached homes (P4)"
)

# Step 3: format for table (one column per variable, in final row order).
# Rows follow the six blocks defined in Section 3 (tenure, income, education,
# climate, dwelling, household), so that Table 3, Table A.1 and the Shapley
# attribution all group the variables in the same way. Adoption rate and
# sample size head the table; GMM uncertainty closes it as a diagnostic.
kbl3_df <- gmm_char %>%
  arrange(profile_gmm) %>%
  transmute(
    Profile = paste0("P", profile_gmm, ": ", profile_labels[as.character(profile_gmm)]),
    
    # Header block: size and outcome
    `Unweighted survey n (weighted %)` = sprintf("%s (%.1f%%)", format(n, big.mark = ","), weighted_pct),
    `Adoption rate`                   = scales::percent(adopt,        accuracy = 0.1),
    
    # A. Tenure
    `% Ownership`                     = scales::percent(pct_own,      accuracy = 0.1),
    
    # B. Income
    `Income (log)`                    = round(inc, 2),
    
    # C. Education
    `% Highly educated`               = scales::percent(pct_edu,      accuracy = 0.1),
    
    # D. Climatic exposure
    `HDD (mean '17-23)`               = round(HDD, 2),
    `CDD (mean '17-23)`               = round(CDD, 2),
    
    # E. Dwelling characteristics
    `Number of rooms`                 = round(rooms, 2),
    `% Detached dwelling`             = scales::percent(pct_det,      accuracy = 0.1),
    `% Moisture / structural damage`  = scales::percent(pct_moisture, accuracy = 0.1),
    `% Polluted environment`          = scales::percent(pct_poll,     accuracy = 0.1),
    
    # F. Household composition and location
    `Mean age of adults`              = round(age, 2),
    `% Children present`              = scales::percent(pct_kids,     accuracy = 0.1),
    `% Highly urbanised municipality` = scales::percent(pct_highpop,  accuracy = 0.1),
    
    # Diagnostic
    `Uncertainty`                     = round(uncert, 2)
  )

# Step 4: transpose (variables as rows, profiles as columns) and render
kbl3_df_t           <- as.data.frame(t(kbl3_df), stringsAsFactors = FALSE)
colnames(kbl3_df_t) <- as.character(unlist(kbl3_df_t[1, ]))  # profile labels -> headers
kbl3_df_t           <- kbl3_df_t[-1, , drop = FALSE] # drop the Profile row

# Row positions of each block after transposing. Rows 1-2 are the header block
# (n and adoption rate) and are left without a section label.
r_tenure    <- 3
r_income    <- 4
r_education <- 5
r_climate   <- 6:7
r_dwelling  <- 8:11
r_household <- 12:14
r_quality   <- 15

kbl3 <- kbl(kbl3_df_t,
            format  = "html",
            caption = "Table 3. GMM profile characterisation",
            align   = "l") %>%
  kable_classic(full_width = FALSE, html_font = "Times New Roman") %>%
  kable_styling(font_size = 10) %>%
  row_spec(0, bold = TRUE) %>%
  column_spec(1, width = "150px") %>%
  pack_rows("A. Tenure",                             r_tenure,          r_tenure) %>%
  pack_rows("B. Income",                             r_income,          r_income) %>%
  pack_rows("C. Education",                          r_education,       r_education) %>%
  pack_rows("D. Climatic exposure",                  min(r_climate),    max(r_climate)) %>%
  pack_rows("E. Dwelling characteristics",           min(r_dwelling),   max(r_dwelling)) %>%
  pack_rows("F. Household composition and location", min(r_household),  max(r_household)) %>%
  pack_rows("GMM quality",                           r_quality,         r_quality) %>%
  footnote(
    general           = paste("All figures are weighted using survey weights,",
                              "except the sample count n and the GMM uncertainty."),
    general_title     = "Notes:",
    footnote_as_chunk = TRUE
  )

print(kbl3)
save_table(kbl3, "Table_3.html")

# weighted_pct is a column of gmm_char, not a standalone object
remove(gmm_char, kbl3_df, kbl3_df_t, kbl3,
       r_tenure, r_income, r_education, r_climate, r_dwelling,
       r_household, r_quality)
gc()

# ══════════════════════════════════════════════════════════════════════════════
# SAVE INTERMEDIATE OBJECTS FOR THE NEXT STAGE
# ══════════════════════════════════════════════════════════════════════════════

save(data_clean, logit_baseline, formula_baseline, ame_df,
     active_vars, var_labels, outcome_var, profile_labels, profile_labels_short,
     famd_result,
     file = "03_famd_clustering.RData")

cat(sprintf("\nStage 03 complete. Runtime: %.1f minutes\n",
            as.numeric(difftime(Sys.time(), start_time, units = "mins"))))
