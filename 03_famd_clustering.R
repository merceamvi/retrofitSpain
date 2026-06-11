###############################################################################
# Project       : [rehabSpain] [3] FAMD & GMM CLUSTERING
# Creation date : 02/12/2024
# Last update   : 09/06/2026
# Author        : Mercè Amich (merce.amich@ehu.eus)
# Institution   : UPV/EHU, BC3
# Last run time : —  (heavy: Ward + GMM bootstrap stability)
# Script Overview:
#   Factor Analysis of Mixed Data (FAMD), multi-criteria cluster assessment
#   (Ward + GMM-VEV bootstrap stability over candidate dimensions), final
#   GMM fit and household-profile characterisation.
# Requirements:
#   - This script must be in the same directory as:
#       a) _setup.R
#       b) 02_descriptives_logit.RData   (output of 02_descriptives_logit.R)
# Output  : Publication-ready tables (kableExtra) + scientific ggplots (TIFF)
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

# Load objects from the previous stage
#   (data_clean, logit_baseline, formula_baseline, ame_df,
#    active_vars, var_labels, outcome_var)
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
                                ncp     = 10, # Should comment
                                sup.var = idx_sup,
                                row.w   = data_clean$weight_final,
                                graph   = FALSE)

# Extract eigenvalues
eig_df           <- as.data.frame(famd_result$eig)
colnames(eig_df) <- c("Eigenvalue", "Variance_%", "Cumulative_%")
eig_df$Dimension <- seq_len(nrow(eig_df))

# Check Kaiser criterion (eigenvalue > 1) results (informative)
n_factors <- sum(eig_df$Eigenvalue > 1)

cat(sprintf("  Kaiser criterion (informative): keep %d factors (%.1f%% cumulative variance)\n",
            n_factors, eig_df$`Cumulative_%`[n_factors]))

# TABLE S4: Eigenvalues
kblS4 <- eig_df %>%
  filter(Dimension <= 7) %>%
  mutate(Retained = ifelse(Dimension %in% 3:5, "Candidate", 
                           ifelse(Dimension <= 5, "Yes", "—")),
         across(c(Eigenvalue, `Variance_%`, `Cumulative_%`), ~ round(., 2))) %>%
  select(Dimension, Eigenvalue, `Variance_%`, `Cumulative_%`, Retained) %>%
  kbl(format = "html",
      caption = "Table S4. FAMD eigenvalues",
      col.names = c("Dim.", "Eigenvalue", "Var. (%)", "Cumul. (%)", "Retained"),
      align = c("c", "r", "r", "r", "c")) %>%
  kable_classic(full_width = FALSE, html_font = "Times New Roman") %>%
  kable_styling(font_size = 11) %>%
  row_spec(0, bold = TRUE) %>%
  row_spec(1:5, bold = TRUE, background = "#EBF3FB") %>%
  footnote(general = paste0(
    "Dimensions 1-5 retained based on scree plot elbow (Cattell, 1966). ",
    "Dimensions 3-5 mark candidate retention points."
  ),
  general_title = "Notes:", footnote_as_chunk = TRUE)

print(kblS4)
save_table(kblS4, "TableS4_FAMD_eigenvalues.html")

# FIGURE 2: Scree plot
fig2 <- ggplot(eig_df %>% filter(Dimension <= 8),
               aes(x = Dimension, y = `Variance_%`)) +
  geom_col(aes(fill = Eigenvalue > 1),
           width = 0.65, color = "gray20", linewidth = 0.3) +
  geom_line(linewidth = 0.6, color = "gray40") +
  geom_point(size = 2.5, color = "black", fill = "white", shape = 21, stroke = 0.5) +
  geom_hline(yintercept = 100 / length(active_vars),
             linetype = "dashed", color = "gray40", linewidth = 0.4) +
  scale_fill_manual(values = c("TRUE" = "gray50", "FALSE" = "gray85"),
                    labels = c("TRUE" = "Retained", "FALSE" = "Not retained"),
                    name = NULL) +
  scale_x_continuous(breaks = 1:8) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title    = NULL,
       x        = "Dimension", 
       y        = "Variance explained (%)",
       caption  = "Notes: Dashed line indicates average variance under uniform distribution (8.3%).\nDimensions 1-5 retained based on scree plot elbow.") +
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

print(fig2)
ggsave("Figure2_FAMD_scree.pdf", fig2, width = 6.5, height = 3.5, device = cairo_pdf)
ggsave("Figure2_FAMD_screeplot.tiff",
       plot   = fig2,
       device = "tiff",
       width  = 14,
       height = 9,
       units  = "cm",
       dpi    = 900)

# Extract FAMD components
var_contrib   <- famd_result$var$contrib    [, 1:n_factors, drop = FALSE]
var_coord     <- famd_result$var$coord      [,  1:n_factors, drop = FALSE]
sup_coord     <- famd_result$quali.sup$coord[, 1:n_factors, drop = FALSE]
factor_scores <- famd_result$ind$coord      [, 1:n_factors, drop = FALSE]

# Extract retrofit (== 1) coordinates
ret_row <- grepl("^Yes$", rownames(sup_coord))
ret_yes <- as.numeric(sup_coord[ret_row, ])

# TABLE 2: Dimension interpretation
dim2 <- data.frame(
  Dim = 1:n_factors,
  `Var (%)` = round(eig_df$`Variance_%`[1:n_factors], 2),
  `Top 3` = sapply(1:n_factors, function(d) {
    cc <- var_contrib[, d, drop = TRUE]
    paste(names(sort(cc[!is.na(cc)], decreasing = TRUE)[1:3]), collapse = ", ")
  }),
  `|lambda|>0.4` = sapply(1:n_factors, function(d) {
    cl <- var_coord[, d, drop = TRUE]
    cl <- cl[!is.na(cl)]
    s  <- names(cl[abs(cl) > 0.4])
    if (length(s) == 0) "—" else paste(s, collapse = ", ")
  }),
  `Retrofit` = round(ret_yes[1:n_factors], 2),
  check.names = FALSE
)
kbl2 <- kbl(dim2, format = "html",
            caption = "Table 2. FAMD dimension interpretation",
            align = c("c", "r", "l", "l", "r")) %>%
  kable_classic(full_width = FALSE, html_font = "Times New Roman") %>%
  kable_styling(font_size = 11) %>%
  row_spec(0, bold = TRUE) %>%
  footnote(general = "Retrofit = G_s(Adopt = Yes). Variables with |G_s(·)| > 0.4 shown.",
           general_title = "Notes:", footnote_as_chunk = TRUE)

print(kbl2)
save_table(kbl2, "Table2_FAMD_dimensions.html")

# Clean up
remove(eig_df, var_contrib, var_coord, sup_coord, ret_row, ret_yes, 
       dim2, kbl2, fig2, idx_sup)
gc()


# ══════════════════════════════════════════════════════════════════════════════
# MULTI-CRITERIA ASSESSMENT (WARD + GMM-VEV on FAMD=3:5)
# ══════════════════════════════════════════════════════════════════════════════

dims     <- c(3, 4, 5) # Loop over 3, 4 and 5 retained dimensions
n_bs     <- 200        # 200 bootstrap iterations
G_values <- 3:6        # For k=3:6 in GMM-VEV
n_boot   <- 50         # Number of bootstrap iterations

# Initialise vectors for storing results
optimal_k_vector      <- numeric(length(dims))
silhouette_vector     <- numeric(length(dims))
ward_stability_vector <- numeric(length(dims))
all_quality_list      <- vector("list", length(dims))

# Loop
for (i in seq_along(dims)) {
  
  set.seed(123) # Set seed for reproducibility
  
  dim <- dims[i] # For each dim,
  cat(sprintf("=== Processing dimension: %d ===\n", dim)) # Print control message,
  
  # Compute factor scores with n. dimensions retained = dims[i]
  factor_scores <- famd_result$ind$coord[, 1:dim, drop = FALSE]
  
  # ── Ward ────────────────────────────────────────────────────────────────────
  dist_mat    <- dist(factor_scores, method = "euclidean")
  hclust_ward <- hclust(dist_mat,    method = "ward.D2")
  
  # Silhouette (k=3:5)
  sil_df <- data.frame(k = 3:5, avg_sil = NA_real_)
  for (k in 3:5) {
    cl  <- cutree(hclust_ward, k = k)
    sil <- silhouette(cl, dist_mat)
    sil_df$avg_sil[sil_df$k == k] <- mean(sil[, "sil_width"])
  }
  
  best_sil  <- max(sil_df$avg_sil)                 # Best silhouette
  optimal_k <- sil_df$k[which.max(sil_df$avg_sil)] # Optimal k
  optimal_k_vector[i]  <- optimal_k # Extract and store
  silhouette_vector[i] <- best_sil  # Extract and store
  cat(sprintf("  Ward optimal k: %d (silhouette = %.3f)\n", optimal_k, best_sil))
  
  data_clean[[paste0("profile_ward_", dim)]] <- cutree(hclust_ward, k = optimal_k)
  
  # ── Ward bootstrap stability ─────────────────────────────────────────────────
  cat("  Running Ward bootstrap stability...\n")
  bstab_clusters <- matrix(NA, nrow(data_clean), n_bs)
  
  for (b in seq_len(n_bs)) {
    
    # Extract indices
    idx_b <- sample(nrow(data_clean), replace = TRUE)
    
    # Run FAMD on ncp = dim
    fb  <- FactoMineR::FAMD(data_famd[idx_b, ],
                            ncp     = dim,
                            sup.var = match(outcome_var, names(data_famd)),
                            row.w   = data_clean$weight_final[idx_b],
                            graph   = FALSE)
    
    hcb <- hclust(dist(fb$ind$coord[, 1:dim]), method = "ward.D2")
    bstab_clusters[idx_b, b] <- cutree(hcb, k = optimal_k)
    
    if (b %% 20 == 0) cat(sprintf("    Bootstrap %d/%d\n", b, n_bs))
    rm(fb, hcb); gc()
  }
  
  ward_stab <- apply(bstab_clusters, 1, function(x) {
    mc <- as.numeric(names(sort(table(x), decreasing = TRUE)[1]))
    mean(x == mc, na.rm = TRUE)
  })
  
  ward_stability_vector[i] <- mean(ward_stab, na.rm = TRUE)
  cat(sprintf("  Ward cluster stability: %.2f%%\n", 100 * ward_stability_vector[i]))
  rm(bstab_clusters, ward_stab, hclust_ward); gc()
  
  # ── LOOP 1: GMM Jaccard over G = 3:6 ────────────────────────────────────────
  cat("  GMM bootstrap stability (G = 3:6)...\n")
  stab_results <- data.frame(G = 3:6, mean_jaccard = NA_real_)
  
  for (g in 3:6) {
    cat(sprintf("    G = %d...\n", g))
    
    gmm_g <- tryCatch(
      Mclust(factor_scores, G = g, modelNames = "VEV", verbose = FALSE),
      error = function(e) NULL
    )
    if (is.null(gmm_g)) next
    
    jacc_vec <- numeric(n_boot)
    for (b in seq_len(n_boot)) {
      idx_b <- sample(nrow(factor_scores), replace = TRUE)
      
      gmm_b <- tryCatch(
        Mclust(factor_scores[idx_b, ], G = g, modelNames = "VEV", verbose = FALSE),
        error = function(e) NULL
      )
      if (is.null(gmm_b)) { jacc_vec[b] <- NA; next }
      
      cl_full <- gmm_g$classification[idx_b]
      cl_boot <- gmm_b$classification
      cont    <- table(cl_full, cl_boot)
      jacc_vec[b] <- mean(apply(cont, 1, function(r) {
        j <- max(r) / (sum(r) + max(cont[, which.max(r)]) - max(r))
        ifelse(is.nan(j), NA, j)
      }), na.rm = TRUE)
      
      rm(gmm_b, cl_full, cl_boot, cont); gc()
    }
    
    stab_results$mean_jaccard[stab_results$G == g] <- mean(jacc_vec, na.rm = TRUE)
    cat(sprintf("      Jaccard = %.3f\n", mean(jacc_vec, na.rm = TRUE)))
    rm(gmm_g, jacc_vec); gc()
  }
  
  # ── LOOP 2: Multi-criteria quality assessment G = 3:6 ────────────────────────
  cat("  Multi-criteria quality assessment (G = 3:6)...\n")
  cluster_quality <- data.frame()
  dist_mat        <- dist(factor_scores, method = "euclidean")
  
  for (g in 3:6) {
    gmm_g <- tryCatch(
      Mclust(factor_scores, G = g, modelNames = "VEV", verbose = FALSE),
      error = function(e) NULL
    )
    if (is.null(gmm_g)) next
    
    jaccard        <- stab_results$mean_jaccard[stab_results$G == g]
    sil            <- silhouette(gmm_g$classification, dist_mat)
    avg_silhouette <- mean(sil[, "sil_width"])
    
    z <- gmm_g$z
    z[z < 1e-10] <- 1e-10
    norm_entropy     <- mean(-rowSums(z * log(z))) / log(g)
    mean_uncertainty <- mean(1 - apply(z, 1, max))
    
    cluster_quality <- rbind(cluster_quality, data.frame(
      Dimension   = dim,
      G           = g,
      Jaccard     = round(jaccard,          3),
      Entropy     = round(norm_entropy,     3),
      Silhouette  = round(avg_silhouette,   3),
      Uncertainty = round(mean_uncertainty, 3)
    ))
    rm(gmm_g, sil, z); gc()
  }
  
  print(cluster_quality)
  cat(sprintf("=== Dimension %d done ===\n\n", dim))
  
  all_quality_list[[i]] <- cluster_quality   
  
}

names(optimal_k_vector)      <- dims
names(silhouette_vector)     <- dims
names(ward_stability_vector) <- dims

cat("Optimal k per dimension:       ", optimal_k_vector,                       "\n")
cat("Best silhouette per dimension: ", round(silhouette_vector,        3),     "\n")
cat("Ward bootstrap stability (%):  ", round(100 * ward_stability_vector, 2), "\n")

# CHECK:
all_quality_list

# TABLE S5: Full cluster quality assessment
all_quality  <- do.call(rbind, all_quality_list)
ward_jaccard <- setNames(round(ward_stability_vector * 100, 2), as.character(dims))

selected_k <- c("3" = 3, "4" = 4, "5" = 4)

highlight_rows <- which(
  mapply(function(d, k) all_quality$Dimension == d & all_quality$G == k,
         as.integer(names(selected_k)),
         selected_k) |>
    apply(1, any)
)

rows_d3 <- which(all_quality$Dimension == 3)
rows_d4 <- which(all_quality$Dimension == 4)
rows_d5 <- which(all_quality$Dimension == 5)
lcss    <- "background-color: #e8e8e8; font-weight: bold;"

kblS5 <- kbl(
  all_quality[, -1],
  format    = "html",
  caption   = "Table S5. GMM-VEV cluster quality diagnostics and Ward stability assessment",
  col.names = c("k", "Jaccard (> 0.60)", "Entropy (↓)",
                "Silhouette (↑)", "Mean uncertainty (↓)"),
  align     = c("c", "r", "r", "r", "r"),
  row.names = FALSE,
  digits    = 3
) %>%
  kable_classic(full_width = FALSE, html_font = "Times New Roman") %>%
  kable_styling(font_size = 11) %>%
  row_spec(0, bold = TRUE) %>%
  row_spec(highlight_rows, background = "#d1ecf1", bold = TRUE) %>%
  pack_rows(
    sprintf("FAMD = 3 retained dimensions  |  Ward optimal k = %d  |  Ward Jaccard = %.2f%%",
            optimal_k_vector[1], ward_jaccard["3"]),
    min(rows_d3), max(rows_d3), label_row_css = lcss
  ) %>%
  pack_rows(
    sprintf("FAMD = 4 retained dimensions  |  Ward optimal k = %d  |  Ward Jaccard = %.2f%%",
            optimal_k_vector[2], ward_jaccard["4"]),
    min(rows_d4), max(rows_d4), label_row_css = lcss
  ) %>%
  pack_rows(
    sprintf("FAMD = 5 retained dimensions  |  Ward optimal k = %d  |  Ward Jaccard = %.2f%%",
            optimal_k_vector[3], ward_jaccard["5"]),
    min(rows_d5), max(rows_d5), label_row_css = lcss
  ) %>%
  footnote(
      general = paste0(
        "Highlighted rows: best k per FAMD specification."
      ),
      general_title     = "Notes:",
      footnote_as_chunk = TRUE
    )


print(kblS5)
save_table(kblS5, "TableS5_Cluster_quality_full.html")

remove(kblS5, all_quality, all_quality_list, highlight_rows,
       rows_d3, rows_d4, rows_d5)
gc()



# ── Fit final GMM at selected_G ────
set.seed(123)
factor_scores <- famd_result$ind$coord[, 1:4, drop = FALSE]
selected_G <- 4

gmm_model <- Mclust(factor_scores, 
                    G          = selected_G, 
                    modelNames = "VEV", 
                    verbose    = FALSE)

data_clean$profile_gmm     <- gmm_model$classification
optimal_k                  <- selected_G
post                       <- gmm_model$z
data_clean$profile_uncertainty <- 1 - apply(post, 1, max)

cat(sprintf("  Profile sizes: %s\n",
            paste(table(data_clean$profile_gmm), collapse = ", ")))
cat(sprintf("  Mean uncertainty: %.4f\n", mean(data_clean$profile_uncertainty)))

remove(post, gmm_model, cluster_quality)
gc()

# ══════════════════════════════════════════════════════════════════════════════
# 6. PROFILE CHARACTERIZATION
# ══════════════════════════════════════════════════════════════════════════════

# Extract profile characteristics
gmm_char <- data_clean %>%
  group_by(profile_gmm) %>%
  summarise(
    n            = n(),
    adopt        = mean(retrofit, na.rm = TRUE),
    inc_z        = mean(ln_income_z, na.rm = TRUE),
    pct_own      = mean(ownership, na.rm = TRUE),
    pct_edu      = mean(allhigheducation, na.rm = TRUE),
    pct_det      = mean(detached, na.rm = TRUE),
    pct_moisture = mean(moisturedamage, na.rm = TRUE),
    rooms_z      = mean(rooms_z, na.rm = TRUE),
    pct_poll     = mean(pollution, na.rm = TRUE),
    HDD_z        = mean(HDD_1723_z, na.rm = TRUE),
    CDD_z        = mean(CDD_1723_z, na.rm = TRUE),
    age_z        = mean(meanageadults_z, na.rm = TRUE),
    pct_kids     = mean(kids, na.rm = TRUE),
    pct_highpop  = mean(highpopulation, na.rm = TRUE),
    uncert       = mean(profile_uncertainty, na.rm = TRUE),
    .groups = "drop"
  )

# Calculate weighted percentages
weighted_pct <- data_clean %>%
  group_by(profile_gmm) %>%
  summarise(weighted_n = sum(weight_final), .groups = "drop") %>%
  mutate(weighted_pct = 100 * weighted_n / sum(weighted_n))

gmm_char <- gmm_char %>%
  left_join(weighted_pct, by = "profile_gmm")

# Define profile labels
profile_labels <- c(
  "1" = "Cold-climate semi-rural homeowners in large detached homes",
  "2" = "Educated urban homeowners in dense apartment blocks",
  "3" = "Warm-climate homeowners in detached homes",
  "4" = "Precarious young renters in small damaged urban dwellings"
)

# Step 1: compute profile characteristics
gmm_char <- data_clean %>%
  group_by(profile_gmm) %>%
  summarise(
    n            = n(),
    adopt        = mean(retrofit,            na.rm = TRUE),
    inc          = mean(ln_income,           na.rm = TRUE),
    pct_own      = mean(ownership,           na.rm = TRUE),
    pct_edu      = mean(allhigheducation,    na.rm = TRUE),
    pct_det      = mean(detached,            na.rm = TRUE),
    pct_moisture = mean(moisturedamage,      na.rm = TRUE),
    rooms        = mean(rooms,               na.rm = TRUE),
    pct_poll     = mean(pollution,           na.rm = TRUE),
    HDD          = mean(HDD_1723,            na.rm = TRUE),
    CDD          = mean(CDD_1723,            na.rm = TRUE),
    age          = mean(meanageadults,       na.rm = TRUE),
    pct_kids     = mean(kids,                na.rm = TRUE),
    pct_highpop  = mean(highpopulation,      na.rm = TRUE),
    uncert       = mean(profile_uncertainty, na.rm = TRUE),
    .groups = "drop"
  )

# Step 2: add weighted percentages
weighted_pct <- data_clean %>%
  group_by(profile_gmm) %>%
  summarise(weighted_n = sum(weight_final), .groups = "drop") %>%
  mutate(weighted_pct = 100 * weighted_n / sum(weighted_n))
gmm_char <- gmm_char %>%
  left_join(weighted_pct, by = "profile_gmm")

# Step 3: format for table
kbl3_df <- gmm_char %>%
  mutate(
    Profile          = paste0("P", profile_gmm, ": ", profile_labels[as.character(profile_gmm)]),
    `n (%)*`         = sprintf("%s (%.1f%%)", format(n, big.mark = ","), weighted_pct),
    `Adoption rate`  = scales::percent(adopt,        accuracy = 0.1),
    `% owners`       = scales::percent(pct_own,      accuracy = 0.1),
    `% high-educ`    = scales::percent(pct_edu,      accuracy = 0.1),
    `% detached`     = scales::percent(pct_det,      accuracy = 0.1),
    `% moisture`     = scales::percent(pct_moisture, accuracy = 0.1),
    `% polluted`     = scales::percent(pct_poll,     accuracy = 0.1),
    `% kids`         = scales::percent(pct_kids,     accuracy = 0.1),
    `% high-pop`     = scales::percent(pct_highpop,  accuracy = 0.1),
    across(c(inc, rooms, HDD, CDD, age, uncert), ~ round(., 2))
  ) %>%
  rename(
    `Income (log)`      = inc,
    `Rooms`             = rooms,
    `HDD (mean '17-23)` = HDD,
    `CDD (mean '17-23)` = CDD,
    `Mean age`          = age,
    `Uncertainty`       = uncert
  ) %>%
  select(
    Profile, `n (%)*`, `Adoption rate`, `Income (log)`, `% owners`,
    `% high-educ`, `% detached`, `% moisture`, `Rooms`, `% polluted`,
    `HDD (mean '17-23)`, `CDD (mean '17-23)`, `Mean age`,
    `% kids`, `% high-pop`, `Uncertainty`
  )

# Step 4: transpose and render
kbl3_df_t <- as.data.frame(t(kbl3_df))
kbl3 <- kbl(kbl3_df_t,
            format  = "html",
            caption = "Table 3. GMM profile characterisation",
            align   = "l") %>%
  kable_classic(full_width = FALSE, html_font = "Times New Roman") %>%
  kable_styling(font_size = 10) %>%
  row_spec(0, bold = TRUE) %>%
  column_spec(1, width = "150px") %>%
  footnote(
    general = "* Weighted percentages. Continuous variables: unweighted profile means. Uncertainty = 1 - max(posterior probability).",
    general_title     = "Notes:",
    footnote_as_chunk = TRUE
  )

print(kbl3)
save_table(kbl3, "Table3_GMM_profiles_transposed.html")

remove(gmm_char, weighted_pct, kbl3_df, kbl3_df_t, kbl3)
gc()


# ══════════════════════════════════════════════════════════════════════════════
# SAVE INTERMEDIATE OBJECTS FOR THE NEXT STAGE
# ══════════════════════════════════════════════════════════════════════════════

save(data_clean, logit_baseline, formula_baseline, ame_df,
     active_vars, var_labels, outcome_var, profile_labels, famd_result,
     file = "03_famd_clustering.RData")

cat(sprintf("\nStage 03 complete. Runtime: %.1f minutes\n",
            as.numeric(difftime(Sys.time(), start_time, units = "mins"))))
