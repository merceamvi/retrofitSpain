###############################################################################
# Project       : [rehabSpain] [4] PROFILE-SPECIFIC ANALYSIS
# Creation date : 02/12/2024
# Last update   : 23/11/2026 (profiles numbered by adoption in stage 03; Figure C.2
#                 axis labels taken from stage 03; Table D.2 note as in the paper;
#                 English object names; LR console check removed; bootstrap seed
#                 set per profile; AMEs of continuous variables scaled by the
#                 weighted standard deviation; glm warnings silenced)
# Author        : Merce Amich (merce.amich@ehu.eus)
# Institution   : UPV/EHU, BC3
# Last run time : ~5 min.

# Script Overview:
#   Profile-specific adoption analysis on the robust GMM partition:
#     - profile-specific HC1-robust logit coefficients (Table D.1)
#     - profile-specific bootstrapped AMEs (Table D.2)
#     - observed vs predicted adoption by profile, a calibration check of the
#       baseline logit (Figure C.2)

# Requirements: _setup.R and 03_famd_clustering.RData in the same directory
# Output  : paper/ Table_D1, Table_D2 (HTML) and Figure_C2 (PDF, TIFF)
###############################################################################

# ══════════════════════════════════════════════════════════════════════════════
# 0. SETUP
# ══════════════════════════════════════════════════════════════════════════════

rm(list = ls(all = TRUE))
Sys.setenv(LANG = "en")
start_time <- Sys.time()

source("_setup.R")
load("03_famd_clustering.RData")

# ══════════════════════════════════════════════════════════════════════════════
# PROFILE-SPECIFIC ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════

# Get predicted probabilities from baseline model
data_clean$pred_prob <- predict(logit_baseline, type = "response")

# Calculate observed vs predicted adoption by profile (survey-weighted)
prof_adopt <- data_clean %>%
  group_by(profile_gmm) %>%
  summarise(n    = n(),
            obs  = weighted.mean(retrofit,  w = weight_final, na.rm = TRUE),
            pred = weighted.mean(pred_prob, w = weight_final, na.rm = TRUE),
            .groups = "drop")
gc()

# ── Profile-specific HC1-robust coefficients (Table D.1) ─────────────────────
# One weighted logit per profile subsample with heteroskedasticity-robust (HC1)
# covariance, mirroring the pooled Table B.1 estimator

profiles          <- sort(unique(data_clean$profile_gmm))
profile_coef_list <- lapply(profiles, function(p) {
  dp <- data_clean %>% filter(profile_gmm == p)
  m  <- suppressWarnings(glm(formula_baseline, family = binomial("logit"),
                             data = dp, weights = weight_final))
  ct <- coeftest(m, vcov = vcovHC(m, type = "HC1"))
  data.frame(
    profile  = p,
    n        = nobs(m),
    variable = rownames(ct),
    estimate = ct[, 1],
    se       = ct[, 2],
    p_value  = ct[, 4],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
})

profile_coef <- bind_rows(profile_coef_list) %>%
  mutate(
    stars = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      p_value < 0.10  ~ ".",
      TRUE            ~ ""
    ),
    label = ifelse(variable == "(Intercept)", "Intercept", var_labels[variable])
  )

prof_n <- profile_coef %>% distinct(profile, n) %>% arrange(profile)

# ── Profile-specific AMEs (Table D.2) ────
cat("  Profile-specific AMEs (1,000 bootstrap iterations per profile)...\n")
n_boot_p     <- 1000
prf_ame_list <- lapply(profiles, function(p) {
  cat(sprintf("    Profile %d/%d\n", p, max(profiles)))
  
  # Same seed for every profile, so each profile's intervals do not depend on
  # the order in which the profiles are processed
  set.seed(123)
  
  dp <- data_clean %>% filter(profile_gmm == p)
  
  valid <- active_vars[sapply(active_vars, function(v)
    length(unique(dp[[v]][!is.na(dp[[v]])])) >= 2)]
  
  if (!length(valid)) return(NULL)
  
  bfn <- function(data, idx) {
    d <- as.data.frame(data[idx, ])
    m <- tryCatch(
      suppressWarnings(glm(formula_baseline, family = binomial("logit"),
                           data = d, weights = weight_final)),
      error = function(e) NULL
    )
    if (is.null(m) || !inherits(m, "glm")) return(rep(NA_real_, length(valid)))
    
    # Survey weights for population-consistent averaging of the AME operator.
    # The glm is weighted, so the AME average must be weighted too
    w <- d$weight_final
    
    a <- setNames(numeric(length(valid)), valid)
    for (v in valid) {
      uv  <- unique(d[[v]][!is.na(d[[v]])])
      bin <- length(uv) == 2 && all(sort(uv) %in% c(0, 1))
      
      if (bin) {
        # Binary variable: discrete first difference
        d0 <- d; d0[[v]] <- 0L
        d1 <- d; d1[[v]] <- 1L
        a[v] <- weighted.mean(predict(m, newdata = d1, type = "response") -
                                predict(m, newdata = d0, type = "response"),
                              w = w, na.rm = TRUE)
      } else {
        # Continuous variable: marginal effect scaled by the weighted SD
        pp     <- predict(m, newdata = d, type = "response")
        sd_x   <- sqrt(Hmisc::wtd.var(d[[v]], weights = w, normwt = TRUE, na.rm = TRUE))
        a[v]   <- coef(m)[v] * weighted.mean(pp * (1 - pp), w = w, na.rm = TRUE) * sd_x
      }
    }
    a * 100
  }
  
  br <- boot(dp, bfn, R = n_boot_p)
  
  # Point estimate on the full profile subsample, same reasoning as script 02
  ame_point_p <- bfn(dp, seq_len(nrow(dp)))
  
  data.frame(profile = p, variable = valid,
             ame   = ame_point_p,
             se    = apply(br$t, 2, sd, na.rm = TRUE),
             ci_lo = apply(br$t, 2, quantile, probs = 0.025, na.rm = TRUE),
             ci_hi = apply(br$t, 2, quantile, probs = 0.975, na.rm = TRUE))
})

prf_ame <- bind_rows(prf_ame_list) %>%
  mutate(sig   = !(ci_lo < 0 & ci_hi > 0),
         label = var_labels[variable])

# ── Row structure for Tables D.1 and D.2 ─────────────────────────────────────
# The row order and the block headers are derived from var_order and block_map,
# both defined in _setup.R, so that Tables D.1 and D.2 group the variables in
# the same six blocks used in Section 3 and in Tables 3, 4, A.1, B.1 and B.2.
# Editing the blocks there propagates here without touching this script.

# Variable labels in canonical order
canonical_labels <- unname(var_labels[var_order])

# Block of each variable, as an ordered factor
var_block <- factor(block_map[var_order], levels = block_order)

# Row structure for Table D.2 (AMEs): six blocks, no intercept
tblD2_structure <- split(canonical_labels, var_block)
tblD2_structure <- tblD2_structure[lengths(tblD2_structure) > 0]

ordered_labels_d2 <- unlist(tblD2_structure, use.names = FALSE)

d2_sizes  <- lengths(tblD2_structure)
d2_ends   <- cumsum(d2_sizes)
d2_starts <- d2_ends - d2_sizes + 1L

# Row structure for Table D.1 (HC1 coefficients): the intercept first, then the
# same six blocks
tblD1_structure <- c(list("Intercept" = "Intercept"), tblD2_structure)

ordered_labels_d1 <- unlist(tblD1_structure, use.names = FALSE)

d1_sizes  <- lengths(tblD1_structure)
d1_ends   <- cumsum(d1_sizes)
d1_starts <- d1_ends - d1_sizes + 1L

# Check: every active variable must appear exactly once in each structure
stopifnot(setequal(ordered_labels_d2, unname(var_labels[active_vars])))

# TABLE D.1: Profile-specific HC1-robust logit coefficients
# Cell: estimate (with significance stars) over robust SE in parentheses.
d1_lcss <- "background-color: #f0f0f0; font-weight: bold;"

tblD1_df <- profile_coef %>%
  mutate(cell = sprintf("%.3f%s<br>(%.3f)", estimate, stars, se)) %>%
  select(profile, label, cell) %>%
  pivot_wider(names_from = profile, values_from = cell, names_prefix = "Profile ") %>%
  rename(Variable = label) %>%
  mutate(Variable = factor(Variable, levels = ordered_labels_d1)) %>%
  arrange(Variable) %>%
  mutate(Variable = as.character(Variable)) %>%
  mutate(across(-Variable, ~ ifelse(is.na(.) | grepl("NA", .), "-", .)))

tblD1_html <- kbl(tblD1_df,
                  format  = "html",
                  escape  = FALSE,
                  caption = "Table D.1. Profile-specific binary logit coefficient estimates with heteroskedasticity-robust standard errors (HC1)",
                  align   = c("l", rep("c", ncol(tblD1_df) - 1))) %>%
  kable_classic(full_width = FALSE, html_font = "Times New Roman") %>%
  kable_styling(font_size = 9) %>%
  row_spec(0, bold = TRUE)

for (i in seq_along(tblD1_structure)) {
  tblD1_html <- pack_rows(tblD1_html, names(tblD1_structure)[i],
                          d1_starts[i], d1_ends[i],
                          label_row_css = d1_lcss)
}

tblD1_html <- tblD1_html %>%
  footnote(general = paste0(
    "Heteroskedasticity-robust (HC1) standard errors in parentheses. ",
    "*** p < 0.001, ** p < 0.01, * p < 0.05, . p < 0.10."),
    general_title = "Notes:", footnote_as_chunk = TRUE)

print(tblD1_html)

writeLines(
  as.character(tblD1_html),
  con = paper_path("Table_D1.html")
)

remove(tblD1_df, tblD1_html, profile_coef, profile_coef_list, prof_n,
       tblD1_structure, ordered_labels_d1, d1_sizes, d1_ends, d1_starts)
gc()

# TABLE D.2: Profile-specific AMEs
# Same compact layout as Table D.1 (profiles in columns).
d2_lcss <- "background-color: #f0f0f0; font-weight: bold;"

tblD2_df <- prf_ame %>%
  mutate(cell = sprintf("%.2f<br>[%.2f, %.2f]%s",
                        ame, ci_lo, ci_hi, ifelse(sig, "*", ""))) %>%
  select(profile, label, cell) %>%
  pivot_wider(names_from = profile, values_from = cell, names_prefix = "Profile ") %>%
  rename(Variable = label) %>%
  mutate(Variable = factor(Variable, levels = ordered_labels_d2)) %>%
  arrange(Variable) %>%
  mutate(Variable = as.character(Variable))

tblD2_html <- kbl(tblD2_df,
                  format  = "html",
                  escape  = FALSE,
                  caption = "Table D.2. Profile-specific AMEs",
                  align   = c("l", rep("c", ncol(tblD2_df) - 1))) %>%
  kable_classic(full_width = FALSE, html_font = "Times New Roman") %>%
  kable_styling(font_size = 9) %>%
  row_spec(0, bold = TRUE)

for (i in seq_along(tblD2_structure)) {
  tblD2_html <- pack_rows(tblD2_html, names(tblD2_structure)[i],
                          d2_starts[i], d2_ends[i],
                          label_row_css = d2_lcss)
}

tblD2_html <- tblD2_html %>%
  footnote(general = paste0("Survey weights trimmed at the 99.5th percentile. ",
                            "AME (percentage points) [bootstrapped 95% CI]. ",
                            "Continuous variables: effect of a one-standard-deviation ",
                            "increase, using the standard deviation within each profile. ",
                            "* = CI excludes zero."),
           general_title = "Note:", footnote_as_chunk = TRUE)

print(tblD2_html)

writeLines(
  as.character(tblD2_html),
  con = paper_path("Table_D2.html")
)

remove(canonical_labels, var_block,
       tblD2_df, tblD2_html, tblD2_structure, ordered_labels_d2,
       d2_sizes, d2_ends, d2_starts)
gc()


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE C.2: Observed vs predicted adoption by profile (calibration check)
# ══════════════════════════════════════════════════════════════════════════════

# Overall weighted adoption mean (dashed reference line)
overall_obs <- weighted.mean(data_clean$retrofit, w = data_clean$weight_final,
                             na.rm = TRUE)

# Short axis labels, defined with the profile labels in stage 03
figC2_df <- prof_adopt %>%
  transmute(profile_gmm,
            Observed  = obs,
            Predicted = pred) %>%
  pivot_longer(c(Observed, Predicted),
               names_to = "type", values_to = "rate") %>%
  mutate(type        = factor(type, levels = c("Observed", "Predicted")),
         profile_gmm = factor(profile_gmm, levels = names(profile_labels_short),
                              labels = unname(profile_labels_short)))

figC2 <- ggplot(figC2_df, aes(x = profile_gmm, y = rate, fill = type)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62,
           color = "gray20", linewidth = 0.3) +
  geom_hline(yintercept = overall_obs, linetype = "dashed",
             color = "gray40", linewidth = 0.4) +
  geom_text(aes(label = scales::percent(rate, accuracy = 0.1)),
            position = position_dodge(width = 0.7),
            vjust = -0.5, size = 2.8, color = "black") +
  scale_fill_manual(values = c("Observed" = "gray50", "Predicted" = "gray80"),
                    name = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, max(figC2_df$rate) * 1.15),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title   = NULL,
       x       = "Household profile (GMM)",
       y       = "Adoption rate (%)")  +
  theme_minimal(base_size = 10) +
  theme(
    plot.caption          = element_text(size = 8, color = "gray30",
                                         hjust = 0, margin = margin(t = 10)),
    plot.caption.position = "plot",
    axis.text.x           = element_text(size = 8, color = "black",
                                         lineheight = 0.95),
    axis.text.y           = element_text(size = 9, color = "black"),
    axis.title.x          = element_text(size = 10, margin = margin(t = 8)),
    axis.title.y          = element_text(size = 10, margin = margin(r = 6)),
    panel.background      = element_blank(),
    panel.grid.major.x    = element_blank(),
    panel.grid.major.y    = element_line(color = "gray85", linewidth = 0.3),
    panel.grid.minor      = element_blank(),
    axis.ticks.x          = element_line(color = "gray40", linewidth = 0.3),
    axis.ticks.y          = element_line(color = "gray40", linewidth = 0.3),
    axis.ticks.length     = unit(2, "pt"),
    legend.position       = c(0.99, 0.99),
    legend.justification  = c("right", "top"),
    legend.key.size       = unit(0.4, "cm"),
    legend.text           = element_text(size = 8),
    legend.background     = element_rect(fill = "white", color = "gray70", linewidth = 0.3),
    legend.margin         = margin(3, 4, 3, 4),
    plot.margin           = margin(t = 5, r = 10, b = 5, l = 5)
  )

print(figC2)
save_figure(figC2, "Figure_C2",
            w_pdf = 7.5, h_pdf = 4.2, w_tiff = 19, h_tiff = 10.5)

remove(figC2, figC2_df, prof_adopt, overall_obs)
gc()


# (The profile-faceted AME figure was removed: it is not used in the
# manuscript. Tables D.1 and D.2 above carry the profile-specific results.)


# ══════════════════════════════════════════════════════════════════════════════
# SCRIPT COMPLETE
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== SCRIPT COMPLETE ===\n")
cat(sprintf("Total runtime: %.1f minutes\n",
            as.numeric(difftime(Sys.time(), start_time, units = "mins"))))
cat(sprintf("Final memory usage: %.1f MB\n", sum(gc()[, 2])))
