###############################################################################
# Project       : [rehabSpain] [4] PROFILE-SPECIFIC ANALYSIS
# Creation date : 02/12/2024
# Last update   : 26/06/2026
# Author        : Mercè Amich (merce.amich@ehu.eus)
# Institution   : UPV/EHU, BC3
# Last run time : 5 min.

# Script Overview:
#   Profile-specific adoption analysis:
#     - effect heterogeneity: restricted vs profile-interacted logit (LR test)
#     - profile-specific HC1-robust logit coefficients (Table D.1)
#     - profile-specific bootstrapped AMEs (Table D.2)
#     - observed vs predicted adoption by profile (Figure C.2)
#     - profile-faceted AME figure (Figure 3)

# Requirements:
#   - This script must be in the same directory as:
#       a) _setup.R
#       b) 03_famd_clustering.RData   (output of 03_famd_clustering.R)

# Output  : Publication-ready tables (kableExtra) + ggplots (TIFF/PDF)

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
#   data_clean (with profile_gmm / profile_uncertainty), logit_baseline,
#   formula_baseline, ame_df, active_vars, var_labels, outcome_var,
#   profile_labels, famd_result
load("03_famd_clustering.RData")

# ══════════════════════════════════════════════════════════════════════════════
# 7. SECTION 4.4: PROFILE-SPECIFIC ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════

# ── Effect heterogeneity: restricted vs profile-interacted logit (LR test) ────
lr_terms          <- c(all.vars(formula_baseline), "profile_gmm", "weight_final")
df_lr             <- data_clean[complete.cases(data_clean[, lr_terms]), ]
df_lr$profile_gmm <- factor(df_lr$profile_gmm)

# Restricted: common slopes, profile-specific intercepts
#   retrofit ~ active_vars + profile_gmm
form_restricted <- update(formula_baseline, . ~ . + profile_gmm)
m_restricted    <- glm(form_restricted,
                    family  = binomial("logit"),
                    data    = df_lr,
                    weights = weight_final)

# Interacted: retrofit ~ (active_vars) * profile_gmm
#   (profile-specific intercepts and slopes). Interaction terms not identified
#   within a profile (e.g. ownership in the near-owner-free renter profile) are
#   dropped by glm; the test df reflects identified parameters only.
form_interacted <- as.formula(
  paste("retrofit ~ (", paste(active_vars, collapse = " + "), ") * profile_gmm")
)
m_interacted <- glm(form_interacted,
                    family  = binomial("logit"),
                    data    = df_lr,
                    weights = weight_final)

lr_effect_het <- anova(m_restricted, m_interacted, test = "LRT")
print(lr_effect_het)

lr_chisq <- lr_effect_het$Deviance[2]
lr_df    <- lr_effect_het$Df[2]
lr_p     <- lr_effect_het$`Pr(>Chi)`[2]

cat(sprintf("  Effect heterogeneity (restricted vs interacted): LR χ²(%d) = %.2f, p = %.4g\n",
            lr_df, lr_chisq, lr_p))

# Sanity check: df should be 36 (12 predictors x 3 profile contrasts). A lower
# value means an interaction term was dropped for collinearity
if (!is.na(lr_df) && lr_df != 36)
  warning(sprintf("LR test df = %d, expected 36 (an interaction term may be unidentified).", lr_df))

remove(df_lr, m_restricted, m_interacted, form_restricted, form_interacted, lr_terms)
gc()

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
  m  <- glm(formula_baseline, family = binomial("logit"),
            data = dp, weights = weight_final)
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
set.seed(123)
n_boot_p     <- 1000
prf_ame_list <- lapply(profiles, function(p) {
  cat(sprintf("    Profile %d/%d\n", p, max(profiles)))
  
  dp <- data_clean %>% filter(profile_gmm == p)
  
  valid <- active_vars[sapply(active_vars, function(v)
    length(unique(dp[[v]][!is.na(dp[[v]])])) >= 2)]
  
  if (!length(valid)) return(NULL)
  
  bfn <- function(data, idx) {
    d <- as.data.frame(data[idx, ])
    m <- tryCatch(
      glm(formula_baseline, family = binomial("logit"),
          data = d, weights = weight_final),
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
        # Continuous variable: marginal effect scaled by SD
        pp     <- predict(m, newdata = d, type = "response")
        sd_x   <- sd(d[[v]], na.rm = TRUE)
        a[v]   <- coef(m)[v] * weighted.mean(pp * (1 - pp), w = w, na.rm = TRUE) * sd_x
      }
    }
    a * 100
  }
  
  br <- boot(dp, bfn, R = n_boot_p)
  data.frame(profile = p, variable = valid,
             ame   = apply(br$t, 2, mean, na.rm = TRUE),
             se    = apply(br$t, 2, sd, na.rm = TRUE),
             ci_lo = apply(br$t, 2, quantile, probs = 0.025, na.rm = TRUE),
             ci_hi = apply(br$t, 2, quantile, probs = 0.975, na.rm = TRUE))
})

prf_ame <- bind_rows(prf_ame_list) %>%
  mutate(sig   = !(ci_lo < 0 & ci_hi > 0),
         label = var_labels[variable])

# ── Row structure for Table D.1 (HC1 coefficients): Intercept + A/B/C ─────────

tblD1_structure <- list(
  "Intercept" = c("Intercept"),
  "A. Climatic variables" = c(
    "Cooling Degree Days (mean '17-'23)",
    "Heating Degree Days (mean '17-'23)"
  ),
  "B. Dwelling and environmental characteristics" = c(
    "Ownership",
    "Detached dwelling",
    "Moisture / structural damage",
    "Number of rooms",
    "Polluted environment"
  ),
  "C. Socio-demographic characteristics" = c(
    "Income (log)",
    "High educated",
    "Mean age of adults",
    "Highly urbanised municipality",
    "Children present"
  )
)

ordered_labels_d1 <- unlist(tblD1_structure, use.names = FALSE)

d1_sizes  <- lengths(tblD1_structure)
d1_ends   <- cumsum(d1_sizes)
d1_starts <- d1_ends - d1_sizes + 1L

# ── Row structure for Table D.2 (AMEs): A/B/C only ────────────────────────────
tblD2_structure <- list(
  "A. Climatic variables" = c(
    "Cooling Degree Days (mean '17-'23)",
    "Heating Degree Days (mean '17-'23)"
  ),
  "B. Dwelling and environmental characteristics" = c(
    "Ownership",
    "Detached dwelling",
    "Moisture / structural damage",
    "Number of rooms",
    "Polluted environment"
  ),
  "C. Socio-demographic characteristics" = c(
    "Income (log)",
    "High educated",
    "Mean age of adults",
    "Highly urbanised municipality",
    "Children present"
  )
)

ordered_labels_d2 <- unlist(tblD2_structure, use.names = FALSE)

d2_sizes  <- lengths(tblD2_structure)
d2_ends   <- cumsum(d2_sizes)
d2_starts <- d2_ends - d2_sizes + 1L

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
  mutate(across(-Variable, ~ ifelse(is.na(.) | grepl("NA", .), "—", .)))

tblD1_html <- kbl(tblD1_df,
                  format  = "html",
                  escape  = FALSE,
                  caption = "Table D1. Profile-specific binary logit coefficient estimates with heteroskedasticity-robust standard errors (HC1)",
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
  footnote(general = sprintf(paste0(
    "N = %s. Maximum likelihood estimation with survey weights capped at the 99.5th percentile, ",
    "applied separately within each profile subsample. Coefficients are log-odds; continuous ",
    "coefficients are per unit of the predictor and are not SD-scaled, consistent with Table B.1. ",
    "Heteroskedasticity-robust (HC1) standard errors in parentheses. ",
    "*** p < 0.001, ** p < 0.01, * p < 0.05, . p < 0.10."),
    paste(sprintf("%s (P%d)", format(prof_n$n, big.mark = ","), prof_n$profile),
          collapse = ", ")),
    general_title = "Notes:", footnote_as_chunk = TRUE)

print(tblD1_html)

writeLines(
  as.character(tblD1_html),
  con = file.path(getwd(), "TableD1_Profile_coefficients.html")
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
                  caption = "Table D2. Profile-specific AMEs",
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
  footnote(general = paste0("Survey weights capped at the 99.5th percentile. ",
                            "AME (percentage points) [bootstrapped 95% CI]. ",
                            "* = CI excludes zero."),
           general_title = "Note:", footnote_as_chunk = TRUE)

print(tblD2_html)

writeLines(
  as.character(tblD2_html),
  con = file.path(getwd(), "TableD2_Profile_AME.html")
)

remove(tblD2_df, tblD2_html, tblD2_structure, ordered_labels_d2,
       d2_sizes, d2_ends, d2_starts)
gc()


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE C.2: Observed vs predicted adoption by profile (validation)
# ══════════════════════════════════════════════════════════════════════════════

# Overall weighted adoption mean (dashed reference line)
overall_obs <- weighted.mean(data_clean$retrofit, w = data_clean$weight_final,
                             na.rm = TRUE)

profile_x_labels <- c(
  "1" = "Cold-climate\nsemi-rural homeowners\nin large detached\nhomes (P1)",
  "2" = "Educated urban\nhomeowners in dense\napartment blocks (P2)",
  "3" = "Warm-climate\nhomeowners in\ndetached homes (P3)",
  "4" = "Low-income young\nrenters in small,\ndeteriorated urban\ndwellings (P4)"
)

figC2_df <- prof_adopt %>%
  transmute(profile_gmm,
            Observed  = obs,
            Predicted = pred) %>%
  pivot_longer(c(Observed, Predicted),
               names_to = "type", values_to = "rate") %>%
  mutate(type        = factor(type, levels = c("Observed", "Predicted")),
         profile_gmm = factor(profile_gmm, levels = 1:4,
                              labels = profile_x_labels[as.character(1:4)]))

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
       y       = "Adoption rate (%)",
       caption = sprintf(paste0(
         "Notes: Predicted rates from baseline logit. Dashed line = overall weighted mean (%.1f%%)"),
         100 * overall_obs, chi_adopt$parameter, chi_adopt$statistic)) +
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
ggsave("FigureC2_observed_vs_predicted.pdf", figC2,
       width = 7.5, height = 4.2, device = cairo_pdf)
ggsave("FigureC2_observed_vs_predicted.tiff",
       plot   = figC2,
       device = "tiff",
       width  = 19,
       height = 10.5,
       units  = "cm",
       dpi    = 900)

remove(figC2, figC2_df, prof_adopt, overall_obs, profile_x_labels)
gc()


# FIGURE 3: Profile-specific AMEs (faceted by profile)

block_map_fig3 <- c(
  "CDD_1723"         = "Climatic",
  "HDD_1723"         = "Climatic",
  "ownership"        = "Dwelling and environmental",
  "detached"         = "Dwelling and environmental",
  "moisturedamage"   = "Dwelling and environmental",
  "rooms"            = "Dwelling and environmental",
  "pollution"        = "Dwelling and environmental",
  "ln_income"        = "Socio-demographic",
  "allhigheducation" = "Socio-demographic",
  "meanageadults"    = "Socio-demographic",
  "highlyurbanised"  = "Socio-demographic",
  "kids"             = "Socio-demographic"
)
block_order_fig3 <- c("Climatic",
                      "Dwelling and environmental",
                      "Socio-demographic")

# Okabe-Ito colorblind-safe palette
pal_profiles_named <- c(
  "P1" = "#0072B2",
  "P2" = "#D55E00",
  "P3" = "#009E73",
  "P4" = "#CC79A7"
)

# Profile labels for facet headers
profile_facet_labels <- c(
  "P1" = "P1: Cold-climate semi-rural\nhomeowners (detached)",
  "P2" = "P2: Educated urban homeowners\n(apartment blocks)",
  "P3" = "P3: Warm-climate homeowners\n(detached)",
  "P4" = "P4: Low-income young renters\n(small, deteriorated urban)"
)

# Derive variable order from Figure 1's ame_df (overall AME magnitude),
# pooled across blocks - same ordering logic as Figure 1
var_order_fig1 <- ame_df %>%
  mutate(label = var_labels[variable]) %>%
  arrange(estimate) %>%
  pull(label)

var_order_fig1 <- rev(c(
  "Heating Degree Days (mean '17-'23)",
  "Cooling Degree Days (mean '17-'23)",
  "Ownership",
  "Polluted environment",
  "Detached dwelling",
  "Number of rooms",
  "Moisture / structural damage",
  "High educated",
  "Income (log)",
  "Highly urbanised municipality",
  "Children present",
  "Mean age of adults"
))

prf_ame_plot <- prf_ame %>%
  mutate(
    block = block_map_fig3[variable],
    block = factor(block, levels = block_order_fig3),
    label = factor(label, levels = var_order_fig1),
    prf   = factor(paste0("P", profile), levels = paste0("P", 1:4)),
    # Fill color: profile color if significant, white if not
    point_fill = ifelse(sig, as.character(pal_profiles_named[prf]), "white")
  )

fig3 <- ggplot(prf_ame_plot,
               aes(x = label, y = ame, color = prf)) +
  
  geom_hline(yintercept = 0, linetype = "solid",
             color = "gray40", linewidth = 0.3) +
  
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi, alpha = sig),
                width = 0, linewidth = 0.5) +
  
  geom_point(aes(fill = point_fill),
             shape = 21, size = 2.2, stroke = 0.7) +
  
  scale_color_manual(values = pal_profiles_named, guide = "none") +
  scale_fill_identity(guide = "none") +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5), guide = "none") +
  
  facet_wrap(~ prf, ncol = 4, labeller = labeller(prf = profile_facet_labels)) +
  
  coord_flip() +
  
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 5),
    labels = number_format(accuracy = 1),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  
  labs(
    title    = NULL,
    x        = NULL,
    y        = "Average marginal effect (pp)",
    caption  = "Notes: AMEs re-estimated within each profile subsample. 
    Filled points indicate 95% CI excludes zero;\nhollow points indicate non-significant effects. 
    Bootstrapped CIs (1,000 iterations per profile)."
  ) +
  
  theme_minimal(base_size = 10) +
  theme(
    plot.caption          = element_text(size = 8, color = "gray30",
                                         hjust = 0, margin = margin(t = 10)),
    plot.caption.position = "plot",
    axis.text.y           = element_text(size = 9, color = "black"),
    axis.text.x           = element_text(size = 8.5, color = "black"),
    axis.title.x          = element_text(size = 10, margin = margin(t = 6)),
    panel.spacing.x       = unit(0.9, "lines"),
    panel.background      = element_blank(),
    panel.grid.major.x    = element_line(color = "gray88", linewidth = 0.3),
    panel.grid.minor.x    = element_blank(),
    panel.grid.major.y    = element_line(color = "gray95", linewidth = 0.2),
    strip.text            = element_text(face = "bold", size = 8.5,
                                         color = "black", lineheight = 0.95,
                                         margin = margin(t = 4, b = 8)),
    strip.background      = element_blank(),
    strip.clip            = "off",
    axis.ticks.x          = element_line(color = "gray40", linewidth = 0.3),
    axis.ticks.length.x   = unit(2, "pt"),
    axis.ticks.y          = element_blank(),
    plot.margin           = margin(t = 10, r = 12, b = 5, l = 5)
  )

print(fig3)
ggsave("Figure3_Profile_AME_FAMD4_k4.pdf",
       plot = fig3, device = cairo_pdf,
       width = 10.5, height = 5.5, units = "in")


# ══════════════════════════════════════════════════════════════════════════════
# SCRIPT COMPLETE
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== SCRIPT COMPLETE ===\n")
cat(sprintf("Total runtime: %.1f minutes\n",
            as.numeric(difftime(Sys.time(), start_time, units = "mins"))))
cat(sprintf("Final memory usage: %.1f MB\n", sum(gc()[, 2])))
