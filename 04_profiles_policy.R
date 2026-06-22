###############################################################################
# Project       : [rehabSpain] [4] PROFILE-SPECIFIC ANALYSIS
# Creation date : 02/12/2024
# Last update   : 22/06/2026
# Author        : Mercè Amich (merce.amich@ehu.eus)
# Institution   : UPV/EHU, BC3
# Last run time : 5 min.
# Script Overview:
#   Profile-specific adoption analysis:
#     - adoption-rate heterogeneity across profiles (Pearson chi-square)
#     - effect heterogeneity: restricted vs profile-interacted logit (LR test)
#     - profile-specific bootstrapped AMEs with significance classification
#     - profile-faceted AME figure (Figure 3)
# Requirements:
#   - This script must be in the same directory as:
#       a) _setup.R
#       b) 03_famd_clustering.RData   (output of 03_famd_clustering.R)
# Output  : Publication-ready table (kableExtra) + ggplots (TIFF/PDF)
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

# ── Adoption rate heterogeneity ────
chi_adopt <- chisq.test(data_clean$profile_gmm, data_clean$retrofit)
cat(sprintf("  Profile × Retrofit: χ²(%d) = %.2f, p = %.4f\n",
            chi_adopt$parameter, chi_adopt$statistic, chi_adopt$p.value))

# ── Effect heterogeneity: restricted vs profile-interacted logit (LR test) ────
lr_terms <- c(all.vars(formula_baseline), "profile_gmm", "weight_final")
df_lr    <- data_clean[complete.cases(data_clean[, lr_terms]), ]
df_lr$profile_gmm <- factor(df_lr$profile_gmm)

# Restricted: common slopes, profile-specific intercepts
#   retrofit ~ active_vars + profile_gmm
form_restricted <- update(formula_baseline, . ~ . + profile_gmm)
m_restricted <- glm(form_restricted,
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
# value means an interaction term was dropped for collinearity (most likely
# ownership in P4), which would then need to be noted.
if (!is.na(lr_df) && lr_df != 36)
  warning(sprintf("LR test df = %d, expected 36 (an interaction term may be unidentified).", lr_df))

remove(df_lr, m_restricted, m_interacted, form_restricted, form_interacted, lr_terms)
gc()

# Get predicted probabilities from baseline model
data_clean$pred_prob <- predict(logit_baseline, type = "response")

# Calculate observed vs predicted adoption by profile
prof_adopt <- data_clean %>%
  group_by(profile_gmm) %>%
  summarise(n    = n(),
            obs  = mean(retrofit, na.rm = TRUE),
            pred = mean(pred_prob, na.rm = TRUE),
            .groups = "drop")
gc()

# ── Profile-specific AMEs ────
cat("  Profile-specific AMEs (1,000 bootstrap iterations per profile)...\n")
profiles <- sort(unique(data_clean$profile_gmm))
set.seed(123)
n_boot_p <- 1000
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
    
    a <- setNames(numeric(length(valid)), valid)
    for (v in valid) {
      uv  <- unique(d[[v]][!is.na(d[[v]])])
      bin <- length(uv) == 2 && all(sort(uv) %in% c(0, 1))
      
      if (bin) {
        # Binary variable: discrete first difference
        d0 <- d; d0[[v]] <- 0L
        d1 <- d; d1[[v]] <- 1L
        a[v] <- mean(predict(m, newdata = d1, type = "response") -
                       predict(m, newdata = d0, type = "response"), na.rm = TRUE)
      } else {
        # Continuous variable: marginal effect scaled by SD
        pp     <- predict(m, newdata = d, type = "response")
        sd_x   <- sd(d[[v]], na.rm = TRUE)
        a[v]   <- coef(m)[v] * mean(pp * (1 - pp), na.rm = TRUE) * sd_x
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
  mutate(sig = !(ci_lo < 0 & ci_hi > 0),
         label = var_labels[variable])

# TABLE S6: Profile-specific AMEs
# Build the data frame first
tblS6_df <- prf_ame %>%
  mutate(cell = sprintf("%.2f<br>[%.2f, %.2f]%s",
                        ame, ci_lo, ci_hi, ifelse(sig, "*", ""))) %>%
  select(profile, label, cell) %>%
  pivot_wider(names_from = profile, values_from = cell, names_prefix = "Profile ") %>%
  rename(Variable = label)

# Then build the kable object from the data frame
tblS6_html <- kbl(tblS6_df,
                  format  = "html",
                  escape  = FALSE,
                  caption = "Table S6. Profile-specific AMEs",
                  align   = c("l", rep("c", ncol(tblS6_df) - 1))) %>%
  kable_classic(full_width = FALSE, html_font = "Times New Roman") %>%
  kable_styling(font_size = 9) %>%
  row_spec(0, bold = TRUE) %>%
  footnote(general = "AME (percentage points) [bootstrapped 95% CI]. * = CI excludes zero.",
           general_title = "Note:", footnote_as_chunk = TRUE)

print(tblS6_html)

writeLines(
  as.character(tblS6_html),
  con = file.path(getwd(), "TableS6_Profile_AME.html")
)

remove(tblS6_df, tblS6_html)
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
  "highpopulation"   = "Socio-demographic",
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
  "High population municipality",
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
    caption  = "Notes: AMEs re-estimated within each profile subsample. Filled points indicate 95% CI excludes zero;\nhollow points indicate non-significant effects. Bootstrapped CIs (1,000 iterations per profile)."
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
       width = 7.5, height = 5.5, units = "in")

# FIGURE 3: Profile-specific AMEs (block-structured)

block_map_fig3 <- c(
  "CDD_1723"         = "A. Climatic",
  "HDD_1723"         = "A. Climatic",
  "ownership"        = "B. Dwelling and environmental",
  "detached"         = "B. Dwelling and environmental",
  "moisturedamage"   = "B. Dwelling and environmental",
  "rooms"            = "B. Dwelling and environmental",
  "pollution"        = "B. Dwelling and environmental",
  "ln_income"        = "C. Socio-demographic",
  "allhigheducation" = "C. Socio-demographic",
  "meanageadults"    = "C. Socio-demographic",
  "highpopulation"   = "C. Socio-demographic",
  "kids"             = "C. Socio-demographic"
)
block_order_fig3 <- c("A. Climatic",
                      "B. Dwelling and environmental",
                      "C. Socio-demographic")

pal_profiles_distinct <- c("#e74c3c", "#3498db", "#2ecc71", "#f39c12")
pal_profiles_named    <- setNames(pal_profiles_distinct, paste0("P", 1:4))
pal_with_gray         <- c("ns" = "gray75", pal_profiles_named)

prf_ame_plot <- prf_ame %>%
  mutate(
    block     = block_map_fig3[variable],
    block     = factor(block, levels = block_order_fig3)
  ) %>%
  group_by(block) %>%
  mutate(label = fct_reorder(factor(label), ame, .fun = mean)) %>%
  ungroup() %>%
  mutate(
    prf       = factor(paste0("P", profile)),
    color_sig = ifelse(sig, paste0("P", profile), "ns"),
    color_sig = factor(color_sig, levels = c("ns", paste0("P", 1:4)))
  )

fig3 <- ggplot(prf_ame_plot,
               aes(x = label, y = ame, color = color_sig, group = prf)) +
  
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "firebrick", linewidth = 0.4) +
  
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                position = position_dodge(0.7),
                width = 0.2, linewidth = 0.35) +
  
  geom_point(position = position_dodge(0.7),
             size = 1.8, stroke = 0.5) +
  
  scale_color_manual(
    values = pal_with_gray,
    breaks = paste0("P", 1:4),
    labels = paste0("P", 1:4),
    name   = "Profile"
  ) +
  
  facet_grid(block ~ ., scales = "free_y", space = "free_y") +
  
  coord_flip() +
  
  scale_y_continuous(
    breaks = seq(-6, 16, by = 2),
    labels = number_format(accuracy = 0.1),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  
  labs(
    title    = "Figure 3. Profile-specific average marginal effects on retrofit adoption",
    subtitle = "FAMD 4 dims / GMM k = 4. Bootstrapped 95% CIs (1,000 iterations per profile).",
    x        = NULL,
    y        = "AME (percentage points)",
    caption = paste0(
      "Notes:\n",
      "Re-estimated within each profile subsample.\n",
      "For continuous variables, the AME represents the effect of a one standard deviation increase.\n",
      "For binary variables, the AME represents the average difference ",
      "in predicted probability between categories.\n",
      "Gray = not significant (CI includes zero).\n",
      "Survey weights applied throughout."
    )
  ) +
  
  theme_minimal(base_size = 11, base_family = "Times New Roman") +
  theme(
    plot.title.position = "plot",
    plot.title          = element_text(face = "bold", size = 12,
                                       hjust = 0, margin = margin(b = 4)),
    plot.subtitle       = element_text(size = 10, hjust = 0,
                                       color = "grey30", margin = margin(b = 10)),
    plot.caption        = element_text(size = 8, color = "grey30",
                                       hjust = 0, margin = margin(t = 10)),
    plot.caption.position = "plot",
    axis.text.y         = element_text(size = 9,   color = "black"),
    axis.text.x         = element_text(size = 8.5, color = "black"),
    axis.title.x        = element_text(size = 10,  margin = margin(t = 8)),
    panel.spacing.y     = unit(0.8, "lines"),
    panel.background    = element_rect(fill = "grey98", color = NA),
    panel.grid.major.x  = element_line(color = "grey90", linewidth = 0.25),
    panel.grid.minor.x  = element_blank(),
    panel.grid.major.y  = element_blank(),
    strip.text.y        = element_text(face = "bold", size = 9, angle = 90,
                                       hjust = 0.5, vjust = 0.5, color = "black"),
    strip.background    = element_rect(fill = "grey93", color = NA),
    axis.ticks.y        = element_blank(),
    legend.position     = "right",
    legend.direction    = "vertical",
    legend.title        = element_text(face = "bold", size = 9),
    legend.text         = element_text(size = 9),
    plot.margin         = margin(t = 15, r = 10, b = 15, l = 10)
  ) +
  
  guides(color = guide_legend(
    override.aes = list(size = 3, shape = 19, linetype = 0)
  ))

print(fig3)
ggsave("Figure3_Profile_AME_FAMD4_k4.tiff",
       plot = fig3, device = "tiff",
       width = 20, height = 20, units = "cm", dpi = 900)

remove(prf_ame_plot, fig3, pal_profiles_distinct, pal_profiles_named,
       pal_with_gray, block_map_fig3, block_order_fig3)
gc()


# ══════════════════════════════════════════════════════════════════════════════
# SCRIPT COMPLETE
# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== SCRIPT COMPLETE ===\n")
cat(sprintf("Total runtime: %.1f minutes\n",
            as.numeric(difftime(Sys.time(), start_time, units = "mins"))))
cat(sprintf("Final memory usage: %.1f MB\n", sum(gc()[, 2])))
