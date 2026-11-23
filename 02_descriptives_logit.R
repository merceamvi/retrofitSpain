###############################################################################
# Project       : [rehabSpain] [2] DESCRIPTIVES & BASELINE LOGIT
# Creation date : 02/12/2024
# Last update   : 23/11/2026 (captions, labels and notes of Tables 1, A.1, B.1 and
#                 B.2 as in the paper; education label "Highly educated"; Figure
#                 B.1 labels as in Table B.2; AMEs of continuous variables scaled
#                 by the weighted standard deviation; glm warnings silenced)
# Author        : Mercè Amich (merce.amich@ehu.eus)
# Institution   : UPV/EHU, BC3
# Last run time : ~11 min.

# Script Overview:
#   Data preparation, outcome distribution and descriptive statistics,
#   baseline weighted binary logit, HC1-robust estimates and bootstrapped
#   average marginal effects (AMEs).

# Requirements:
#   - This script must be in the same directory as:
#       a) _setup.R
#       b) data.RData   (output of 01_build_data.R)

# Output  : paper/ Table_1, Table_A1, Table_B1, Table_B2 (HTML) and
#           Figure_B1 (PDF, TIFF)
#           + 02_descriptives_logit.RData (objects for the next stage)

###############################################################################

# ══════════════════════════════════════════════════════════════════════════════
# 0. SETUP
# ══════════════════════════════════════════════════════════════════════════════

rm(list = ls(all = TRUE))
Sys.setenv(LANG = "en")
start_time <- Sys.time()

# Shared packages, working directory, theme, palettes and save helpers
source("_setup.R")

# Load data.Rdata (contains 'data' and 'sample_info' objects)
load("data.RData")

cat(sprintf("After Script 1 cleaning (survey sample): n = %s\n", 
            format(nrow(data), big.mark = ",")))
cat(sprintf("After Script 1 cleaning (weighted pop):  N = %s\n", 
            format(round(sum(data$weight)), big.mark = ",")))

# ══════════════════════════════════════════════════════════════════════════════
# 1. DATA PREPARATION
# ══════════════════════════════════════════════════════════════════════════════

# Define active variables and labels

active_vars <- c(
  "CDD_1723", "HDD_1723",  # Cooling and Heating Degree Days (mean 2017-2023)
  "ownership",             # Ownership
  "detached",              # Detached dwelling
  "moisturedamage",        # Self-reported moisture damage in dwelling
  "rooms",                 # Number of rooms
  "pollution",             # Self-reported pollution exposure
  "ln_income",             # Income (log)
  "allhigheducation",      # All members have high education
  "meanageadults",         # Mean age adults (>= 18 years old)
  "highlyurbanised",       # Highly urbanised municipality
  "kids"                   # Presence of children in the household
)

var_labels <- c(
  "CDD_1723"         = "Cooling Degree Days (mean '17-'23)",
  "HDD_1723"         = "Heating Degree Days (mean '17-'23)",
  "ownership"        = "Ownership",
  "detached"         = "Detached dwelling",
  "moisturedamage"   = "Moisture / structural damage",
  "rooms"            = "Number of rooms",
  "pollution"        = "Polluted environment",
  "ln_income"        = "Income (log)",
  "allhigheducation" = "Highly educated",
  "meanageadults"    = "Mean age of adults",
  "highlyurbanised"  = "Highly urbanised municipality",
  "kids"             = "Children present"
)

outcome_var   <- "retrofit" # Dummy (0, 1)
required_vars <- c(outcome_var, active_vars, "weight", "region") # Intersect

# Data cleaning: listwise deletion on required variables
n_total    <- nrow(data)
data_clean <- data[complete.cases(data[, required_vars]), ] # Keep if no NA's
n_removed  <- n_total - nrow(data_clean) # Keep track of removed observations

# Report of removed observations (both in this and in the first script)
cat(sprintf("  Removed %d obs (%.1f%%) due to missing values\n  
  Total removed in 1st script:   %.2f%% (survey), %.2f%% (weighted)\n  
  Total removed in both scripts: %.2f%% (survey), %.2f%% (weighted)\n  
  Final sample: n = %s\n",
            n_removed, 
            100 * n_removed / n_total,
            sample_info$pct_removed_cumulative_unweighted,
            sample_info$pct_removed_cumulative_weighted,
            100 * (sample_info$n_original_unweighted - nrow(data_clean)) / sample_info$n_original_unweighted,
            100 * (sample_info$n_original_weighted - sum(data_clean$weight)) / sample_info$n_original_weighted,
            format(nrow(data_clean), big.mark = ",")))

# Normalise and trim survey weights at 99.5th percentile
w_raw                   <- data_clean$weight
w_norm                  <- w_raw / mean(w_raw, na.rm = TRUE)
w_cap                   <- pmin(w_norm, quantile(w_norm, 0.995, na.rm = TRUE))
data_clean$weight_final <- w_cap / mean(w_cap)

# Percentage of adopters (retrofit == 1) in the final sample
cat(sprintf("  Adopters: %d (%.1f%%)\n",
            sum(data_clean$retrofit), 100 * mean(data_clean$retrofit)))

# Clean up
remove(data, w_raw, w_norm, w_cap, n_total, n_removed)
gc()


# TABLE 1: DISTRIBUTION OF THE OUTCOME VARIABLE (original categories)

#  Compute row-level stats for each original category
outcome_dist_raw <- data_clean %>%
  mutate(
    count_label = case_when(
      count == 0 ~ "No measure",
      count == 1 ~ "One measure",
      count == 2 ~ "Two measures",
      count == 3 ~ "Three or more measures",
      TRUE       ~ NA_character_
    ),
    count_label = factor(count_label,
                         levels = c("No measure",
                                    "One measure",
                                    "Two measures",
                                    "Three or more measures"))
  ) %>%
  filter(!is.na(count_label)) %>%
  group_by(count_label) %>%
  summarise(
    n_unweighted = n(),
    n_weighted   = sum(weight_final, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pct_unweighted = 100 * n_unweighted / sum(n_unweighted),
    pct_weighted   = 100 * n_weighted   / sum(n_weighted)
  )

# Build summary rows for binary groups
non_adopter_summary <- outcome_dist_raw %>%
  filter(count_label == "No measure") %>%
  summarise(
    Category       = "Non-adopter",
    Outcome        = "0",
    n_unweighted   = sum(n_unweighted),
    pct_unweighted = sum(pct_unweighted),
    pct_weighted   = sum(pct_weighted)
  )

adopter_summary <- outcome_dist_raw %>%
  filter(count_label != "No measure") %>%
  summarise(
    Category       = "Adopter",
    Outcome        = "1",
    n_unweighted   = sum(n_unweighted),
    pct_unweighted = sum(pct_unweighted),
    pct_weighted   = sum(pct_weighted)
  )

# Build disaggregated adopter rows (indented labels)
adopter_detail <- outcome_dist_raw %>%
  filter(count_label != "No measure") %>%
  mutate(
    Category  = paste0("\u00a0\u00a0\u00a0\u00a0", as.character(count_label)),
    Outcome   = ""
  ) %>%
  select(Category, Outcome, n_unweighted, pct_unweighted, pct_weighted)

# Combine into final table
tbl1 <- bind_rows(
  non_adopter_summary,
  adopter_summary,
  adopter_detail
) %>%
  mutate(
    `n (unweighted)`  = format(round(n_unweighted), big.mark = ","),
    `% (unweighted)`  = sprintf("%.1f%%", pct_unweighted),
    `% (weighted)`    = sprintf("%.1f%%", pct_weighted)
  ) %>%
  select(Category, Outcome,
         `n (unweighted)`, `% (unweighted)`, `% (weighted)`)

# Row indices for pack_rows
#   Row 1: non-adopter summary
#   Row 2: adopter summary
#   Rows 3-5: adopter detail (three sub-categories)

kbl1 <- kbl(
  tbl1,
  format    = "html",
  caption   = "Table 1. Distribution of the number of measures and binary outcome definition",
  align     = c("l", "c", "r", "r", "r"),
  row.names = FALSE,
  col.names = c("Category", "Outcome",
                "N", "% (unweighted)", "% (weighted)")
) %>%
  kable_classic(full_width = FALSE, html_font = "Times New Roman") %>%
  kable_styling(font_size = 10) %>%
  row_spec(0, bold = TRUE) %>%
  
  # Bold the two summary rows
  row_spec(1, bold = TRUE) %>%
  row_spec(2, bold = TRUE) %>%
  
  # Light background on detail rows to signal subordination
  row_spec(3:5, background = "#f7f7f7", italic = TRUE) %>%
  footnote(
    general = sprintf(
      "N = %s households. Survey weights trimmed at 99.5th percentile.",
      format(nrow(data_clean), big.mark = ",")
    ),
    general_title     = "Notes:",
    footnote_as_chunk = TRUE
  )

print(kbl1)
save_table(kbl1, "Table_1.html")

# TABLE A.1: DESCRIPTIVE STATISTICS

# Variable descriptions (expanded)
var_descriptions <- c(
  "CDD_1723"         = "CDD (mean 2017-2023)",
  "HDD_1723"         = "HDD (mean 2017-2023)",
  "ownership"        = "Homeownership",
  "detached"         = "Detached dwelling",
  "moisturedamage"   = "Dwelling with moisture or structural damage",
  "rooms"            = "Number of rooms",
  "pollution"        = "Polluted environment",
  "ln_income"        = "Income (log)",
  "allhigheducation" = "All members have tertiary education",
  "meanageadults"    = "Mean age of adults (≥18y)",
  "highlyurbanised"  = "Highly urbanised municipality",
  "kids"             = "Children present"
)

# Variable types
continuous_vars <- c("CDD_1723", "HDD_1723", "rooms", "ln_income", "meanageadults")
binary_vars     <- c("ownership", "detached", "moisturedamage", "pollution", 
                     "allhigheducation", "highlyurbanised", "kids")

# Survey weights used throughout the descriptive block. All descriptive
# statistics are weighted, for consistency with the estimation, the profile
# characterisation in Table 3 and the decomposition. Minima and maxima are the
# exception: they are sample extremes and weighting them has no meaning.
w_desc <- data_clean$weight_final

# Outcome variable
outcome_stats <- data.frame(
  Variable    = "retrofit",
  Description = "Retrofit adoption",
  Type        = "Binary",
  Mean        = sprintf("%.2f", weighted.mean(data_clean$retrofit, w = w_desc, na.rm = TRUE)),
  SD          = "-",
  Min         = "-",
  Median      = "-",
  Max         = "-",
  `Yes (%)`   = sprintf("%.1f%%", 100 * weighted.mean(data_clean$retrofit, w = w_desc, na.rm = TRUE)),
  `No (%)`    = sprintf("%.1f%%", 100 * (1 - weighted.mean(data_clean$retrofit, w = w_desc, na.rm = TRUE))),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Function to compute stats for continuous variables (weighted)
compute_continuous <- function(var_name) {
  x <- data_clean[[var_name]]
  data.frame(
    Variable    = var_name,
    Description = var_descriptions[var_name],
    Type        = "Continuous",
    Mean        = sprintf("%.2f", weighted.mean(x, w = w_desc, na.rm = TRUE)),
    SD          = sprintf("%.2f", sqrt(Hmisc::wtd.var(x, weights = w_desc,
                                                      normwt = TRUE, na.rm = TRUE))),
    Min         = sprintf("%.2f", min(x, na.rm = TRUE)),
    Median      = sprintf("%.2f", Hmisc::wtd.quantile(x, weights = w_desc,
                                                      probs = 0.5, normwt = TRUE)),
    Max         = sprintf("%.2f", max(x, na.rm = TRUE)),
    `Yes (%)`   = "-",
    `No (%)`    = "-",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

# Function to compute stats for binary variables (weighted)
compute_binary <- function(var_name) {
  x <- data_clean[[var_name]]
  pct_yes <- weighted.mean(x, w = w_desc, na.rm = TRUE)
  data.frame(
    Variable    = var_name,
    Description = var_descriptions[var_name],
    Type        = "Binary",
    Mean        = sprintf("%.2f", pct_yes),
    SD          = "-",
    Min         = "-",
    Median      = "-",
    Max         = "-",
    `Yes (%)`   = sprintf("%.1f%%", 100 * pct_yes),
    `No (%)`    = sprintf("%.1f%%", 100 * (1 - pct_yes)),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

# Compute statistics for each variable group.
# The six blocks match the organisation of Section 3 and the blocks used in the
# Shapley attribution, so that the descriptive table, the results tables and the
# decomposition all group the variables in the same way.

# A. DEPENDENT VARIABLE
desc_outcome <- outcome_stats

# B. TENURE
desc_tenure <- compute_binary("ownership")

# C. INCOME
desc_income <- compute_continuous("ln_income")

# D. EDUCATION
desc_education <- compute_binary("allhigheducation")

# E. CLIMATIC EXPOSURE
desc_climate <- rbind(
  compute_continuous("HDD_1723"),
  compute_continuous("CDD_1723")
)

# F. DWELLING CHARACTERISTICS
desc_dwelling <- rbind(
  compute_continuous("rooms"),
  compute_binary("detached"),
  compute_binary("moisturedamage"),
  compute_binary("pollution")
)

# G. HOUSEHOLD COMPOSITION AND LOCATION
desc_household <- rbind(
  compute_continuous("meanageadults"),
  compute_binary("kids"),
  compute_binary("highlyurbanised")
)

# Combine all sections, in the order used throughout the paper
desc_stats_full <- rbind(
  desc_outcome, desc_tenure, desc_income, desc_education,
  desc_climate, desc_dwelling, desc_household
)

# Cumulative end position of each block, computed from the block sizes so that
# adding or removing a variable does not require renumbering by hand
end_outcome   <- nrow(desc_outcome)
end_tenure    <- end_outcome   + nrow(desc_tenure)
end_income    <- end_tenure    + nrow(desc_income)
end_education <- end_income    + nrow(desc_education)
end_climate   <- end_education + nrow(desc_climate)
end_dwelling  <- end_climate   + nrow(desc_dwelling)
end_household <- end_dwelling  + nrow(desc_household)

# Create kableExtra table with section headers
kblA1 <- kbl(
  desc_stats_full %>% select(-Variable),
  format    = "html",
  caption   = "Table A.1. Descriptive statistics",
  align     = c("l", "c", "r", "r", "r", "r", "r", "r", "r", "r"),
  row.names = FALSE,
  col.names = c("Description", "Type", "Mean", "SD", "Min", "Median", "Max",
                "Yes", "No"),
  digits    = 2
) %>%
  kable_classic(full_width = FALSE, html_font = "Times New Roman") %>%
  kable_styling(font_size = 10) %>%
  row_spec(0, bold = TRUE) %>%
  
  pack_rows("A. Dependent variable", 1, end_outcome,
            label_row_css = "background-color: #f0f0f0; font-weight: bold;") %>%
  pack_rows("B. Tenure", end_outcome + 1, end_tenure,
            label_row_css = "background-color: #f0f0f0; font-weight: bold;") %>%
  pack_rows("C. Income", end_tenure + 1, end_income,
            label_row_css = "background-color: #f0f0f0; font-weight: bold;") %>%
  pack_rows("D. Education", end_income + 1, end_education,
            label_row_css = "background-color: #f0f0f0; font-weight: bold;") %>%
  pack_rows("E. Climatic exposure", end_education + 1, end_climate,
            label_row_css = "background-color: #f0f0f0; font-weight: bold;") %>%
  pack_rows("F. Dwelling characteristics", end_climate + 1, end_dwelling,
            label_row_css = "background-color: #f0f0f0; font-weight: bold;") %>%
  pack_rows("G. Household composition and location", end_dwelling + 1, end_household,
            label_row_css = "background-color: #f0f0f0; font-weight: bold;") %>%
  
  footnote(
    general = sprintf(paste("N = %s households. Binary variables coded as",
                            "1 = Yes, 0 = No. Means, standard deviations,",
                            "medians and percentages are weighted using survey",
                            "weights trimmed at the 99.5th percentile.",
                            "Minima and maxima are sample extremes."),
                      format(nrow(data_clean), big.mark = ",")),
    general_title = "Notes:",
    footnote_as_chunk = TRUE
  )

print(kblA1)
save_table(kblA1, "Table_A1.html")

# Clean up
remove(w_desc,
       desc_outcome, desc_tenure, desc_income, desc_education,
       desc_climate, desc_dwelling, desc_household, desc_stats_full,
       end_outcome, end_tenure, end_income, end_education,
       end_climate, end_dwelling, end_household,
       continuous_vars, binary_vars, var_descriptions, outcome_stats,
       compute_continuous, compute_binary, kbl1, kblA1, sample_info,
       adopter_detail, adopter_summary, non_adopter_summary, outcome_dist_raw)
gc()


# ══════════════════════════════════════════════════════════════════════════════
# 2. BASELINE BINARY LOGIT ADOPTION MODEL
# ══════════════════════════════════════════════════════════════════════════════

# Define formula and estimate weighted logit
formula_baseline <- as.formula(paste("retrofit ~", paste(active_vars, collapse = " + ")))

print(summary(logit_baseline <- suppressWarnings(glm(formula  = formula_baseline,
                                                     family   = binomial(link = "logit"),
                                                     data     = data_clean,
                                                     weights  = weight_final)))
)

cat(sprintf("  Log-likelihood: %.2f | AIC: %.2f | BIC: %.2f\n",
            logLik(logit_baseline), AIC(logit_baseline), BIC(logit_baseline)))

# Extract robust standard errors (HC1)
vcov_robust <- vcovHC(logit_baseline, type = "HC1")
coef_robust <- coeftest(logit_baseline, vcov = vcov_robust)


# TABLE B.1: LOGIT MODEL ESTIMATION RESULTS (HC1 ROBUST)

# Variable labels of Appendix B (Tables B.1 and B.2)
appendix_labels <- c(
  "(Intercept)"      = "Intercept",
  "CDD_1723"         = "Cooling Degree Days (mean 2017-2023)",
  "HDD_1723"         = "Heating Degree Days (mean 2017-2023)",
  "ownership"        = "Ownership",
  "detached"         = "Detached dwelling",
  "moisturedamage"   = "Moisture / structural damage",
  "rooms"            = "Number of rooms",
  "pollution"        = "Polluted environment",
  "ln_income"        = "Income (log)",
  "allhigheducation" = "All adults with tertiary education",
  "meanageadults"    = "Mean age of adults (≥ 18 years)",
  "highlyurbanised"  = "Highly urbanised municipality",
  "kids"             = "Children present"
)

kblB1 <- data.frame(
  term             = rownames(coef_robust),
  estimate         = coef_robust[, 1],
  se_hc1           = coef_robust[, 2],
  z_stat           = coef_robust[, 3],
  p_value          = coef_robust[, 4],
  stringsAsFactors = FALSE,
  row.names        = NULL
) %>%
  mutate(
    Variable  = appendix_labels[term],
    Estimate  = sprintf("%.3f", estimate),
    `Rob. SE` = sprintf("%.3f", se_hc1),
    `z`       = sprintf("%.3f", z_stat),
    `p-value` = ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value)),
    ` `       = case_when(p_value < 0.001 ~ "***",
                          p_value < 0.01  ~ "**",
                          p_value < 0.05  ~ "*",
                          p_value < 0.10  ~ ".",
                          TRUE            ~ "")
  ) %>%
  mutate(term = factor(term, levels = c("(Intercept)", var_order))) %>%
  arrange(term) %>%
  select(Variable, Estimate, `Rob. SE`, `z`, `p-value`, ` `) %>%
  kbl(format    = "html",
      caption   = "Table B.1. Baseline binary logit: coefficient estimates with heteroskedasticity-robust standard errors (HC1)",
      align     = c("l", "r", "r", "r", "r", "c"),
      row.names = FALSE) %>%
  kable_classic(full_width = FALSE, html_font = "Times New Roman") %>%
  kable_styling(font_size = 10) %>%
  row_spec(0, bold = TRUE) %>%
  # Row positions follow var_order (defined in _setup.R): the intercept first,
  # then the twelve covariates in the six blocks of Section 3
  pack_rows("Intercept",                             1,  1,
            label_row_css = "background-color: #f0f0f0; font-weight: bold;") %>%
  pack_rows("A. Tenure",                             2,  2,
            label_row_css = "background-color: #f0f0f0; font-weight: bold;") %>%
  pack_rows("B. Income",                             3,  3,
            label_row_css = "background-color: #f0f0f0; font-weight: bold;") %>%
  pack_rows("C. Education",                          4,  4,
            label_row_css = "background-color: #f0f0f0; font-weight: bold;") %>%
  pack_rows("D. Climatic exposure",                  5,  6,
            label_row_css = "background-color: #f0f0f0; font-weight: bold;") %>%
  pack_rows("E. Dwelling characteristics",           7, 10,
            label_row_css = "background-color: #f0f0f0; font-weight: bold;") %>%
  pack_rows("F. Household composition and location", 11, 13,
            label_row_css = "background-color: #f0f0f0; font-weight: bold;") %>%
  footnote(
    general = sprintf(
      paste0(
        "N = %s households. Maximum likelihood estimation with survey weights trimmed at the 99.5th percentile. ",
        "Log-likelihood = %.2f; AIC = %.2f; BIC = %.2f. ",
        "*** p < 0.001, ** p < 0.01, * p < 0.05, . p < 0.10."
      ),
      format(nrow(data_clean), big.mark = ","),
      as.numeric(logLik(logit_baseline)),
      AIC(logit_baseline),
      BIC(logit_baseline)
    ),
    general_title     = "Notes:",
    footnote_as_chunk = TRUE
  )

print(kblB1)
save_table(kblB1, "Table_B1.html")

# Clean up (appendix_labels is kept for Table B.2)
remove(kblB1, vcov_robust, coef_robust)
gc()


# Baseline logit AMEs: bootstrap (1k iterations)

set.seed(123)
n_boot <- 1000

# Create environment for the bootstrap to loop inside it
.eb       <- new.env()
.eb$iter  <- 0
.eb$t0    <- Sys.time()

boot_ame_fn <- function(data, indices) {
  .eb$iter <- .eb$iter + 1
  
  # Add ETA messages each 100 iterations
  if (.eb$iter %% 100 == 0) {
    el  <- as.numeric(Sys.time() - .eb$t0, units = "secs")
    eta <- (el / .eb$iter) * (n_boot - .eb$iter)
    cat(sprintf("    [%d/%d] ETA: %02d:%02d\n",
                .eb$iter, n_boot, eta %/% 60, round(eta %% 60)))
  }
  
  d <- as.data.frame(data[indices, ])
  m <- suppressWarnings(glm(formula_baseline, family = binomial("logit"),
                            data = d, weights = weight_final))
  
  # Survey weights of this bootstrap sample. The GLM is fitted with weights, so
  # the AMEs must average over the same weighted population, not over the raw
  # sample. Matches the profile-specific AMEs in script 04.
  w_d <- d$weight_final
  
  ame <- setNames(numeric(length(active_vars)), active_vars)
  
  for (v in active_vars) {
    uv  <- unique(d[[v]][!is.na(d[[v]])])
    bin <- length(uv) == 2 && all(sort(uv) %in% c(0, 1))
    
    # For binary variables: discrete change from 0 to 1, refit
    if (bin) {
      d0 <- d; d0[[v]] <- 0L
      d1 <- d; d1[[v]] <- 1L
      ame[v] <- weighted.mean(predict(m, newdata = d1, type = "response") -
                                predict(m, newdata = d0, type = "response"),
                              w = w_d, na.rm = TRUE)
    } else {
      
      # For continuous variables, marginal effect scaled by one standard
      # deviation, weighted as in Table A.1
      pp     <- predict(m, newdata = d, type = "response")
      sd_x   <- sqrt(Hmisc::wtd.var(d[[v]], weights = w_d, normwt = TRUE, na.rm = TRUE))
      ame[v] <- coef(m)[v] * weighted.mean(pp * (1 - pp), w = w_d, na.rm = TRUE) * sd_x
    }
  }
  ame
}

# Run bootstrap
boot_ame <- boot(data      = data_clean, 
                 statistic = boot_ame_fn, 
                 R         = n_boot)

# Point estimate on the full sample. The bootstrap replications give the SE and
# the percentile intervals, but the reported estimate should be the one from the
# actual data, not the mean across resamples.
ame_point <- boot_ame_fn(data_clean, seq_len(nrow(data_clean)))

# Store AME results
ame_df <- data.frame(
  variable = active_vars,
  estimate = ame_point                                     * 100,
  se       = apply(boot_ame$t, 2, sd)                      * 100,
  ci_lo    = apply(boot_ame$t, 2, quantile, probs = 0.025) * 100,
  ci_hi    = apply(boot_ame$t, 2, quantile, probs = 0.975) * 100,
  row.names = NULL
) %>%
  mutate(sig = !(ci_lo < 0 & ci_hi > 0),
         label = var_labels[variable])

cat(sprintf("  AME estimation complete.\n"))

# TABLE B2: Average Marginal Effects (block-structured)
# block_map, block_order and var_order are defined in _setup.R, so that this
# table, Table B.1, Table 3 and Table A.1 all group the variables in the same
# six blocks used in Section 3 and in the Shapley attribution.

tblB2 <- ame_df %>%
  mutate(
    block    = block_map[variable],
    block    = factor(block,    levels = block_order),
    variable = factor(variable, levels = var_order)
  ) %>%
  arrange(block, variable) %>%
  transmute(
    block,
    Variable   = unname(appendix_labels[as.character(variable)]),
    `AME (pp)` = round(estimate, 2),
    `SE`       = round(se, 2),
    `95% CI`   = sprintf("[%.2f, %.2f]", ci_lo, ci_hi),
    `Sig.`     = ifelse(sig, "*", "")
  )

# Compute pack_rows indices from block structure
block_indices <- tblB2 %>%
  mutate(row = row_number()) %>%
  group_by(block) %>%
  summarise(start = min(row), end = max(row), .groups = "drop")

kblB2 <- tblB2 %>%
  select(-block) %>%
  kbl(format  = "html",
      caption = "Table B.2. Baseline binary logit: bootstrapped Average Marginal Effects",
      align   = c("l", "r", "r", "c", "c"),
      row.names = FALSE) %>%
  kable_classic(full_width = FALSE, html_font = "Times New Roman") %>%
  kable_styling(font_size = 11) %>%
  row_spec(0, bold = TRUE)

# Add pack_rows dynamically from block_indices
for (i in seq_len(nrow(block_indices))) {
  kblB2 <- kblB2 %>%
    pack_rows(
      as.character(block_indices$block[i]),
      block_indices$start[i],
      block_indices$end[i],
      label_row_css = "background-color: #f0f0f0; font-weight: bold;"
    )
}

kblB2 <- kblB2 %>%
  footnote(
    general = paste0(
      "Survey weights capped at the 99.5th percentile. Bootstrapped (1,000 iterations) 95% CIs. \n",
      "For continuous variables, the AME represents the effect of a one-standard-deviation increase. \n",
      "For binary variables, the AME represents the average difference in predicted probability associated with a discrete change from 0 to 1. \n",
      "* = CI excludes zero."),
    general_title     = "Notes: \n",
    footnote_as_chunk = TRUE
  )

print(kblB2)
save_table(kblB2, "Table_B2.html")

# FIGURE B.1: Forest plot of AMEs (block-structured)

# Short block labels for the facets, defined in _setup.R. The full block names
# used in the tables do not fit in the facet strips.

# Create dataframe for plotting
ame_plot_df <- ame_df %>%
  mutate(
    block = block_map_short[variable],
    block = factor(block, levels = block_order_short),
    label = unname(appendix_labels[variable])   # same labels as Table B.2
  ) %>%
  group_by(block) %>%
  mutate(label = fct_reorder(label, estimate)) %>%
  ungroup()

# Plot:
fig1 <- ggplot(ame_plot_df, aes(x = label, y = estimate)) +
  
  geom_hline(yintercept = 0, 
             linetype   = "solid",
             color      = "gray40", 
             linewidth  = 0.3) +
  
  geom_errorbar(aes(ymin = ci_lo, 
                    ymax = ci_hi),
                width     = 0, 
                linewidth = 0.5, 
                color     = "gray20") +
  
  geom_point(aes(fill = sig),
             shape = 21,
             size  = 2.5, 
             color = "gray20",
             stroke = 0.5) +
  
  scale_fill_manual(values = c("TRUE" = "gray20", "FALSE" = "white"),
                    guide  = "none") +
  
  facet_grid(block ~ ., scales = "free_y", space = "free_y") +
  
  coord_flip() +
  
  scale_y_continuous(
    breaks = seq(-4, 14, by = 2),
    labels = number_format(accuracy = 1),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  
  labs(
    title    = NULL,
    x        = NULL,
    y        = "Average marginal effect (pp)",
    caption  = paste0(
      "Notes: Filled points indicate a 95% bootstrapped CI (1,000 iterations) ",
      "that excludes zero."
    )
  ) +
  
  theme_minimal(base_size = 10, base_family = "Times New Roman") +
  theme(
    plot.caption          = element_text(size = 8, color = "gray30",
                                         hjust = 0, margin = margin(t = 8)),
    plot.caption.position = "plot",
    axis.text.y           = element_text(size = 9, color = "black"),
    axis.text.x           = element_text(size = 9, color = "black"),
    axis.title.x          = element_text(size = 10, margin = margin(t = 6)),
    panel.spacing.y       = unit(0.6, "lines"),
    panel.background      = element_blank(),
    panel.grid.major.x    = element_line(color = "gray85", linewidth = 0.3),
    panel.grid.minor.x    = element_blank(),
    panel.grid.major.y    = element_blank(),
    strip.text.y          = element_blank(),
    strip.background      = element_blank(),
    strip.placement       = "outside",
    axis.ticks.x          = element_line(color = "gray40", linewidth = 0.3),
    axis.ticks.length.x   = unit(2, "pt"),
    axis.ticks.y          = element_blank(),
    plot.margin           = margin(t = 5, r = 10, b = 5, l = 5)
  )

print(fig1)
save_figure(fig1, "Figure_B1",
            w_pdf = 6.5, h_pdf = 4, w_tiff = 20, h_tiff = 15)

# Clean up
# block_map, block_order and var_order live in _setup.R and are kept
remove(boot_ame, boot_ame_fn, .eb, tblB2, fig1, ame_plot_df,
       block_indices, kblB2, appendix_labels)
gc()



# ══════════════════════════════════════════════════════════════════════════════
# SAVE INTERMEDIATE OBJECTS FOR THE NEXT STAGE
# ══════════════════════════════════════════════════════════════════════════════
# NOTE: ame_df is kept for the next stages.

save(data_clean, logit_baseline, formula_baseline, ame_df,
     active_vars, var_labels, outcome_var,
     file = "02_descriptives_logit.RData")

cat(sprintf("\nStage 02 complete. Runtime: %.1f minutes\n",
            as.numeric(difftime(Sys.time(), start_time, units = "mins"))))
