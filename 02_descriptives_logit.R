###############################################################################
# Project       : [rehabSpain] [2] DESCRIPTIVES & BASELINE LOGIT
# Creation date : 02/12/2024
# Last update   : 22/06/2026
# Author        : Mercè Amich (merce.amich@ehu.eus)
# Institution   : UPV/EHU, BC3
# Last run time : 3.1 min.
# Script Overview:
#   Data preparation, outcome distribution and descriptive statistics,
#   baseline weighted binary logit, HC1-robust estimates and bootstrapped
#   average marginal effects (AMEs).
# Requirements:
#   - This script must be in the same directory as:
#       a) _setup.R
#       b) data.RData   (output of 01_build_data.R)
# Output  : Publication-ready tables (kableExtra) + ggplots (TIFF)
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
  "highpopulation",        # High dense populated area
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
  "allhigheducation" = "High educated",
  "meanageadults"    = "Mean age of adults",
  "highpopulation"   = "High population municipality",
  "kids"             = "Children present"
)

outcome_var   <- "retrofit" # Dummy (0, 1)
required_vars <- c(outcome_var, active_vars, "weight", "region") # Intersect

# Data cleaning: listwise deletion on required variables
n_total    <- nrow(data)
data_clean <- data[complete.cases(data[, required_vars]), ] # Keep if no NA's
n_removed  <- n_total - nrow(data_clean) # Keep track of removed observations

# Report of removed observations (both in this and in the first script)
cat(sprintf("  Removed %d obs (%.1f%%) due to missing values\n  Total removed in 1st script:   %.2f%% (survey), %.2f%% (weighted)\n  Total removed in both scripts: %.2f%% (survey), %.2f%% (weighted)\n  Final sample: n = %s\n",
            n_removed, 
            100 * n_removed / n_total,
            sample_info$pct_removed_cumulative_unweighted,
            sample_info$pct_removed_cumulative_weighted,
            100 * (sample_info$n_original_unweighted - nrow(data_clean)) / sample_info$n_original_unweighted,
            100 * (sample_info$n_original_weighted - sum(data_clean$weight)) / sample_info$n_original_weighted,
            format(nrow(data_clean), big.mark = ",")))

# Normalise and cap survey weights at 99.5th percentile
w_raw  <- data_clean$weight
w_norm <- w_raw / mean(w_raw, na.rm = TRUE)
w_cap  <- pmin(w_norm, quantile(w_norm, 0.995, na.rm = TRUE))
data_clean$weight_final <- w_cap / mean(w_cap)

# Percentage of adopters (retrofit == 1) in the finalsample
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
  select(Category, Outcome,
         n_unweighted, pct_unweighted, pct_weighted)

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
  caption   = "Table 1. Distribution of renovation measures and adoption status",
  align     = c("l", "c", "r", "r", "r"),
  row.names = FALSE,
  col.names = c("Category", "Outcome", 
                "n (unweighted)", "% (unweighted)", "% (weighted)")
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
      "N = %s households. Survey weights capped at 99.5th percentile.",
      format(nrow(data_clean), big.mark = ",")
    ),
    general_title     = "Notes:",
    footnote_as_chunk = TRUE
  )

print(kbl1)
save_table(kbl1, "Table1_Outcome_distribution.html")

# TABLE 2: DESCRIPTIVE STATISTICS

# Variable descriptions (expanded)
var_descriptions <- c(
  "CDD_1723"         = "Cooling Degree Days (mean value 2017-2023)",
  "HDD_1723"         = "Heating Degree Days (mean value 2017-2023)",
  "ownership"        = "Homeownership",
  "detached"         = "Detached dwelling",
  "moisturedamage"   = "Moisture / structurally damaged dwelling",
  "rooms"            = "Number of rooms",
  "pollution"        = "Polluted environment",
  "ln_income"        = "Income (log)",
  "allhigheducation" = "All members have high education",
  "meanageadults"    = "Mean age of adults (≥18y)",
  "highpopulation"   = "High density and populated municipality",
  "kids"             = "Children present in the household"
)

# Variable types
continuous_vars <- c("CDD_1723", "HDD_1723", "rooms", "ln_income", "meanageadults")
binary_vars     <- c("ownership", "detached", "moisturedamage", "pollution", 
                     "allhigheducation", "highpopulation", "kids")

# Outcome variable
outcome_stats <- data.frame(
  Variable    = "retrofit",
  Description = "Retrofit adoption",
  Type        = "Binary",
  Mean        = sprintf("%.2f", mean(data_clean$retrofit)),
  SD          = "-",
  Min         = "-",
  Median      = "-",
  Max         = "-",
  `Yes (%)`   = sprintf("%.1f%%", 100 * mean(data_clean$retrofit)),
  `No (%)`    = sprintf("%.1f%%", 100 * (1 - mean(data_clean$retrofit))),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Function to compute stats for continuous variables
compute_continuous <- function(var_name) {
  x <- data_clean[[var_name]]
  data.frame(
    Variable    = var_name,
    Description = var_descriptions[var_name],
    Type        = "Continuous",
    Mean        = sprintf("%.2f", mean(x,   na.rm = TRUE)),
    SD          = sprintf("%.2f", sd(x,     na.rm = TRUE)),
    Min         = sprintf("%.2f", min(x,    na.rm = TRUE)),
    Median      = sprintf("%.2f", median(x, na.rm = TRUE)),
    Max         = sprintf("%.2f", max(x,    na.rm = TRUE)),
    `Yes (%)`   = "-",
    `No (%)`    = "-",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

# Function to compute stats for binary variables
compute_binary <- function(var_name) {
  x <- data_clean[[var_name]]
  pct_yes <- mean(x, na.rm = TRUE)
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

# Compute statistics for each variable group:

# A. DEPENDENT VARIABLE
desc_outcome <- outcome_stats

# B. CLIMATIC VARIABLES
climatic_vars <- c("CDD_1723", "HDD_1723")
desc_climatic <- do.call(rbind, lapply(climatic_vars, compute_continuous))

# C. BUILDING VARIABLES
building_vars <- c("ownership", "detached", "moisturedamage", "rooms", "pollution")
desc_building <- rbind(
  compute_binary("ownership"),
  compute_binary("detached"),
  compute_binary("moisturedamage"),
  compute_continuous("rooms"),
  compute_binary("pollution")
)

# D. SOCIO-DEMOGRAPHIC VARIABLES
sociodem_vars <- c("ln_income", "allhigheducation", "meanageadults", 
                   "highpopulation", "kids")
desc_sociodem <- rbind(
  compute_continuous("ln_income"),
  compute_binary("allhigheducation"),
  compute_continuous("meanageadults"),
  compute_binary("highpopulation"),
  compute_binary("kids")
)

# Combine all sections
desc_stats_full <- rbind(
  desc_outcome,
  desc_climatic,
  desc_building,
  desc_sociodem
)

# Row indices for section headers
n_outcome  <- 1
n_climatic <- 2
n_building <- 7
n_sociodem <- 12

# Create kableExtra table with section headers
kblS1 <- kbl(
  desc_stats_full %>% select(-Variable),
  format    = "html",
  caption   = "Table 2. Descriptive statistics",
  align     = c("l", "c", "r", "r", "r", "r", "r", "r", "r", "r"),
  row.names = FALSE,
  digits    = 2
) %>%
  kable_classic(full_width = FALSE, html_font = "Times New Roman") %>%
  kable_styling(font_size = 10) %>%
  row_spec(0, bold = TRUE) %>%
  
  # Section A: Outcome variable
  pack_rows("A. Dependent variable", 1, n_outcome, 
            label_row_css = "background-color: #f0f0f0; font-weight: bold;") %>%
  
  # Section B: Climatic variables
  pack_rows("B. Climatic variables", n_outcome + 1, n_outcome + n_climatic, 
            label_row_css = "background-color: #f0f0f0; font-weight: bold;") %>%
  
  # Section C: Building characteristics
  pack_rows("C. Dwelling and environmental characteristics", n_outcome + n_climatic + 1, 
            n_outcome + n_climatic + n_building - n_climatic, 
            label_row_css = "background-color: #f0f0f0; font-weight: bold;") %>%
  
  # Section D: Socio-demographic characteristics
  pack_rows("D. Socio-demographic characteristics", 
            n_outcome + n_climatic + n_building - n_climatic + 1, 
            nrow(desc_stats_full), 
            label_row_css = "background-color: #f0f0f0; font-weight: bold;") %>%
  
  footnote(
    general = sprintf("N = %s households. Binary variables coded as 1 = Yes, 0 = No. Survey weights capped at 99.5th percentile.",
                      format(nrow(data_clean), big.mark = ",")),
    general_title = "Notes:",
    footnote_as_chunk = TRUE
  )

print(kblS1)
save_table(kblS1, "TableS1_Descriptive_statistics.html")

# Clean up
remove(desc_outcome, desc_climatic, desc_building, desc_sociodem, desc_stats_full,
       climatic_vars, building_vars, sociodem_vars, continuous_vars, binary_vars,
       var_descriptions, outcome_stats, n_outcome, n_climatic, n_building, n_sociodem,
       compute_continuous, compute_binary, kbl1, kblS1, sample_info, adopter_detail,
       adopter_summary, non_adopter_summary, outcome_dist_raw)
gc()



# ══════════════════════════════════════════════════════════════════════════════
# 2. BASELINE BINARY LOGIT ADOPTION MODEL
# ══════════════════════════════════════════════════════════════════════════════

# Define formula and estimate weighted logit
formula_baseline <- as.formula(paste("retrofit ~", paste(active_vars, collapse = " + ")))

print(summary(logit_baseline <- glm(formula  = formula_baseline,
                              family   = binomial(link = "logit"),
                              data     = data_clean,
                              weights  = weight_final))
)

cat(sprintf("  Log-likelihood: %.2f | AIC: %.2f | BIC: %.2f\n",
            logLik(logit_baseline), AIC(logit_baseline), BIC(logit_baseline)))

# Extract robust standard errors (HC1)
vcov_robust <- vcovHC(logit_baseline, type = "HC1")
coef_robust <- coeftest(logit_baseline, vcov = vcov_robust)


# TABLE S2: LOGIT MODEL ESTIMATION RESULTS (HC1 ROBUST)

appendix_labels <- c(
  "(Intercept)"      = "Intercept",
  "CDD_1723"         = "Cooling Degree Days (mean 2017-2023)",
  "HDD_1723"         = "Heating Degree Days (mean 2017-2023)",
  "ownership"        = "Homeownership",
  "detached"         = "Detached dwelling",
  "moisturedamage"   = "Moisture / structural damage",
  "rooms"            = "Number of rooms",
  "pollution"        = "Self-reported polluted environment",
  "ln_income"        = "Household income (log)",
  "allhigheducation" = "All adults with tertiary education",
  "meanageadults"    = "Mean age of adults (>=18 years)",
  "highpopulation"   = "High-population municipality (>10,000 inhab.)",
  "kids"             = "Children present in household"
)

kblS2 <- data.frame(
  term    = rownames(coef_robust),
  estimate = coef_robust[, 1],
  se_hc1   = coef_robust[, 2],
  z_stat   = coef_robust[, 3],
  p_value  = coef_robust[, 4],
  stringsAsFactors = FALSE,
  row.names = NULL
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
                          p_value < 0.10  ~ "†",
                          TRUE            ~ "")
  ) %>%
  mutate(term = factor(term, levels = c("(Intercept)", active_vars))) %>%
  arrange(term) %>%
  select(Variable, Estimate, `Rob. SE`, `z`, `p-value`, ` `) %>%
  kbl(format    = "html",
      caption   = "Table S2. Baseline binary logit coefficient estimates with heteroskedasticity-robust standard errors (HC1)",
      align     = c("l", "r", "r", "r", "r", "c"),
      row.names = FALSE) %>%
  kable_classic(full_width = FALSE, html_font = "Times New Roman") %>%
  kable_styling(font_size = 10) %>%
  row_spec(0, bold = TRUE) %>%
  pack_rows("Intercept",                                     1,  1,
            label_row_css = "background-color: #f0f0f0; font-weight: bold;") %>%
  pack_rows("A. Climatic variables",                         2,  3,
            label_row_css = "background-color: #f0f0f0; font-weight: bold;") %>%
  pack_rows("B. Dwelling and environmental characteristics",  4,  8,
            label_row_css = "background-color: #f0f0f0; font-weight: bold;") %>%
  pack_rows("C. Socio-demographic characteristics",          9, 13,
            label_row_css = "background-color: #f0f0f0; font-weight: bold;") %>%
  footnote(
    general = sprintf(
      paste0(
        "N = %s households. Maximum likelihood estimation with survey weights capped at 99.5th percentile. ",
        "Log-likelihood = %.2f; AIC = %.2f; BIC = %.2f. ",
        "*** p < 0.001, ** p < 0.01, * p < 0.05, † p < 0.10."
      ),
      format(nrow(data_clean), big.mark = ","),
      as.numeric(logLik(logit_baseline)),
      AIC(logit_baseline),
      BIC(logit_baseline)
    ),
    general_title     = "Notes:",
    footnote_as_chunk = TRUE
  )

print(kblS2)
save_table(kblS2, "TableS2_Logit_robust.html")

# Clean up
remove(appendix_labels, kblS2, vcov_robust, coef_robust)
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
  m <- glm(formula_baseline, family = binomial("logit"), 
           data = d, weights = weight_final)
  
  ame <- setNames(numeric(length(active_vars)), active_vars)
  
  for (v in active_vars) {
    uv  <- unique(d[[v]][!is.na(d[[v]])])
    bin <- length(uv) == 2 && all(sort(uv) %in% c(0, 1))
    
    # For binary variables: discrete change from 0 to 1, refit
    if (bin) {
      d0 <- d; d0[[v]] <- 0L
      d1 <- d; d1[[v]] <- 1L
      ame[v] <- mean(predict(m, newdata = d1, type = "response") -
                     predict(m, newdata = d0, type = "response"), na.rm = TRUE)
    } else {
      
      # For continuous variables, add 1 standard deviation, refit
      pp     <- predict(m, newdata = d, type = "response")
      sd_x   <- sd(d[[v]], na.rm = TRUE)
      ame[v] <- coef(m)[v] * mean(pp * (1 - pp), na.rm = TRUE) * sd_x
      }
  }
  ame
}

boot_ame <- boot(data      = data_clean, 
                 statistic = boot_ame_fn, 
                 R         = n_boot)

# Store AME results
ame_df <- data.frame(
  variable = active_vars,
  estimate = apply(boot_ame$t, 2, mean)                    * 100,
  se       = apply(boot_ame$t, 2, sd)                      * 100, # is this ok?????
  ci_lo    = apply(boot_ame$t, 2, quantile, probs = 0.025) * 100,
  ci_hi    = apply(boot_ame$t, 2, quantile, probs = 0.975) * 100,
  row.names = NULL
) %>%
  mutate(sig = !(ci_lo < 0 & ci_hi > 0),
         label = var_labels[variable])

cat(sprintf("  AME estimation complete.\n"))

# TABLE S3: Average Marginal Effects (block-structured)
# Assign thematic blocks matching Data section order
block_map <- c(
  "ln_income"        = "C. Socio-demographic characteristics",
  "allhigheducation" = "C. Socio-demographic characteristics",
  "meanageadults"    = "C. Socio-demographic characteristics",
  "highpopulation"   = "C. Socio-demographic characteristics",
  "kids"             = "C. Socio-demographic characteristics",
  "ownership"        = "B. Dwelling and environmental characteristics",
  "detached"         = "B. Dwelling and environmental characteristics",
  "moisturedamage"   = "B. Dwelling and environmental characteristics",
  "rooms"            = "B. Dwelling and environmental characteristics",
  "pollution"        = "B. Dwelling and environmental characteristics",
  "CDD_1723"         = "A. Climatic variables",
  "HDD_1723"         = "A. Climatic variables"
)

block_order <- c(
  "A. Climatic variables",
  "B. Dwelling and environmental characteristics",
  "C. Socio-demographic characteristics"
)

tblS3 <- ame_df %>%
  mutate(
    block = block_map[variable],
    block = factor(block, levels = block_order),
    # Within each block, sort by descending absolute AME:
    abs_est = abs(estimate)
  ) %>%
  arrange(block, desc(abs_est)) %>%
  transmute(
    block,
    Variable   = label,
    `AME (%)` = round(estimate, 2),
    `SE`       = round(se, 2),
    `95% CI`   = sprintf("[%.2f, %.2f]", ci_lo, ci_hi),
    `Sig.`     = ifelse(sig, "*", "")
  )
# Compute pack_rows indices from block structure
block_indices <- tblS3 %>%
  mutate(row = row_number()) %>%
  group_by(block) %>%
  summarise(start = min(row), end = max(row), .groups = "drop")
kbls3 <- tblS3 %>%
  select(-block) %>%
  kbl(format  = "html",
      caption = "Table S3. Bootstrapped average marginal effects",
      align   = c("l", "r", "r", "c", "c"),
      row.names = FALSE) %>%
  kable_classic(full_width = FALSE, html_font = "Times New Roman") %>%
  kable_styling(font_size = 11) %>%
  row_spec(0, bold = TRUE)
# Add pack_rows dynamically from block_indices
for (i in seq_len(nrow(block_indices))) {
  kbls3 <- kbls3 %>%
    pack_rows(
      as.character(block_indices$block[i]),
      block_indices$start[i],
      block_indices$end[i],
      label_row_css = "background-color: #f0f0f0; font-weight: bold;"
    )
}
kbls3 <- kbls3 %>%
  footnote(
    general = paste0(
      "AME in percentage points (pp). Continuous variables scaled by one standard deviation. ",
      "95% CI from 1,000 bootstrap iterations. * = CI excludes zero."
    ),
    general_title     = "Notes:",
    footnote_as_chunk = TRUE
  )

print(kbls3)
save_table(kbls3, "TableS3_AME.html")

# FIGURE 1: Forest plot of AMEs (block-structured)

# Add block assignment to ame_df
block_map <- c(
  "ln_income"        = "Socio-demographic",
  "allhigheducation" = "Socio-demographic",
  "meanageadults"    = "Socio-demographic",
  "highpopulation"   = "Socio-demographic",
  "kids"             = "Socio-demographic",
  "ownership"        = "Dwelling and environmental",
  "detached"         = "Dwelling and environmental",
  "moisturedamage"   = "Dwelling and environmental",
  "rooms"            = "Dwelling and environmental",
  "pollution"        = "Dwelling and environmental",
  "CDD_1723"         = "Climatic",
  "HDD_1723"         = "Climatic"
)

block_order <- c(
  "Climatic",
  "Dwelling and environmental",
  "Socio-demographic"
)

# Create dataframe for plotting
ame_plot_df <- ame_df %>%
  mutate(
    block = block_map[variable],
    block = factor(block, levels = block_order),
    label = var_labels[variable]
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
      "Notes: Filled points indicate 95% CI excludes zero. ",
      "Continuous variables: effect of one standard deviation change. ",
      "Bootstrapped CIs (1,000 iterations)."
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
    strip.text.y          = element_text(face = "bold", size = 8, angle = 270,
                                         hjust = 0.5, vjust = 0.5, color = "black"),
    strip.background      = element_blank(),
    strip.placement       = "outside",
    axis.ticks.x          = element_line(color = "gray40", linewidth = 0.3),
    axis.ticks.length.x   = unit(2, "pt"),
    axis.ticks.y          = element_blank(),
    plot.margin           = margin(t = 5, r = 10, b = 5, l = 5)
  )

print(fig1)
ggsave("Figure1_AME_forest.pdf", fig1, width = 6.5, height = 4.5, device = cairo_pdf)

ggsave("Figure1_AME_forest.tiff",
       plot   = fig1,
       device = "tiff",
       width  = 20,
       height = 17,
       units  = "cm",
       dpi    = 900)

# Clean up
remove(boot_ame, boot_ame_fn, .eb, tblS3, fig1, 
       ame_plot_df, block_map, block_order, 
       block_indices, kbls3)
gc()



# ══════════════════════════════════════════════════════════════════════════════
# SAVE INTERMEDIATE OBJECTS FOR THE NEXT STAGE
# ══════════════════════════════════════════════════════════════════════════════
# NOTE: ame_df kept because the
# profile-specific Figure 3 (PDF) in stage 04 references it.

save(data_clean, logit_baseline, formula_baseline, ame_df,
     active_vars, var_labels, outcome_var,
     file = "02_descriptives_logit.RData")

cat(sprintf("\nStage 02 complete. Runtime: %.1f minutes\n",
            as.numeric(difftime(Sys.time(), start_time, units = "mins"))))
