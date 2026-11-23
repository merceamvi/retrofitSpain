###############################################################################
# Project       : [rehabSpain] [5] DECOMPOSITION OF THE ADOPTION GAP
# Creation date : 10/08/2026
# Last update   : 23/11/2026 (profiles numbered by adoption in stage 03, so the
#                 reference is P1; Table 4 note on unsupported interactions built
#                 from the data; robustness table renamed D.3, as in the paper;
#                 object names, console output and CSV columns in English;
#                 renumbered as stage 05)
# Author        : Mercè Amich (merce.amich@ehu.eus)
# Institution   : UPV/EHU, BC3
# Last run time : ~16 min with N_BOOT_DECOMP = 500; set it to 0 to skip the CIs

# Script Overview:
#   Profiles differ sharply in adoption. That difference can come from two
#   sources: households in each profile respond differently to the same
#   characteristics (different betas), or they simply hold different
#   characteristics (different Xs). This script separates the two.
#     1. Interacted model with the support rule: are the betas different?
#        (joint and per-variable Wald, Table 4)
#     2. Common-slope model (Model A): profile fixed effects, shared beta.
#     3. Aggregate decomposition of each gap vs the reference (Table 5).
#     4. Shapley attribution of the composition component by block, with
#        Fairlie (2005) rank matching (Table 6).
#     5. Model C: separate models by profile, pooled baseline-logit reference;
#        5b. Model B: profile-specific age coefficient. Models A, B and C are
#        compared in Table D.3.
#     6. Stratified bootstrap for the confidence intervals.
#     7-9. Main results table, manuscript tables.

# Requirements:
#   - This script must be in the same directory as:
#       a) _setup.R
#       b) 03_famd_clustering.RData   (output of 03_famd_clustering.R)

# Output  : paper/  Table_4, Table_5, Table_6, Table_D3 (HTML)
#           output/ decomposition tables behind the text (CSV)
#           + 05_decomposition.RData

###############################################################################

# ══════════════════════════════════════════════════════════════════════════════
# 0. SETUP
# ══════════════════════════════════════════════════════════════════════════════

rm(list = ls(all = TRUE))
Sys.setenv(LANG = "en")

source("_setup.R")
load("03_famd_clustering.RData")

rdata_dir <- "."

start_time <- Sys.time(); set.seed(123)

# Bootstrap replications. Read from the environment so it survives the rm()
n_boot <- as.integer(Sys.getenv("N_BOOT_DECOMP", "500"))

# Matching draws in the point estimate
n_match <- as.integer(Sys.getenv("N_MATCH", "50"))

# Matching draws within each bootstrap replication
n_match_b <- 5

# Reference profile: the one with the highest weighted adoption, which is P1
# because stage 03 numbers the profiles by adoption. Override by setting the
# REF_PROFILE environment variable; never fixed by hand.
adopt_w <- sapply(sort(unique(as.character(data_clean$profile_gmm))), function(g) {
  i <- as.character(data_clean$profile_gmm) == g
  weighted.mean(data_clean[[outcome_var]][i],
                w = data_clean$weight_final[i], na.rm = TRUE)
})
ref_profile <- Sys.getenv("REF_PROFILE", names(which.max(adopt_w)))
cat(sprintf("  Reference profile (highest weighted adoption): P%s\n", ref_profile))

# The reference level's delta is zero in the common-slope model.
data_clean$profile_gmm <- factor(data_clean$profile_gmm)
data_clean$profile_gmm <- relevel(data_clean$profile_gmm, ref = ref_profile)
profiles <- levels(data_clean$profile_gmm)

# Profiles compared against the reference
comparisons <- setdiff(profiles, ref_profile)

# Weights, profile (as text) and binary outcome
w <- data_clean$weight_final
grp <- as.character(data_clean$profile_gmm)
y <- data_clean[[outcome_var]]

cat(
  "\n=== [5] DECOMPOSITION OF THE ADOPTION GAP ===\n",
  sprintf("  N = %s households. Reference profile: P%s\n",
          format(nrow(data_clean), big.mark = ","), ref_profile),
  sprintf("  Bootstrap: %d replications\n", n_boot)
)

# Variable blocks for the Shapley attribution. HDD and CDD are SEPARATE blocks
# (21/08/2026) so the climatic contribution reads as heating vs cooling;
# K = 7 blocks -> 2^7 = 128 subsets per comparison, still tractable.
blocks <- list(
  tenure      = c("ownership"),
  income      = c("ln_income"),
  education   = c("allhigheducation"),
  climate_hdd = c("HDD_1723"),
  climate_cdd = c("CDD_1723"),
  dwelling    = c("rooms", "detached", "moisturedamage", "pollution"),
  household   = c("meanageadults", "kids", "highlyurbanised")
)

# The partition must be exhaustive and non-overlapping
stopifnot(setequal(unlist(blocks), active_vars))

block_names <- names(blocks)
K <- length(block_names)

cat(sprintf("  Shapley over %d blocks (%d subsets per comparison)\n", K, 2^K))


# ══════════════════════════════════════════════════════════════════════════════
# 1. INTERACTED MODEL: ARE THE BETAS DIFFERENT?
# ══════════════════════════════════════════════════════════════════════════════
# Weighted logit with profile fixed effects and profile x characteristic
# interactions. SUPPORT RULE (Section 4.3): the interaction of a BINARY
# characteristic with a profile is estimated only if, within that profile, the
# less frequent category contains at least MIN_EVENTS adopters and MIN_EVENTS
# non-adopters; otherwise the profile keeps the reference profile's slope, since
# the interaction would be identified by a handful of observations. Continuous
# characteristics are always interacted. Inference uses HC1.

MIN_EVENTS <- 10

vars_bin <- active_vars[vapply(active_vars, function(v)
  all(na.omit(data_clean[[v]]) %in% 0:1), logical(1))]

# Support of every profile x binary characteristic cell
support <- do.call(rbind, lapply(profiles, function(g) {
  do.call(rbind, lapply(vars_bin, function(v) {
    m  <- grp == g
    x  <- data_clean[[v]][m]
    yy <- y[m]
    tb <- table(factor(x, levels = 0:1))
    mv <- as.numeric(names(tb)[which.min(tb)])   # less frequent category
    k  <- x == mv
    data.frame(profile      = g,
               variable     = v,
               min_category = mv,
               n_minority   = sum(k),
               adopters     = sum(yy[k] == 1),
               non_adopters = sum(yy[k] == 0),
               supported    = sum(yy[k] == 1) >= MIN_EVENTS &
                 sum(yy[k] == 0) >= MIN_EVENTS,
               stringsAsFactors = FALSE)
  }))
}))

# Interaction columns built explicitly, so that unsupported cells are left out
d_int     <- data_clean
int_terms <- data.frame()
excluded  <- data.frame()
for (g in comparisons) {
  d_int[[paste0("P", g)]] <- as.numeric(grp == g)
  for (v in active_vars) {
    ok <- TRUE
    if (v %in% vars_bin) {
      ok <- support$supported[support$profile == g & support$variable == v] &&
        support$supported[support$profile == ref_profile & support$variable == v]
    }
    if (ok) {
      nm <- paste0("int_P", g, "_", v)
      d_int[[nm]] <- d_int[[paste0("P", g)]] * d_int[[v]]
      int_terms <- rbind(int_terms, data.frame(term = nm, profile = g, variable = v,
                                               stringsAsFactors = FALSE))
    } else {
      excluded <- rbind(excluded, data.frame(profile = g, variable = v,
                                             stringsAsFactors = FALSE))
    }
  }
}

formula_int <- as.formula(paste(
  outcome_var, "~",
  paste(c(paste0("P", comparisons), active_vars, int_terms$term), collapse = " + ")))

# Terms: (a) profile intercept shifts relative to the reference, (b) the 12
# characteristic main effects (reference-profile slopes), (c) the supported
# profile x characteristic interactions. The tests use only (c): common SLOPES
# are tested while LEVELS may differ.
m_int <- suppressWarnings(
  glm(formula_int,
      family  = binomial("logit"),
      data    = d_int,
      weights = weight_final))

# HC1 robust variance and coefficient vector
V_int <- sandwich::vcovHC(m_int, type = "HC1")
b_int <- coef(m_int)

cat(sprintf("\n-- 1. Support rule: %d interactions estimated, %d left out --\n",
            nrow(int_terms), nrow(excluded)))
if (nrow(excluded) > 0) print(support[!support$supported, ], row.names = FALSE)
cat(sprintf("   Largest hat value: %.4f\n", max(hatvalues(m_int))))

# ── 1a. Joint Wald over ALL estimated interactions ───────────────────────────
# H0: all interactions are zero (common slopes). Ha: at least one is non-zero.

idx_int <- intersect(int_terms$term,
                     intersect(names(b_int)[!is.na(b_int)],
                               rownames(V_int)))
bb <- b_int[idx_int]
VV <- V_int[idx_int, idx_int, drop = FALSE]

# Wald statistic b' V^-1 b, its df, and the chi-squared p-value
W_global  <- as.numeric(t(bb) %*% solve(VV) %*% bb)
df_global <- length(idx_int)
p_global  <- pchisq(W_global, df = df_global, lower.tail = FALSE)

cat("\n-- 1a. Joint Wald over all estimated interactions (HC1) --\n",
    sprintf("    W = %.2f, df = %d, p = %.4f\n", W_global, df_global, p_global))

# ── 1b. Per-variable Wald ────────────────────────────────────────────────────
# For each characteristic, test whether its estimated interactions are jointly
# zero (three, or fewer where the support rule leaves a profile out).

tbl_wald <- data.frame()

for (v in active_vars) {
  
  idx_v <- intersect(int_terms$term[int_terms$variable == v], idx_int)
  if (length(idx_v) == 0) next
  
  bb_v <- b_int[idx_v]
  VV_v <- V_int[idx_v, idx_v, drop = FALSE]
  
  W_v  <- as.numeric(t(bb_v) %*% solve(VV_v) %*% bb_v)
  df_v <- length(idx_v)
  p_v  <- pchisq(W_v, df = df_v, lower.tail = FALSE)
  
  tbl_wald <- rbind(tbl_wald, data.frame(
    variable = v,
    label    = unname(var_labels[v]),
    W        = round(W_v, 3),
    df       = df_v,
    p_value  = round(p_v, 4),
    p_raw    = p_v,          # unrounded, for the Bonferroni column in Table 4
    stringsAsFactors = FALSE))
}

rownames(tbl_wald) <- NULL
tbl_wald <- tbl_wald[order(tbl_wald$p_value), ]   # smallest p-value first

cat(
  paste0(
    "\n-- 1b. Per-variable Wald over its estimated interactions (HC1) --\n",
    paste(capture.output(print(tbl_wald)), collapse = "\n"),
    "\n"
  )
)

# ── 1c. Collinearity check ───────────────────────────────────────────────────
# Interacting 12 covariates with a four-level factor inflates the design. A high
# VIF does not invalidate the joint test; it is reported here as a console check.

X_int     <- model.matrix(m_int)
cols_vif  <- setdiff(colnames(X_int), "(Intercept)")
vif_vals  <- rep(NA_real_, length(cols_vif))
names(vif_vals) <- cols_vif

# Column-by-column VIF from auxiliary regressions
for (j in seq_along(cols_vif)) {
  y_j <- X_int[, cols_vif[j]]
  X_j <- X_int[, setdiff(cols_vif, cols_vif[j]), drop = FALSE]
  r2  <- summary(lm(y_j ~ X_j))$r.squared
  vif_vals[j] <- 1 / (1 - min(r2, 0.999999))   # cap R2 to avoid dividing by zero
}

tbl_vif <- data.frame(term = names(vif_vals),
                      VIF  = round(unname(vif_vals), 2),
                      stringsAsFactors = FALSE)
tbl_vif <- tbl_vif[order(-tbl_vif$VIF), ]
rownames(tbl_vif) <- NULL

cat("\n-- 1c. Collinearity of the interacted design (15 largest VIF) --\n")
print(head(tbl_vif, 15))
cat(sprintf("    Terms with VIF > 10: %d of %d\n",
            sum(tbl_vif$VIF > 10), nrow(tbl_vif)))

# Collinearity here is mechanical (a profile dummy and its interactions share
# most of their variation), not pathological.


# ══════════════════════════════════════════════════════════════════════════════
# 2. COMMON-SLOPE MODEL
# ══════════════════════════════════════════════════════════════════════════════
#   Pr(Y_i = 1 | X_i, g_i) = Lambda(alpha + delta_g + X_i'beta), with delta_ref = 0
#   beta    : common association structure
#   delta_g : level differences not explained by X
# The restriction is licensed by Section 1: the slope differences are not
# detectable, jointly or one variable at a time. "cannot be rejected", not
# "the slopes are equal".

formula_common <- update(formula_baseline, . ~ . + profile_gmm)

m_common <- suppressWarnings(
  glm(formula_common,
      family  = binomial("logit"),
      data    = data_clean,
      weights = weight_final))

V_common <- sandwich::vcovHC(m_common, type = "HC1")

cat("\n-- 2. Common-slope model with HC1 errors --\n")
print(lmtest::coeftest(m_common, vcov. = V_common))

# Profile fixed effects, extracted for readability
ct_common <- lmtest::coeftest(m_common, vcov. = V_common)
print(ct_common[grep("^profile_gmm", rownames(ct_common)), ])

# If no profile fixed effect is distinguishable from zero once characteristics
# are controlled for, the gaps are compositional. Report marginal deltas as
# "not distinguishable from zero", not "zero".


# ══════════════════════════════════════════════════════════════════════════════
# 3. AGGREGATE DECOMPOSITION VS THE REFERENCE PROFILE
# ══════════════════════════════════════════════════════════════════════════════
#   mu_ref        : mean predicted adoption in the reference profile
#   mu_g_own      : mean predicted adoption in g, with ITS intercept
#   mu_g_refalpha : mean predicted adoption in g, with its X but the reference intercept
#
#   Delta(ref,g) = [mu_ref - mu_g_refalpha] + [mu_g_refalpha - mu_g_own]
#                      composition             profile residual

coef_common <- coef(m_common)
alpha <- coef_common["(Intercept)"]   # intercept = the reference profile
beta  <- coef_common[active_vars]     # common slopes, in active_vars order
Xmat <- as.matrix(data_clean[, active_vars])

# Linear index WITHOUT the profile shift (everyone on the reference intercept)
eta_refalpha <- as.numeric(alpha + Xmat %*% beta)

# Profile shifts (reference stays at zero)
delta_profile <- setNames(numeric(length(profiles)), profiles)
for (g in profiles) {
  coef_name <- paste0("profile_gmm", g)
  if (coef_name %in% names(coef_common)) delta_profile[g] <- coef_common[coef_name]
}

# Linear index with each household's own intercept
eta_own <- eta_refalpha + delta_profile[grp]

# Mean predicted adoption of the reference
mu_ref <- weighted.mean(plogis(eta_refalpha[grp == ref_profile]),
                        w = w[grp == ref_profile])

tbl_decomp <- data.frame()

for (g in comparisons) {
  
  i <- grp == g
  
  # Predicted adoption with the own vs the reference intercept
  mu_g_own      <- weighted.mean(plogis(eta_own[i]),      w = w[i])
  mu_g_refalpha <- weighted.mean(plogis(eta_refalpha[i]), w = w[i])
  
  gap_pred    <- mu_ref - mu_g_own          # total predicted gap
  composition <- mu_ref - mu_g_refalpha     # attributable to characteristics
  residual    <- mu_g_refalpha - mu_g_own   # left unexplained
  
  # Observed gap, to check the model reproduces the profile means
  gap_obs <- weighted.mean(y[grp == ref_profile], w = w[grp == ref_profile]) -
    weighted.mean(y[i], w = w[i])
  
  tbl_decomp <- rbind(tbl_decomp, data.frame(
    comparison      = paste0("P", ref_profile, " vs P", g),
    observed_gap    = round(gap_obs     * 100, 2),
    predicted_gap   = round(gap_pred    * 100, 2),
    composition_pp  = round(composition * 100, 2),
    residual_pp     = round(residual    * 100, 2),
    # % share suppressed when the predicted gap is under 2 pp: a ratio against
    # a near-zero denominator is not interpretable
    composition_pct  = ifelse(abs(gap_pred) < 0.02, NA_real_,
                              round(100 * composition / gap_pred, 1)),
    stringsAsFactors = FALSE))
}

cat("\n-- 3. Aggregate decomposition (pp). gap = composition + residual --\n")
print(tbl_decomp)
cat("   NOTE: for a very small gap the split is unstable; check in section 7\n")
cat("   whether both intervals cross zero, and do not interpret that split.\n")


# ══════════════════════════════════════════════════════════════════════════════
# 4. SHAPLEY ATTRIBUTION OF THE COMPOSITION COMPONENT
# ══════════════════════════════════════════════════════════════════════════════
# Substituting variables one at a time is order-dependent. Shapley averages each
# block's marginal contribution over ALL orders, so the contributions are
# order-invariant and sum exactly to the composition C.
#
# v(S) = mean predicted adoption of profile g when the blocks in S take the
# reference profile's values and the rest keep g's, holding beta fixed and
# evaluating every household at the REFERENCE intercept (alpha, delta_g = 0).
# Hence v(0) = mu_g|ref and v(K) = mu_ref, and the contributions sum to the
# composition C of Section 3. The counterfactual is built by matching, following
# Fairlie (2005): samples of equal size are drawn from profile g and from the
# reference profile (size = the smaller profile), each household is ranked by its
# predicted adoption under the common coefficients, and households are paired by
# rank. Fairlie's sequential substitution is order-dependent; here every order is
# averaged through Shapley weights. With survey weights, BOTH samples are drawn
# with replacement with probability proportional to the weight, so each sample
# reproduces its profile's weighted distribution and the means over matched pairs
# are unweighted. Results are averaged over n_match draws.

# All 2^K subsets of blocks, as a TRUE/FALSE matrix
subsets <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), K)))
colnames(subsets) <- block_names

# Shapley weight of each subset: |S|! (K-|S|-1)! / K!
shapley_weights <- numeric(nrow(subsets))
for (r in seq_len(nrow(subsets))) {
  size_s <- sum(subsets[r, ])
  # The full set is never a base for adding a block; leave its weight at zero
  if (size_s == K) { shapley_weights[r] <- 0; next }
  shapley_weights[r] <- factorial(size_s) * factorial(K - size_s - 1) / factorial(K)
}

# For each subset and absent block, the row holding that subset plus the block.
# Precomputed so the bootstrap loops do not repeat the search.
row_with_block <- matrix(NA_integer_, nrow = nrow(subsets), ncol = K)
for (r in seq_len(nrow(subsets))) {
  for (k in seq_len(K)) {
    if (subsets[r, k]) next
    target <- subsets[r, ]
    target[k] <- TRUE
    for (r2 in seq_len(nrow(subsets))) {
      if (all(subsets[r2, ] == target)) { row_with_block[r, k] <- r2; break }
    }
  }
}

# ── Fairlie (2005) rank matching with Shapley attribution over blocks ───────
# X_g, w_g     : characteristics and weights of profile g
# X_ref, w_ref : characteristics and weights of the reference profile
# alpha, beta  : reference intercept and common coefficients
# n_draws      : number of matched samples averaged
# Returns the block contributions (phi) and v(0), v(K), averaged over draws.
fairlie_shapley <- function(X_g, w_g, X_ref, w_ref, alpha, beta, n_draws) {
  n_m     <- min(nrow(X_g), nrow(X_ref))
  eta_g   <- as.numeric(alpha + X_g   %*% beta)
  eta_ref <- as.numeric(alpha + X_ref %*% beta)
  phi_sum <- setNames(numeric(K), block_names)
  v0_sum  <- 0
  vK_sum  <- 0
  
  for (d in seq_len(n_draws)) {
    # Weighted samples of equal size from each profile
    s_g   <- sample.int(nrow(X_g),   n_m, replace = TRUE, prob = w_g)
    s_ref <- sample.int(nrow(X_ref), n_m, replace = TRUE, prob = w_ref)
    # Rank each sample by predicted adoption (random tie-break), pair by rank
    s_g   <- s_g[order(eta_g[s_g],     runif(n_m))]
    s_ref <- s_ref[order(eta_ref[s_ref], runif(n_m))]
    X_g_m   <- X_g[s_g,     , drop = FALSE]
    X_ref_m <- X_ref[s_ref, , drop = FALSE]
    
    v_vals <- numeric(nrow(subsets))
    for (r in seq_len(nrow(subsets))) {
      X_cf <- X_g_m
      for (bn in block_names[subsets[r, ]]) {
        X_cf[, blocks[[bn]]] <- X_ref_m[, blocks[[bn]]]
      }
      v_vals[r] <- mean(plogis(as.numeric(alpha + X_cf %*% beta)))
    }
    
    for (k in seq_len(K)) {
      for (r in seq_len(nrow(subsets))) {
        if (subsets[r, k]) next
        r2 <- row_with_block[r, k]
        phi_sum[k] <- phi_sum[k] + shapley_weights[r] * (v_vals[r2] - v_vals[r])
      }
    }
    v0_sum <- v0_sum + v_vals[1]
    vK_sum <- vK_sum + v_vals[nrow(subsets)]
  }
  
  list(phi = phi_sum / n_draws, v0 = v0_sum / n_draws, vK = vK_sum / n_draws)
}

# Reference households (used by the matching counterfactual below)
i_ref <- which(grp == ref_profile)

# ── Fairlie matching implementation ──────────────────────────────────────────

phi_match <- list()
v0_match  <- setNames(numeric(length(comparisons)), comparisons)
vK_match  <- setNames(numeric(length(comparisons)), comparisons)

for (g in comparisons) {
  i_g <- which(grp == g)
  fs  <- fairlie_shapley(X_g   = Xmat[i_g,   , drop = FALSE], w_g   = w[i_g],
                         X_ref = Xmat[i_ref, , drop = FALSE], w_ref = w[i_ref],
                         alpha = alpha, beta = beta, n_draws = n_match)
  phi_match[[g]] <- fs$phi
  v0_match[g]    <- fs$v0
  vK_match[g]    <- fs$vK
}

# ══════════════════════════════════════════════════════════════════════════════
# 5. FAIRLIE ROBUSTNESS WITH SEPARATE MODELS
# ══════════════════════════════════════════════════════════════════════════════
# Sections 3 and 4 impose one beta. Fairlie lets each profile keep its own and
# splits the gap into explained (different X, same coefficients) and unexplained
# (same X, different coefficients). The split needs a reference coefficient
# vector; we use only the pooled vector of the baseline logit, which does not
# privilege either profile (Model C in the paper).

profile_models <- list()
for (g in profiles) {
  d_g <- data_clean[grp == g, ]
  profile_models[[g]] <- suppressWarnings(
    glm(formula_baseline, family = binomial("logit"),
        data = d_g, weights = weight_final))
}

# Pooled (baseline logit) model: the single coefficient reference for Model C
pooled_model <- suppressWarnings(
  glm(formula_baseline, family = binomial("logit"),
      data = data_clean, weights = weight_final))

tbl_modelC <- data.frame()

for (g in comparisons) {
  
  i_A <- which(grp == ref_profile)   # reference group (A)
  i_B <- which(grp == g)            # compared group (B)
  
  # Own coefficients of each group give the total gap (reference-invariant)
  cf_A <- coef(profile_models[[ref_profile]])
  cf_B <- coef(profile_models[[g]])
  
  eta_A_A <- as.numeric(cf_A["(Intercept)"] + Xmat[i_A, , drop = FALSE] %*% cf_A[active_vars])
  eta_B_B <- as.numeric(cf_B["(Intercept)"] + Xmat[i_B, , drop = FALSE] %*% cf_B[active_vars])
  mu_A_A <- weighted.mean(plogis(eta_A_A), w = w[i_A])
  mu_B_B <- weighted.mean(plogis(eta_B_B), w = w[i_B])
  total  <- mu_A_A - mu_B_B
  
  # Explained component at the pooled (baseline logit) coefficient vector
  cf_r <- coef(pooled_model)
  eta_A_r <- as.numeric(cf_r["(Intercept)"] + Xmat[i_A, , drop = FALSE] %*% cf_r[active_vars])
  eta_B_r <- as.numeric(cf_r["(Intercept)"] + Xmat[i_B, , drop = FALSE] %*% cf_r[active_vars])
  mu_A_r <- weighted.mean(plogis(eta_A_r), w = w[i_A])
  mu_B_r <- weighted.mean(plogis(eta_B_r), w = w[i_B])
  
  explained   <- mu_A_r - mu_B_r   # same beta, different X
  unexplained <- total - explained
  
  tbl_modelC <- rbind(tbl_modelC, data.frame(
    comparison        = paste0("P", ref_profile, " vs P", g),
    coefficients_from = "pooled",
    total_pp          = round(total       * 100, 2),
    explained_pp      = round(explained   * 100, 2),
    unexplained_pp    = round(unexplained * 100, 2),
    explained_pct     = ifelse(abs(total) < 0.02, NA_real_,
                               round(100 * explained / total, 1)),
    stringsAsFactors = FALSE))
}

rownames(tbl_modelC) <- NULL

cat("\n-- 5. Fairlie decomposition with separate models, pooled reference (pp) --\n")
print(tbl_modelC)

# Closeness of the pooled-reference explained component to the Model A
# composition indicates the common-slope restriction is not driving the result.
# Small gaps produce unstable shares (suppressed as NA) and are not interpreted.


# ══════════════════════════════════════════════════════════════════════════════
# 5b. MODEL B: PROFILE-SPECIFIC AGE COEFFICIENT
# ══════════════════════════════════════════════════════════════════════════════
# Intermediate specification (Section 4.6): each profile has its own age
# coefficient; the other eleven characteristics keep common coefficients.
#   Pr(Y_i = 1 | X_i, g_i = g) = Lambda(alpha + delta_g + beta_age,g Age_i
#                                        + X_i,-age' beta_-age)
# The compositional component is evaluated at the reference profile's
# intercept AND age coefficient, as in Model A (Section 4.4):
#   composition_B = mu_ref - mu_g|ref, where mu_g|ref evaluates the households
#   of g with alpha, beta_age,ref and beta_-age.
# Because P1's age range is narrow, applying its age coefficient to the other
# profiles extrapolates; this is why Model A is the main specification.

age_var    <- "meanageadults"
other_vars <- setdiff(active_vars, age_var)
formula_B  <- as.formula(paste(outcome_var, "~ profile_gmm + profile_gmm:", age_var,
                               "+", paste(other_vars, collapse = " + ")))

m_B  <- suppressWarnings(glm(formula_B, family = binomial("logit"),
                             data = data_clean, weights = weight_final))
cf_B <- coef(m_B)

# Age coefficient of each profile (no age main effect: one slope per level)
beta_age <- sapply(profiles, function(g) {
  nm <- grep(paste0("^profile_gmm", g, ":", age_var, "$"), names(cf_B), value = TRUE)
  unname(cf_B[nm])
})
delta_B <- sapply(profiles, function(g) {
  nm <- paste0("profile_gmm", g)
  if (nm %in% names(cf_B)) unname(cf_B[nm]) else 0
})

X_other   <- as.matrix(data_clean[, other_vars])
age       <- data_clean[[age_var]]
eta_B_ref <- as.numeric(cf_B["(Intercept)"] + beta_age[ref_profile] * age +
                          X_other %*% cf_B[other_vars])
eta_B_own <- as.numeric(cf_B["(Intercept)"] + delta_B[grp] + beta_age[grp] * age +
                          X_other %*% cf_B[other_vars])

i_ref_B  <- which(grp == ref_profile)
mu_ref_B <- weighted.mean(plogis(eta_B_ref[i_ref_B]), w = w[i_ref_B])

tbl_modelB <- data.frame()
for (g in comparisons) {
  i_g <- which(grp == g)
  mu_g_ref <- weighted.mean(plogis(eta_B_ref[i_g]), w = w[i_g])
  mu_g_own <- weighted.mean(plogis(eta_B_own[i_g]), w = w[i_g])
  tbl_modelB <- rbind(tbl_modelB, data.frame(
    comparison     = paste0("P", ref_profile, " vs P", g),
    gap_pp         = round(100 * (mu_ref_B - mu_g_own), 2),
    composition_pp = round(100 * (mu_ref_B - mu_g_ref), 2),
    residual_pp    = round(100 * (mu_g_ref - mu_g_own), 2),
    stringsAsFactors = FALSE))
}

cat("\n-- 5b. Model B: profile-specific age coefficient (pp) --\n")
print(tbl_modelB)
cat(sprintf("   Age coefficients by profile: %s\n",
            paste(sprintf("P%s = %.4f", profiles, beta_age), collapse = ", ")))

remove(age_var, other_vars, formula_B, cf_B, beta_age, delta_B, X_other, age,
       eta_B_ref, eta_B_own, i_ref_B, mu_ref_B)


# ══════════════════════════════════════════════════════════════════════════════
# 6. STRATIFIED BOOTSTRAP
# ══════════════════════════════════════════════════════════════════════════════
# Households are resampled WITHIN each profile, everything is re-estimated, and
# the aggregate decomposition and Shapley contributions are recomputed (the
# Shapley step uses Fairlie rank matching, as in the point estimate).

boot_mat <- NULL

if (n_boot > 0) {
  
  cat(sprintf("\n-- 6. Stratified bootstrap (%d replications, Fairlie matching) --\n",
              n_boot))
  
  idx_by_profile <- split(seq_len(nrow(data_clean)), data_clean$profile_gmm)
  
  # Names of the quantities stored in each replication
  keys <- c()
  for (g in comparisons) {
    keys <- c(keys, paste0("comp_P", g), paste0("resid_P", g))
    for (bn in block_names) keys <- c(keys, paste0("phi_", bn, "_P", g))
  }
  
  boot_mat <- matrix(NA_real_, nrow = n_boot, ncol = length(keys),
                     dimnames = list(NULL, keys))
  
  t0 <- Sys.time()
  
  for (b in seq_len(n_boot)) {
    
    # Resample with replacement within each profile
    idx <- unlist(lapply(idx_by_profile,
                         function(ii) sample(ii, length(ii), replace = TRUE)))
    d_b <- data_clean[idx, ]
    
    # Whole computation inside a try in case a replication fails to converge
    result_b <- try({
      
      m_b <- suppressWarnings(
        glm(formula_common, family = binomial("logit"),
            data = d_b, weights = weight_final))
      cf_b    <- coef(m_b)
      alpha_b <- cf_b["(Intercept)"]
      beta_b  <- cf_b[active_vars]
      
      X_b   <- as.matrix(d_b[, active_vars])
      w_b   <- d_b$weight_final
      grp_b <- as.character(d_b$profile_gmm)
      
      eta_refalpha_b <- as.numeric(alpha_b + X_b %*% beta_b)
      
      delta_b <- setNames(numeric(length(profiles)), profiles)
      for (g in profiles) {
        coef_name <- paste0("profile_gmm", g)
        if (coef_name %in% names(cf_b)) delta_b[g] <- cf_b[coef_name]
      }
      eta_own_b <- eta_refalpha_b + delta_b[grp_b]
      
      i_ref_b  <- which(grp_b == ref_profile)
      mu_ref_b <- weighted.mean(plogis(eta_refalpha_b[i_ref_b]), w = w_b[i_ref_b])
      
      vals_b <- c()
      
      for (g in comparisons) {
        
        i_g   <- which(grp_b == g)
        X_g_b <- X_b[i_g, , drop = FALSE]
        w_g_b <- w_b[i_g]
        
        mu_g_own_b      <- weighted.mean(plogis(eta_own_b[i_g]),      w = w_g_b)
        mu_g_refalpha_b <- weighted.mean(plogis(eta_refalpha_b[i_g]), w = w_g_b)
        
        vals_b[paste0("comp_P",  g)] <- mu_ref_b - mu_g_refalpha_b
        vals_b[paste0("resid_P", g)] <- mu_g_refalpha_b - mu_g_own_b
        
        # Fairlie rank matching, fewer draws than in the point estimate
        fs_b  <- fairlie_shapley(X_g   = X_g_b,
                                 w_g   = w_g_b,
                                 X_ref = X_b[i_ref_b, , drop = FALSE],
                                 w_ref = w_b[i_ref_b],
                                 alpha = alpha_b, beta = beta_b,
                                 n_draws = n_match_b)
        phi_b <- fs_b$phi
        
        for (bn in block_names) {
          vals_b[paste0("phi_", bn, "_P", g)] <- phi_b[bn]
        }
      }
      
      vals_b
    }, silent = TRUE)
    
    # Store the replication if it did not fail
    if (!inherits(result_b, "try-error")) {
      boot_mat[b, names(result_b)] <- result_b
    }
    
    # Progress every 25 replications
    if (b %% 25 == 0) {
      elapsed   <- as.numeric(Sys.time() - t0, units = "secs")
      remaining <- (elapsed / b) * (n_boot - b)
      cat(sprintf("    [%d/%d] approx. %02d:%02d left\n",
                  b, n_boot, remaining %/% 60, round(remaining %% 60)))
    }
  }
}

# Count failed replications (guarded in case no bootstrap was run)
if (!is.null(boot_mat)) cat(sprintf("\n    Failed replications: %d of %d\n",
                                    sum(is.na(boot_mat[, 1])), n_boot))


# ══════════════════════════════════════════════════════════════════════════════
# 7. MAIN RESULT WITH BOOTSTRAP INTERVALS
# ══════════════════════════════════════════════════════════════════════════════

# Point estimates (Fairlie matching)
phi_point <- phi_match

tbl_ci <- data.frame()

for (g in comparisons) {
  
  comp_label <- paste0("P", ref_profile, " vs P", g)
  
  # ── Total composition ──
  key <- paste0("comp_P", g)
  if (!is.null(boot_mat) && key %in% colnames(boot_mat)) {
    q <- quantile(boot_mat[, key] * 100, probs = c(0.025, 0.975), na.rm = TRUE)
  } else {
    q <- c(NA_real_, NA_real_)
  }
  tbl_ci <- rbind(tbl_ci, data.frame(
    comparison = comp_label,
    component  = "Composition (total)",
    estimate   = tbl_decomp$composition_pp[tbl_decomp$comparison == comp_label],
    ci_low     = round(unname(q[1]), 2), ci_high = round(unname(q[2]), 2),
    stringsAsFactors = FALSE))
  
  # ── Profile residual ──
  key <- paste0("resid_P", g)
  if (!is.null(boot_mat) && key %in% colnames(boot_mat)) {
    q <- quantile(boot_mat[, key] * 100, probs = c(0.025, 0.975), na.rm = TRUE)
  } else {
    q <- c(NA_real_, NA_real_)
  }
  tbl_ci <- rbind(tbl_ci, data.frame(
    comparison = comp_label,
    component  = "Residual",
    estimate   = tbl_decomp$residual_pp[tbl_decomp$comparison == comp_label],
    ci_low     = round(unname(q[1]), 2), ci_high = round(unname(q[2]), 2),
    stringsAsFactors = FALSE))
  
  # ── Shapley contributions, block by block ──
  for (bn in block_names) {
    key <- paste0("phi_", bn, "_P", g)
    if (!is.null(boot_mat) && key %in% colnames(boot_mat)) {
      q <- quantile(boot_mat[, key] * 100, probs = c(0.025, 0.975), na.rm = TRUE)
    } else {
      q <- c(NA_real_, NA_real_)
    }
    tbl_ci <- rbind(tbl_ci, data.frame(
      comparison = comp_label,
      component  = paste0("  Shapley: ", bn),
      estimate   = round(phi_point[[g]][bn] * 100, 2),
      ci_low     = round(unname(q[1]), 2), ci_high = round(unname(q[2]), 2),
      stringsAsFactors = FALSE))
  }
}

rownames(tbl_ci) <- NULL

cat("\n-- 7. Decomposition with bootstrap intervals (pp) --\n")
print(tbl_ci)

# Read the headline comparison (the largest gap) from this table. A residual
# whose interval excludes zero is small-but-not-null: describe it as such.


# ══════════════════════════════════════════════════════════════════════════════
# 8. MANUSCRIPT TABLES
# ══════════════════════════════════════════════════════════════════════════════
# Body tables (4, 5, 6) and the Appendix D robustness table (D.3),
# in the same format as scripts 02 and 03.

# ── TABLE 4. Per-variable Wald tests ─────────────────────────────────────────
# Rows in canonical var_order, grouped by the six blocks of Section 3.

# Bonferroni adjustment over the family of characteristic-specific tests
n_tests_t4 <- nrow(tbl_wald)

tbl4_df <- tbl_wald %>%
  mutate(block = block_map[variable],
         block = factor(block, levels = block_order),
         variable = factor(variable, levels = var_order),
         p_bonf = pmin(1, p_raw * n_tests_t4)) %>%
  arrange(block, variable) %>%
  transmute(
    Characteristic = label,
    `Wald statistic` = sprintf("%.2f", W),
    df = df,
    `p-value` = sprintf("%.3f", p_value),
    `Bonferroni correction` = sprintf("%.3f", p_bonf)
  )

# Block positions, derived from the data so they are not fixed by hand
b4 <- tbl_wald %>%
  mutate(block = factor(block_map[variable], levels = block_order),
         variable = factor(variable, levels = var_order)) %>%
  arrange(block, variable) %>%
  pull(block)

idx4 <- lapply(levels(b4), function(b) which(b4 == b))
names(idx4) <- levels(b4)
idx4 <- idx4[sapply(idx4, length) > 0]

kbl_t4 <- kbl(tbl4_df,
              format    = "html",
              caption   = "Table 4. Wald tests of common slopes, by characteristic",
              align     = c("l", "r", "r", "r", "r"),
              row.names = FALSE) %>%
  kable_classic(full_width = FALSE, html_font = "Times New Roman") %>%
  kable_styling(font_size = 10) %>%
  row_spec(0, bold = TRUE)

for (b in names(idx4)) {
  kbl_t4 <- pack_rows(kbl_t4, b, min(idx4[[b]]), max(idx4[[b]]),
                      label_row_css = "background-color: #f0f0f0; font-weight: bold;")
}

# Note on the interactions left out by the support rule, built from the data,
# e.g. "Detached dwelling and children present (P2) lack support ..."
note_excl <- ""
if (nrow(excluded) > 0) {
  # Characteristics left out in each profile, in canonical order
  by_profile <- split(excluded$variable, excluded$profile)
  phrases <- vapply(names(by_profile), function(g) {
    vs <- unname(var_labels[intersect(var_order, by_profile[[g]])])
    # Lower-case initial for every label but the first
    if (length(vs) > 1)
      vs[-1] <- paste0(tolower(substr(vs[-1], 1, 1)), substring(vs[-1], 2))
    items <- if (length(vs) == 1) vs else
      paste(paste(vs[-length(vs)], collapse = ", "), "and", vs[length(vs)])
    sprintf("%s (P%s)", items, g)
  }, character(1))
  # Lower-case initial for every phrase but the first
  if (length(phrases) > 1)
    phrases[-1] <- paste0(tolower(substr(phrases[-1], 1, 1)), substring(phrases[-1], 2))
  n_excl <- nrow(excluded)
  note_excl <- paste0(" ", paste(phrases, collapse = "; "),
                      if (n_excl == 1) " lacks" else " lack",
                      " support and keep",
                      if (n_excl == 1) "s" else "",
                      " the common slope (Section 4.3), so ",
                      if (n_excl == 1) "its test has" else "their tests have",
                      " fewer degrees of freedom.")
}

kbl_t4 <- kbl_t4 %>%
  footnote(
    general = paste0("Each row is a joint Wald test (HC1) that the profile ",
                     "interactions of that characteristic are zero; the last ",
                     "column gives the Bonferroni-adjusted p-value.", note_excl),
    general_title     = "Notes:",
    footnote_as_chunk = TRUE
  )

print(kbl_t4)
save_table(kbl_t4, "Table_4.html")


# ── TABLE 5. Aggregate decomposition ─────────────────────────────────────────

tbl5_df <- tbl_decomp %>%
  transmute(
    Comparison         = comparison,
    `Difference (pp)`  = sprintf("%.2f", observed_gap),
    `Composition (pp)` = sprintf("%.2f", composition_pp),
    `Residual (pp)`    = sprintf("%.2f", residual_pp)
  )

kbl_t5 <- kbl(tbl5_df,
              format    = "html",
              caption   = paste0("Table 5. Decomposition of the adoption differences ",
                                 "(relative to P", ref_profile, ")"),
              align     = c("l", "r", "r", "r"),
              row.names = FALSE) %>%
  kable_classic(full_width = FALSE, html_font = "Times New Roman") %>%
  kable_styling(font_size = 10) %>%
  row_spec(0, bold = TRUE) %>%
  footnote(
    general = paste("Percentage points of adoption relative to the",
                    "highest-adopting profile. Composition is the part of the gap",
                    "accounted for by observable characteristics under a common",
                    "coefficient vector; the residual is what they leave",
                    "unexplained. Bootstrap intervals are in Table 6."),
    general_title     = "Notes:",
    footnote_as_chunk = TRUE
  )

print(kbl_t5)
save_table(kbl_t5, "Table_5.html")


# ── TABLE 6. Shapley attribution with intervals ──────────────────────────────
# One row per block, one column per comparison, with the estimate and its
# bootstrap interval. The first two rows are the totals.

block_labels <- c(tenure      = "Tenure",
                  income      = "Income",
                  education   = "Education",
                  climate_hdd = "Climatic exposure: heating (HDD)",
                  climate_cdd = "Climatic exposure: cooling (CDD)",
                  dwelling    = "Dwelling characteristics",
                  household   = "Household composition and location")

# Estimate and interval in a single cell
# Two decimals without a negative zero ("-0.00" prints as "0.00")
f2 <- function(x) sub("^-0\\.00$", "0.00", sprintf("%.2f", x))
cell_fmt <- function(est, lo, hi) {
  if (is.na(lo)) return(f2(est))
  sprintf("%s [%s, %s]", f2(est), f2(lo), f2(hi))
}

# Rows: the two totals first, then the seven blocks
rows_t6    <- c("Composition (total)", "Residual",
                paste0("  Shapley: ", block_names))
row_labels <- c("Composition (total)", "Residual",
                unname(block_labels[block_names]))

tbl6_df <- data.frame(Component = row_labels, stringsAsFactors = FALSE)

for (g in comparisons) {
  comp_label <- paste0("P", ref_profile, " vs P", g)
  col_vals <- character(length(rows_t6))
  for (k in seq_along(rows_t6)) {
    row_k <- tbl_ci[tbl_ci$comparison == comp_label &
                      tbl_ci$component == rows_t6[k], ]
    col_vals[k] <- if (nrow(row_k) == 0) "" else
      cell_fmt(row_k$estimate[1], row_k$ci_low[1], row_k$ci_high[1])
  }
  tbl6_df[[comp_label]] <- col_vals
}

kbl_t6 <- kbl(tbl6_df,
              format    = "html",
              caption   = paste0("Table 6. Shapley attribution of the ",
                                 "compositional component"),
              align     = c("l", rep("r", length(comparisons))),
              row.names = FALSE) %>%
  kable_classic(full_width = FALSE, html_font = "Times New Roman") %>%
  kable_styling(font_size = 10) %>%
  row_spec(0, bold = TRUE) %>%
  pack_rows("Aggregate decomposition", 1, 2,
            label_row_css = "background-color: #f0f0f0; font-weight: bold;") %>%
  pack_rows("Contribution of each block to the composition component",
            3, 2 + length(block_names),
            label_row_css = "background-color: #f0f0f0; font-weight: bold;") %>%
  footnote(
    general = paste("Percentage points of adoption, with 95% bootstrap intervals",
                    "from", n_boot, "stratified replications."),
    general_title     = "Notes:",
    footnote_as_chunk = TRUE
  )

print(kbl_t6)
save_table(kbl_t6, "Table_6.html")


# ── TABLE D.3. Composition and residual under alternative coefficient specs ──
# One block per comparison (observed difference in the group label); within each
# block, three rows: Model A (common coefficients), Model B (profile-specific
# age) and Model C (separate models by profile, Fairlie 2005, with the pooled
# baseline-logit vector as the reference). Composition + residual equals the
# observed difference in every row.

spec_labels <- c("Model A (common coefficients)",
                 "Model B (profile-specific age)",
                 "Model C (baseline logit vector)")

tblD3_df  <- data.frame()
gap_label <- character(0)   # one group label per comparison, with the difference

for (g in comparisons) {
  comp_label <- paste0("P", ref_profile, " vs P", g)
  a  <- tbl_decomp[tbl_decomp$comparison == comp_label, ]
  b  <- tbl_modelB[tbl_modelB$comparison == comp_label, ]
  cc <- tbl_modelC[tbl_modelC$comparison == comp_label, ]
  
  # Composition and residual of each specification (pp); C is the explained /
  # unexplained split at the pooled baseline-logit coefficient vector
  comp <- c(a$composition_pp,
            b$composition_pp,
            cc$explained_pp[cc$coefficients_from == "pooled"])
  resid <- c(a$residual_pp,
             b$residual_pp,
             cc$unexplained_pp[cc$coefficients_from == "pooled"])
  
  tblD3_df <- rbind(tblD3_df, data.frame(
    Specification      = spec_labels,
    `Composition (pp)` = sprintf("%.2f", comp),
    `Residual (pp)`    = sprintf("%.2f", resid),
    check.names = FALSE, stringsAsFactors = FALSE))
  
  gap_label[comp_label] <- sprintf("%s (observed difference: %.2f pp)",
                                   comp_label, a$observed_gap)
}

kbl_d3 <- kbl(tblD3_df,
              format    = "html",
              caption   = paste0("Table D.3. Composition and residual under ",
                                 "alternative coefficient specifications"),
              align     = c("l", "r", "r"),
              row.names = FALSE) %>%
  kable_classic(full_width = FALSE, html_font = "Times New Roman") %>%
  kable_styling(font_size = 10) %>%
  row_spec(0, bold = TRUE)

# One grouped block per comparison, three rows each
for (j in seq_along(comparisons)) {
  ini <- (j - 1) * length(spec_labels) + 1
  fin <- j * length(spec_labels)
  kbl_d3 <- pack_rows(kbl_d3, gap_label[j], ini, fin,
                      label_row_css = "background-color: #f0f0f0; font-weight: bold;")
}

kbl_d3 <- kbl_d3 %>%
  footnote(
    general = paste("Percentage points of adoption; in each row composition +",
                    "residual equals the observed difference."),
    general_title     = "Notes:",
    footnote_as_chunk = TRUE
  )

print(kbl_d3)
save_table(kbl_d3, "Table_D3.html")

# Cleanup
remove(tbl4_df, tbl5_df, tbl6_df, b4, idx4, kbl_t4, kbl_t5, kbl_t6,
       rows_t6, row_labels, cell_fmt, block_labels, n_tests_t4,
       tblD3_df, kbl_d3, spec_labels, gap_label)


# ══════════════════════════════════════════════════════════════════════════════
# 9. SAVE
# ══════════════════════════════════════════════════════════════════════════════

# Tables behind the numbers cited in the text go to output/ as CSV.
save_output(support,    "05_support_rule_cells.csv")
save_output(tbl_wald,   "05_wald_by_characteristic.csv")
save_output(tbl_decomp, "05_aggregate_decomposition.csv")
save_output(tbl_modelC, "05_modelC_fairlie_separate_models.csv")
save_output(tbl_modelB, "05_modelB_age_heterogeneous.csv")
save_output(tbl_ci,     "05_decomposition_with_intervals.csv")
save_output(data.frame(W = W_global, df = df_global, p = p_global,
                       estimated_interactions = nrow(int_terms),
                       excluded_interactions = nrow(excluded)),
            "05_wald_joint.csv")

save(W_global, df_global, p_global, support, int_terms, excluded,
     tbl_wald, tbl_decomp, tbl_modelC, tbl_ci, tbl_modelB,
     ref_profile, comparisons,
     file = file.path(rdata_dir, "05_decomposition.RData"))

cat(sprintf("\n=== [5] DONE. Time: %.1f min ===\n",
            as.numeric(Sys.time() - start_time, units = "mins")))
