# =============================================================================
# Bayesian Network Meta-Analysis: DAPT Antiplatelet Therapies After ACS
# Code for: Revisiting P2Y12 Inhibitors After Acute Coronary Syndromes -
#   Additional Insights from a Bayesian Network Meta-Analysis.
#   Can J Cardiol 2026 (accepted; CJC-D-26-02373R1).
# Updated from Kutcher et al. (2023) medRxiv:2023.08.12.23294021
#
# Trials added since the 2023 analysis (20 trials, 77,714 patients in total):
#   - TC4               (Kutcher et al., CMAJ 2025;197:E309-18)
#   - TUXEDO-2          (Bangalore et al., JAMA Cardiol 2026;11(4):369-377)
#   - SWITCH-SWEDEHEART (Omerovic et al., N Engl J Med 2026) - stepped-wedge
#     cluster RCT entered as "effective counts" that reproduce its published
#     cluster- and time-adjusted ORs (Section 1c-bis)
#
# Primary model:  brms | binomial cloglog | offset(log_time) | random effects
# Sensitivity analyses (SA1-SA9):
#   SA1 logit (log-OR) model           SA6 tight heterogeneity prior
#   SA2 PCI-only trials                SA7 wide heterogeneity prior
#   SA3 exclude TC4                    SA8 exclude SWITCH-SWEDEHEART
#   SA4 exclude TUXEDO-2               SA9 SWITCH-SWEDEHEART crude counts
#   SA5 exclude PRAGUE-18                  (MACE and major bleeding)
#
# Statistical workflow: Gelman et al. Bayesian Workflow (2020) arxiv:2011.01808
#   Step 1:  Data & reconciliation
#   Step 2:  Prior specification & justification
#   Step 3:  Prior predictive checks
#   Step 4:  Fake-data simulation & parameter recovery
#   Step 5:  Main model fitting
#   Step 6:  Convergence diagnostics
#   Step 7:  Posterior predictive checks
#   Step 8:  LOO-CV model comparison (RE vs FE)
#   Step 9:  Posterior summaries
#   Step 10: Sensitivity analyses and network consistency (node-splitting)
#   Step 11: Figures
#   Step 12: Discussion text, tables and CSV exports read by the manuscript
#   Step 13: Session info
#
# Run from the repository root:  source("dapt_bayesian_nma_oct2026.R")
# Outputs: outputs/{figures,tables,models}/oct/  (see README.md)
#
# Authors: James Brophy + original Kutcher et al. team
# Date: 2026-10-18 (analysis run); header updated 2026-10-28
# R version: >= 4.3.0 (run with R 4.6.1, brms 2.23.0, cmdstanr 0.9.0, CmdStan 2.39.0)
# =============================================================================

# -----------------------------------------------------------------------------
# 0. PACKAGES & SETUP
# -----------------------------------------------------------------------------
pkgs <- c(
  "brms",       # Bayesian regression via Stan
  "cmdstanr",   # CmdStan backend
  "tidyverse",  # Data wrangling + ggplot2
  "posterior",  # Posterior summaries
  "bayesplot",  # MCMC diagnostic plots
  "loo",        # LOO-CV
  "patchwork",  # Plot composition
  "scales",     # Plot formatting
  "gt",         # Publication-quality tables
  "webshot2",   # PDF export for gt (requires Chrome)
  "here"        # File paths
)



lapply(pkgs, library, character.only = TRUE)

# CmdStan backend — run cmdstanr::install_cmdstan() once if not installed
options(brms.backend = "cmdstanr")
options(mc.cores     = parallel::detectCores())
# Refit a cached model (file = ...) whenever its data, formula or priors change
options(brms.file_refit = "on_change")
set.seed(20260529)

# Output directories
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/models",  recursive = TRUE, showWarnings = FALSE)
for (d in c("outputs/figures/oct", "outputs/tables/oct", "outputs/models/oct"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# STEP 1: DATA
# =============================================================================
# -----------------------------------------------------------------------------
# 1a. RECONCILIATION NOTES
# -----------------------------------------------------------------------------
# Sources:
#   Kutcher 2023 medRxiv: 17 low-RoB trials, n=57,814
#   Maqsood 2026 JAMA Cardiol: 15 trials, n=48,904 (PCI-only, different RoB)
#   Bangalore 2026 JAMA Cardiol (TUXEDO-2): P vs T, DM+MVD, n=1800
#   Kutcher 2025 CMAJ (TC4): T vs C, North American, n=1005
#
# NEW TRIALS ADDED:
#
# TC4 (Kutcher 2025 CMAJ):
#   - Single-centre, Montreal; ACS undergoing PCI; T vs C; 12 months
#   - MACE (death+MI+stroke): 50/450 T, 64/555 C  [Table 2, primary outcome]
#   - Bleeding (GI bleed or haemorrhagic stroke): 20/450 T, 28/555 C
#   - Design: time-cluster RCT; LOO favoured pooled over hierarchical model
#     → treat as standard RCT arm in NMA; flag in sensitivity (SA3)
#   - Important: provides the only dedicated North American evidence
#   - RoB: low (investigator-initiated, CIHR-funded, no industry support)
#   - pci_only = TRUE (all patients had cardiac catheterisation)
#
# TUXEDO-2 (Bangalore 2026 JAMA Cardiol):
#   - Multicenter India; DM + multivessel PCI; P vs T; 12 months
#   - Hard MACE = "death, MI, or stroke" secondary endpoint (Table 2):
#     90/901 T, 77/899 P
#   - Bleeding (BARC major): 55/901 T, 43/899 P
#   - NOTE: diabetes + multivessel only; 21.8% chronic coronary syndrome
#     Declared COI: Bangalore received advisory fees from Viatris
#     (formerly Mylan, markets generic ticagrelor) — not declared in JAMA NMA
#
# PRAGUE-18:
#   - Using 1-year outcomes (Motovska 2018 JACC 71:371-381)
#   - MACE: 56/634 T, 44/596 P; Bleeding: 57/634 T, 42/596 P
#   - 39.1% crossover to clopidogrel by 1 year → sensitivity SA5
#
# SWITCH-SWEDEHEART (NEJM 2026):
#   - Omerovic 2026 NEJM (DOI 10.1056/NEJMoa2609057); Sweden; ACS undergoing PCI
#   - registry-based, open-label, stepped-wedge, cluster-randomized trial
#     (7 regions in 3 clusters; default T policy -> default P policy); 12 months
#   - MACE (death any cause, fatal/nonfatal MI, fatal/nonfatal stroke) [Table 2]:
#     1113/9444 T, 848/7651 P; adjusted OR 0.90 (0.77-1.06)
#   - Major bleeding, BROAD ICD-10 definition [Table 2]: 413/9444 T, 321/7651 P;
#     adjusted OR 0.80 (0.64-0.99). (Narrow definition: 225 T, 153 P)
#   - PRIMARY analysis uses effective counts reproducing the trial's adjusted
#     ORs (cluster + calendar time + covariates); see 1c-bis. Crude counts
#     ignore clustering and secular trends: crude MACE SE ~0.048 vs ~0.082
#     implied by the reported CI (~2.8x over-weighting) -> crude = SA9
#   - Appendix Table S9: sensitivity ORs 0.82-0.94 (GEE 0.88, CI unreliable
#     with 3 clusters; log-link RR 0.91 (0.78-1.05)) -> primary GLMM retained
#   - Policy (ITT) estimand: at discharge P-policy 63% prasugrel / 10% T /
#     27% clopidogrel; T-policy 77% T / 22% clopidogrel (Table S7)
#     -> dilution toward null
#   - Registry-ascertained, non-adjudicated outcomes; OAC + prior stroke retained
#   - Sensitivity: SA8 (exclude trial), SA9 (crude counts)
#
# GASECKA: 0 MACE events; retained (contributes to bleeding)
# REDUCE-MVI: small n=110 but retained 
# TRILOGY ACS: 0% PCI; excluded from PCI-only sensitivity
# ALPHEUS: 1-month follow-up < 6-month minimum → excluded
# Tang 2016, Wang 2016: "some concerns" RoB → excluded from main
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# 1b. HARD MACE DATA
# treat levels: C=clopidogrel (reference), T=ticagrelor, P=prasugrel
# -----------------------------------------------------------------------------
mace_raw <- tribble(
  ~study,           ~arm, ~treat, ~r,    ~n,     ~time_yr, ~pci_only,

  # ── T vs C ────────────────────────────────────────────────────────────────
  # TREAT (Berwanger 2019): STEMI + fibrinolysis; 12m
  "TREAT",           1,    "C",   164,  1886,   1.0,      FALSE,
  "TREAT",           2,    "T",   154,  1913,   1.0,      FALSE,
  # Gasecka 2020: STEMI/NSTEMI PCI; 0 MACE; 6m
  "Gasecka",         1,    "C",     0,    28,   0.5,      TRUE,
  "Gasecka",         2,    "T",     0,    27,   0.5,      TRUE,
  # POPular AGE (Gimbel 2020): NSTEMI >70yr; 12m
  "POPular AGE",     1,    "C",    79,   500,   1.0,      TRUE,
  "POPular AGE",     2,    "T",    81,   502,   1.0,      TRUE,
  # PHILO (Goto 2015): ACS Asia; 12m
  "PHILO",           1,    "C",    31,   400,   1.0,      TRUE,
  "PHILO",           2,    "T",    43,   401,   1.0,      TRUE,
  # He 2021: PCI diabetic; 6m
  "He",              1,    "C",     4,   133,   0.5,      TRUE,
  "He",              2,    "T",     5,   133,   0.5,      TRUE,
  # Mohareb 2020: PCI ACS; 12m
  "Mohareb",         1,    "C",     9,   472,   1.0,      TRUE,
  "Mohareb",         2,    "T",     4,   471,   1.0,      TRUE,
  # TICAKOREA (Park 2019): ACS invasive; 12m
  "TICAKOREA",       1,    "C",    31,   400,   1.0,      TRUE,
  "TICAKOREA",       2,    "T",    42,   400,   1.0,      TRUE,
  # PLATO (Wallentin 2009): ACS international; 12m
  "PLATO",           1,    "C",  1205,  9291,   1.0,      TRUE,
  "PLATO",           2,    "T",  1028,  9333,   1.0,      TRUE,
  # Wu 2021: PCI; 0 MACE; 12m
  "Wu",              1,    "C",     0,   174,   1.0,      TRUE,
  "Wu",              2,    "T",     0,   176,   1.0,      TRUE,
  # TC4 (Kutcher 2025 CMAJ): ACS PCI Montreal; T vs C; 12m  ← NEW
  # MACE (death+MI+stroke): T=50/450, C=64/555 [Table 2]
  "TC4",             1,    "C",    64,   555,   1.0,      TRUE,
  "TC4",             2,    "T",    50,   450,   1.0,      TRUE,

  # ── P vs C ────────────────────────────────────────────────────────────────
  # TRILOGY ACS (Roe 2012): UA/NSTEMI no PCI; 30m
  "TRILOGY ACS",     1,    "C",   854,  4663,   2.5,      FALSE,
  "TRILOGY ACS",     2,    "P",   808,  4663,   2.5,      FALSE,
  # PRASFIT-ACS (Saito 2014): scheduled PCI Japan; 12m
  "PRASFIT-ACS",     1,    "C",    89,   678,   1.0,      TRUE,
  "PRASFIT-ACS",     2,    "P",    75,   685,   1.0,      TRUE,
  # Elderly ACS II (Savonitto 2018): STEMI/NSTEMI PCI; 12m
  "Elderly ACS II",  1,    "C",    60,   730,   1.0,      TRUE,
  "Elderly ACS II",  2,    "P",    57,   713,   1.0,      TRUE,
  # TRITON-TIMI 38 (Wiviott 2007): scheduled PCI ACS; 15m
  "TRITON-TIMI 38",  1,    "C",   877,  6795,   1.25,     TRUE,
  "TRITON-TIMI 38",  2,    "P",   724,  6813,   1.25,     TRUE,
  # Yabe 2022: PCI Japan; 12m
  "Yabe",            1,    "C",     2,    39,   1.0,      TRUE,
  "Yabe",            2,    "P",     2,    37,   1.0,      TRUE,

  # ── P vs T ────────────────────────────────────────────────────────────────
  # PRAGUE-18 1-year (Motovska 2018 JACC): AMI PCI; 12m
  "PRAGUE-18",       1,    "T",    56,   634,   1.0,      TRUE,
  "PRAGUE-18",       2,    "P",    44,   596,   1.0,      TRUE,
  # ISAR-REACT 5 (Schupke 2019): ACS invasive; 12m
  "ISAR-REACT 5",    1,    "T",   208,  2012,   1.0,      TRUE,
  "ISAR-REACT 5",    2,    "P",   152,  2006,   1.0,      TRUE,
  # REDUCE-MVI (van der Hoeven 2020): STEMI PCI; 12m
  "REDUCE-MVI",      1,    "T",     7,    54,   1.0,      TRUE,
  "REDUCE-MVI",      2,    "P",     5,    56,   1.0,      TRUE,
  # TUXEDO-2 (Bangalore 2026 JAMA Cardiol): DM+MVD PCI; 12m  ← NEW
  "TUXEDO-2",        1,    "T",    90,   901,   1.0,      TRUE,
  "TUXEDO-2",        2,    "P",    77,   899,   1.0,      TRUE,
  # SWITCH-SWEDEHEART (Omerovic 2026 NEJM): AMI PCI; 12m  ← NEW
  "SWITCH-SWEDEHEART",        1,    "T",    1113,   9444,   1.0,      TRUE,
  "SWITCH-SWEDEHEART",        2,    "P",    848,    7651,   1.0,      TRUE
)

# -----------------------------------------------------------------------------
# 1c. MAJOR BLEEDING DATA
# -----------------------------------------------------------------------------
bleed_raw <- tribble(
  ~study,           ~arm, ~treat, ~r,    ~n,     ~time_yr, ~pci_only,

  # T vs C
  "TREAT",           1,    "C",    23,  1886,   1.0,      FALSE,
  "TREAT",           2,    "T",    20,  1913,   1.0,      FALSE,
  "Gasecka",         1,    "C",     1,    28,   0.5,      TRUE,
  "Gasecka",         2,    "T",     1,    27,   0.5,      TRUE,
  "POPular AGE",     1,    "C",     9,   500,   1.0,      TRUE,
  "POPular AGE",     2,    "T",    21,   502,   1.0,      TRUE,
  "PHILO",           1,    "C",    26,   400,   1.0,      TRUE,
  "PHILO",           2,    "T",    40,   401,   1.0,      TRUE,
  "He",              1,    "C",     2,   133,   0.5,      TRUE,
  "He",              2,    "T",     3,   133,   0.5,      TRUE,
  "Mohareb",         1,    "C",    22,   472,   1.0,      TRUE,
  "Mohareb",         2,    "T",    36,   471,   1.0,      TRUE,
  "TICAKOREA",       1,    "C",    16,   400,   1.0,      TRUE,
  "TICAKOREA",       2,    "T",    29,   400,   1.0,      TRUE,
  "PLATO",           1,    "C",   929,  9235,   1.0,      TRUE,
  "PLATO",           2,    "T",   961,  9333,   1.0,      TRUE,
  "Wu",              1,    "C",     2,   174,   1.0,      TRUE,
  "Wu",              2,    "T",     2,   176,   1.0,      TRUE,
  # TC4: bleeding = GI bleed or haemorrhagic stroke requiring admission
  # [Table 2: T=20/450, C=28/555]
  # NOTE: narrower bleeding definition than most trials (hospitalised only)
  "TC4",             1,    "C",    28,   555,   1.0,      TRUE,
  "TC4",             2,    "T",    20,   450,   1.0,      TRUE,

  # P vs C
  "TRILOGY ACS",     1,    "C",    48,  4617,   2.5,      FALSE,
  "TRILOGY ACS",     2,    "P",    58,  4623,   2.5,      FALSE,
  "PRASFIT-ACS",     1,    "C",    15,   678,   1.0,      TRUE,
  "PRASFIT-ACS",     2,    "P",    13,   685,   1.0,      TRUE,
  "Elderly ACS II",  1,    "C",    28,   713,   1.0,      TRUE,
  "Elderly ACS II",  2,    "P",    28,   730,   1.0,      TRUE,
  "TRITON-TIMI 38",  1,    "C",   111,  6716,   1.25,     TRUE,
  "TRITON-TIMI 38",  2,    "P",   146,  6741,   1.25,     TRUE,
  "Yabe",            1,    "C",     0,    39,   1.0,      TRUE,
  "Yabe",            2,    "P",     0,    37,   1.0,      TRUE,

  # P vs T
  "PRAGUE-18",       1,    "T",    57,   634,   1.0,      TRUE,
  "PRAGUE-18",       2,    "P",    42,   596,   1.0,      TRUE,
  "ISAR-REACT 5",    1,    "T",    95,  1989,   1.0,      TRUE,
  "ISAR-REACT 5",    2,    "P",    80,  1773,   1.0,      TRUE,
  "REDUCE-MVI",      1,    "T",     1,    51,   1.0,      TRUE,
  "REDUCE-MVI",      2,    "P",     1,    51,   1.0,      TRUE,
  "TUXEDO-2",        1,    "T",    55,   901,   1.0,      TRUE,
  "TUXEDO-2",        2,    "P",    43,   899,   1.0,      TRUE,
  "SWITCH-SWEDEHEART",            1,    "T",    413, 9444,   1.0,      TRUE,
  "SWITCH-SWEDEHEART",            2,    "P",    321, 7651,   1.0,      TRUE
)

# -----------------------------------------------------------------------------
# 1c-bis. SWITCH-SWEDEHEART: ADJUSTED EFFECTIVE COUNTS (PRIMARY ANALYSIS)
# -----------------------------------------------------------------------------
# The crude counts in 1b/1c are NOT used in the primary analysis. In a
# stepped-wedge design, policy is confounded with calendar time by
# construction (all clusters start on ticagrelor, end on prasugrel), and
# patients within a region are correlated. The trial's prespecified primary
# analysis (GLMM, logit link, cluster random intercept, restricted cubic
# spline for calendar time, baseline covariates, multiple imputation) gives:
#   MACE:     adjusted OR 0.90 (95% CI 0.77-1.06)   vs crude OR 0.93 (0.85-1.03)
#   Bleeding: adjusted OR 0.80 (95% CI 0.64-0.99)   vs crude OR 0.96 (0.83-1.11)
# The crude analysis is biased by secular trends and too precise
# (MACE log-OR SE 0.048 crude vs 0.082 adjusted).
#
# To enter the adjusted result into this arm-based binomial NMA, we convert it
# into integer "effective counts" (generalising the Cochrane Handbook
# effective-sample-size approach for cluster trials, sec 23.1.4):
#   1. keep the observed ticagrelor-arm risk p0;
#   2. set the prasugrel-arm risk p1 so that odds(p1)/odds(p0) = adjusted OR;
#   3. scale both arms by k (preserving the 9444:7651 allocation) so that the
#      binomial SE of the log-OR equals the SE implied by the reported 95% CI.
# 1/k is the implied design effect (~2.8 MACE, ~1.9 bleeding). The resulting
# counts reproduce OR 0.898 (0.766-1.054) for MACE and 0.802 (0.645-0.998)
# for bleeding; the corresponding cloglog HRs are 0.904 and 0.806.
# The crude counts are retained for sensitivity analysis SA9.
adj_to_counts <- function(r0, n0, n1, or, lo, hi) {
  p0     <- r0 / n0
  odds1  <- p0 / (1 - p0) * or
  p1     <- odds1 / (1 + odds1)
  se_adj <- (log(hi) - log(lo)) / (2 * qnorm(0.975))
  k      <- (1 / (n0 * p0 * (1 - p0)) + 1 / (n1 * p1 * (1 - p1))) / se_adj^2
  N0 <- round(k * n0); N1 <- round(k * n1)
  list(r0 = as.integer(round(N0 * p0)), n0 = as.integer(N0),
       r1 = as.integer(round(N1 * p1)), n1 = as.integer(N1),
       design_effect = 1 / k)
}

set_ss_counts <- function(df, cc) {
  df %>%
    mutate(
      r = case_when(study == "SWITCH-SWEDEHEART" & treat == "T" ~ cc$r0,
                    study == "SWITCH-SWEDEHEART" & treat == "P" ~ cc$r1,
                    TRUE ~ as.integer(r)),
      n = case_when(study == "SWITCH-SWEDEHEART" & treat == "T" ~ cc$n0,
                    study == "SWITCH-SWEDEHEART" & treat == "P" ~ cc$n1,
                    TRUE ~ as.integer(n))
    )
}

ss_mace_adj  <- adj_to_counts(r0 = 1113, n0 = 9444, n1 = 7651,
                              or = 0.90, lo = 0.77, hi = 1.06)
ss_bleed_adj <- adj_to_counts(r0 = 413,  n0 = 9444, n1 = 7651,
                              or = 0.80, lo = 0.64, hi = 0.99)

mace_raw_crude  <- mace_raw            # crude counts -> SA9
bleed_raw_crude <- bleed_raw
mace_raw  <- set_ss_counts(mace_raw,  ss_mace_adj)    # PRIMARY
bleed_raw <- set_ss_counts(bleed_raw, ss_bleed_adj)   # PRIMARY

cat("\nSWITCH-SWEDEHEART effective counts (primary):\n")
cat(sprintf("  MACE:  T %d/%d, P %d/%d (design effect %.2f)\n",
    ss_mace_adj$r0, ss_mace_adj$n0, ss_mace_adj$r1, ss_mace_adj$n1,
    ss_mace_adj$design_effect))
cat(sprintf("  Bleed: T %d/%d, P %d/%d (design effect %.2f)\n",
    ss_bleed_adj$r0, ss_bleed_adj$n0, ss_bleed_adj$r1, ss_bleed_adj$n1,
    ss_bleed_adj$design_effect))

# -----------------------------------------------------------------------------
# 1d. DATA PREPARATION
# KEY: log_time computed as column, used via offset() on RHS of formula
# Integer continuity correction for 0-event cells: r=0 → r=1, n → n+1
# Required: binomial family demands integer responses (0.5 would fail)
# Affected trials: Gasecka (0 MACE both arms) and Wu (0 MACE both arms)
# Effect is minimal: rates become 1/29 and 1/175 respectively
# Sensitivity: these trials contribute negligibly to efficacy estimates
# -----------------------------------------------------------------------------
prepare_nma_data <- function(df, pci_filter = FALSE) {
  d <- df
  if (pci_filter) d <- filter(d, pci_only == TRUE)

  d %>%
    mutate(
      treat    = factor(treat, levels = c("C", "T", "P")),
      study    = factor(study),
      log_time = log(time_yr),
      # Integer continuity correction for 0-event cells
      # 0 events -> 1 event, n -> n+1: gives rate 1/(n+1) near zero
      # Required because binomial family requires integer responses
      # For cloglog, near-zero rates are well-behaved (log-hazard -> -Inf)
      # For logit SA, this avoids log(0) issues
      r_cc     = as.integer(ifelse(r == 0L, 1L, r)),
      n_cc     = as.integer(ifelse(r == 0L, n + 1L, n))
    ) %>%
    arrange(study, arm)
}

mace_data  <- prepare_nma_data(mace_raw)
bleed_data <- prepare_nma_data(bleed_raw)
mace_pci   <- prepare_nma_data(mace_raw,  pci_filter = TRUE)
bleed_pci  <- prepare_nma_data(bleed_raw, pci_filter = TRUE)
mace_data_crude  <- prepare_nma_data(mace_raw_crude)    # SA9
bleed_data_crude <- prepare_nma_data(bleed_raw_crude)   # SA9

# Network summary
cat("\n=== NETWORK SUMMARY ===\n")
cat("Main analysis:\n")
# Patient count uses actual (crude) numbers, not SWITCH-SWEDEHEART effective n
cat("  MACE  — trials:", n_distinct(mace_data$study),
    "| patients:", sum(mace_raw_crude$n), "\n")
cat("  Bleed — trials:", n_distinct(bleed_data$study), "\n")
cat("PCI-only sensitivity:\n")
cat("  MACE  — trials:", n_distinct(mace_pci$study),
    "(TREAT + TRILOGY ACS excluded)\n")
cat("TC4 included:", "TC4" %in% mace_data$study, "\n")
cat("TUXEDO-2 included:", "TUXEDO-2" %in% mace_data$study, "\n")
cat("SWITCH-SWEDEHEART included:", "SWITCH-SWEDEHEART" %in% mace_data$study, "\n")

# =============================================================================
# STEP 2: PRIOR SPECIFICATION
# =============================================================================
# -----------------------------------------------------------------------------
# Model: binomial(cloglog); logclog(pi_ij) = alpha_i + beta_j + log(t_ij)
#   alpha_i = study-specific baseline log-hazard (reference = C arm)
#   beta_j  = population-level log-HR for treatment j vs C
#   tau     = between-study SD (heterogeneity)
# complementary log-log link is the best link
# exponentiate a difference in cloglog linear predictors between two treatments -> hazard ratio
# preferred over logit as it gives HR and not OR and 
# allows incorporation of the time offset
#
# PRIOR JUSTIFICATION:
#
# Intercept (baseline log-hazard for C arm):
#   Observed MACE rates 5–15%/year → cloglog(0.10) ≈ -2.25
#   Normal(-2.0, 1.0) covers event rates 1–30% — appropriately diffuse
#
# Treatment effects (log-HR vs C):
#   Clinically implausible for any antiplatelet HR < 0.4 or > 2.5
#   Normal(0, 0.5) places 95% mass on HR [0.37, 2.70] — weakly informative
#
# Between-study SD (tau):
#   Specified as normal(0, 0.5) for class = "sd" in brms.
#   brms automatically constrains sd parameters to be positive (lb = 0),
#   so normal(0, 0.5) on sd class IS a half-normal in the Stan code —
#   do NOT use the string "half_normal" which Stan does not recognise.
#   This puts 95% mass on tau < 0.98 (on log-HR scale).
#   Empirically, cardiovascular drug NMAs show tau ≈ 0.1–0.3
#   Tight SA: normal(0, 0.25) | Wide SA: normal(0, 1.0) — both on sd class
#
# TC4 PRIOR CONNECTION:
#   TC4 (Kutcher 2025 CMAJ) used our medRxiv BNMA as their "summary prior"
#   (RR_MACE = 0.95, 95%CrI 0.81–1.14). Now TC4 data enters the likelihood
#   directly, completing the Bayesian updating cycle.
# -----------------------------------------------------------------------------

priors_primary <- c(
  prior(normal(-2.0, 1.0),    class = "Intercept"),
  prior(normal(0, 0.5),       class = "b"),
  prior(normal(0, 0.5),  class = "sd")
)

priors_tight <- c(
  prior(normal(-2.0, 1.0),    class = "Intercept"),
  prior(normal(0, 0.5),       class = "b"),
  prior(normal(0, 0.25), class = "sd")
)

priors_wide <- c(
  prior(normal(-2.0, 1.0),    class = "Intercept"),
  prior(normal(0, 0.5),       class = "b"),
  prior(normal(0, 1.0),  class = "sd")
)

priors_logor <- c(
  prior(normal(-2.2, 1.0),    class = "Intercept"),
  prior(normal(0, 0.5),       class = "b"),
  prior(normal(0, 0.5),  class = "sd")
)

priors_fe <- c(
  prior(normal(-2.0, 1.0),    class = "Intercept"),
  prior(normal(0, 0.5),       class = "b")
)

# =============================================================================
# STEP 3: PRIOR PREDICTIVE CHECKS
# A prior predictive check asks: before seeing any data, what event rates 
# does your model consider plausible? 
# Fit the model with sample_prior = "only",
# which tells Stan to sample from the prior distributions and push them 
# through the likelihood to generate predicted event counts — 
# but to ignore the actual observed data entirely.
# =============================================================================
cat("\n=== STEP 3: PRIOR PREDICTIVE CHECKS ===\n")

# Use subset of well-powered trials for speed
#
ppc_data <- mace_data %>%
  filter(study %in% c("PLATO", "ISAR-REACT 5", "TRITON-TIMI 38", "TUXEDO-2"))

fit_prior_ppc <- brm(
  formula = r_cc | trials(n_cc) ~ treat + offset(log_time) + (treat | study),
  data    = ppc_data,
  family  = binomial(link = "cloglog"),
  prior   = priors_primary,
  sample_prior = "only",
  chains  = 4,
  iter    = 2000,
  warmup  = 1000,
  seed      = 20260529,
  save_pars = save_pars(all = TRUE),
  refresh = 0,
  file      = "outputs/models/fit_prior_ppc"
)

# Prior predictive event rates
pp_prior <- posterior_predict(fit_prior_ppc)
pp_rate  <- sweep(pp_prior, 2, ppc_data$n_cc, "/")
prob70 <- quantile2(pp_rate, c(0.05, 0.70))

p_ppc_prior <- ggplot(
  data.frame(rate = as.vector(pp_rate[, ppc_data$n_cc > 50])),
  aes(x = rate)
) +
  geom_histogram(bins = 60, fill = "#2166AC", alpha = 0.7, color = "white") +
  geom_vline(xintercept = c(0.05, 0.24),
             linetype = "dashed", color = "darkred", linewidth = 0.8) +
  scale_x_continuous(labels = percent, limits = c(0, 0.6),
                     name = "Predicted annual event rate") +
  scale_y_continuous(name = "Count") +
  annotate("text", x = 0.05, y = Inf, label = " 5%",
           vjust = 2, color = "darkred", size = 3.5) +
  annotate("text", x = 0.24, y = Inf, label = " 24%",
           vjust = 2, color = "darkred", size = 3.5) +
  labs(
    title    = "Prior Predictive Check: Predicted MACE Rates",
    subtitle = "Most mass (65%) falls between dashed lines"
  ) +
  theme_bw(base_size = 12)

ggsave("outputs/figures/oct/s_01_prior_predictive.pdf", p_ppc_prior,
       width = 7, height = 4.5)
ggsave("outputs/figures/oct/s_01_prior_predictive.png", p_ppc_prior,
       width = 7, height = 4.5)


# Prior HR summaries
pp_draws <- as_draws_df(fit_prior_ppc)
cat("Prior log-HR T vs C: mean =",
    round(mean(pp_draws$b_treatT), 3),
    "SD =", round(sd(pp_draws$b_treatT), 3),
    "→ HR median:", round(median(exp(pp_draws$b_treatT)), 2),
    "95%:", round(quantile(exp(pp_draws$b_treatT), 0.025), 2),
    "to", round(quantile(exp(pp_draws$b_treatT), 0.975), 2), "\n")
cat("Prior log-HR P vs C: mean =",
    round(mean(pp_draws$b_treatP), 3),
    "SD =", round(sd(pp_draws$b_treatP), 3),
    "→ HR median:", round(median(exp(pp_draws$b_treatP)), 2),
    "95%:", round(quantile(exp(pp_draws$b_treatP), 0.025), 2),
    "to", round(quantile(exp(pp_draws$b_treatP), 0.975), 2), "\n")

# =============================================================================
# STEP 4: FAKE-DATA SIMULATION & PARAMETER RECOVERY
# =============================================================================
cat("\n=== STEP 4: PARAMETER RECOVERY ===\n")

# Assume some true parameters and see if the model can recapture them
# Generate fake data from your model using those values
# Fit your model to the fake data
# Check whether the posterior credible intervals contain the true values
#

true_pars <- list(
  log_hr_T = -0.05,   # T vs C (small benefit ~ HR 0.95)
  log_hr_P = -0.15,   # P vs C (modest benefit ~ HR 0.86)
  tau       = 0.15,   # between-study SD
  baseline  = -2.20   # ≈ 11% annual MACE rate in C arm
)

set.seed(42)
n_sim_studies <- 20   # matches our 20-trial network
comparison_assign <- c(rep("TC", 10), rep("PC", 5), rep("PT", 5))

sim_data <- map_dfr(seq_len(n_sim_studies), function(s) {
  comp  <- comparison_assign[s]
  n_arm <- sample(200:10000, 1)
  t_yr  <- sample(c(0.5, 1.0, 1.25, 2.5), 1)
  re_T  <- rnorm(1, 0, true_pars$tau)
  re_P  <- rnorm(1, 0, true_pars$tau)

  arms <- switch(comp,
    TC = tibble(arm = 1:2, treat = c("C", "T"),
                lhr = c(0, true_pars$log_hr_T + re_T)),
    PC = tibble(arm = 1:2, treat = c("C", "P"),
                lhr = c(0, true_pars$log_hr_P + re_P)),
    PT = tibble(arm = 1:2, treat = c("T", "P"),
                lhr = c(true_pars$log_hr_T + re_T,
                        true_pars$log_hr_P + re_P))
  )

  arms %>%
    mutate(
      study    = factor(s),
      n_cc     = n_arm,
      log_time = log(t_yr),
      prob     = 1 - exp(-exp(true_pars$baseline + lhr + log_time)),
      r_cc     = as.integer(rbinom(n(), n_cc, prob)),
      treat    = factor(treat, levels = c("C", "T", "P")),
      pci_only = TRUE
    )
})

fit_recovery <- brm(
  formula = r_cc | trials(n_cc) ~ treat + offset(log_time) + (treat | study),
  data    = sim_data,
  family  = binomial(link = "cloglog"),
  prior   = priors_primary,
  chains  = 4,
  iter    = 3000,
  warmup  = 1500,
  control = list(adapt_delta = 0.95),
  seed      = 20260529,
  save_pars = save_pars(all = TRUE),
  file      = "outputs/models/oct/fit_recovery"
)

recov_draws <- as_draws_df(fit_recovery)
tau_cols    <- grep("^sd_study__treat", names(recov_draws), value = TRUE)

cat("\nParameter recovery (true → recovered):\n")
cat(sprintf("  log-HR T vs C: %.3f → %.3f (95%%CrI: %.3f to %.3f)\n",
    true_pars$log_hr_T,
    mean(recov_draws$b_treatT),
    quantile(recov_draws$b_treatT, 0.025),
    quantile(recov_draws$b_treatT, 0.975)))
cat(sprintf("  log-HR P vs C: %.3f → %.3f (95%%CrI: %.3f to %.3f)\n",
    true_pars$log_hr_P,
    mean(recov_draws$b_treatP),
    quantile(recov_draws$b_treatP, 0.025),
    quantile(recov_draws$b_treatP, 0.975)))
cat(sprintf("  tau:           %.3f → %.3f (95%%CrI: %.3f to %.3f)\n",
    true_pars$tau,
    mean(recov_draws[[tau_cols[1]]]),
    quantile(recov_draws[[tau_cols[1]]], 0.025),
    quantile(recov_draws[[tau_cols[1]]], 0.975)))

# =============================================================================
# STEP 5: MAIN MODEL FITTING
# =============================================================================
cat("\n=== STEP 5: MAIN MODEL FITTING ===\n")

mcmc_ctrl <- list(adapt_delta = 0.95, max_treedepth = 12)

# Shared formula — note offset on RHS throughout
nma_formula_cloglog <- bf(
  r_cc | trials(n_cc) ~ treat + offset(log_time) + (treat | study)
)
nma_formula_logit <- bf(
  r_cc | trials(n_cc) ~ treat + (treat | study)   # no time offset for logOR
)
nma_formula_fe <- bf(
  r_cc | trials(n_cc) ~ treat + offset(log_time) + (1 | study)
)

# To avoid using cached results uncomment the command below 
# file.remove(list.files("outputs/models", pattern="\\.rds$", full.names=TRUE))

# ── 5a. PRIMARY: MACE, cloglog RE ────────────────────────────────────────────
fit_mace_re <- brm(
  formula = nma_formula_cloglog,
  data    = mace_data,
  family  = binomial(link = "cloglog"),
  prior   = priors_primary,
  chains  = 4, iter = 4000, warmup = 2000,
  control   = mcmc_ctrl, seed = 20260529,
  save_pars = save_pars(all = TRUE),
  file      = "outputs/models/oct/fit_mace_re"
)

# ── 5b. MACE, cloglog FE (LOO comparison) ────────────────────────────────────
fit_mace_fe <- brm(
  formula = nma_formula_fe,
  data    = mace_data,
  family  = binomial(link = "cloglog"),
  prior   = priors_fe,
  chains  = 4, iter = 4000, warmup = 2000,
  control   = mcmc_ctrl, seed = 20260529,
  save_pars = save_pars(all = TRUE),
  file      = "outputs/models/oct/fit_mace_fe"
)

# ── 5c. BLEEDING, cloglog RE ──────────────────────────────────────────────────
fit_bleed_re <- brm(
  formula = nma_formula_cloglog,
  data    = bleed_data,
  family  = binomial(link = "cloglog"),
  prior   = priors_primary,
  chains  = 4, iter = 4000, warmup = 2000,
  control   = mcmc_ctrl, seed = 20260529,
  save_pars = save_pars(all = TRUE),
  file      = "outputs/models/oct/fit_bleed_re"
)

# =============================================================================
# STEP 6: CONVERGENCE DIAGNOSTICS
# =============================================================================
cat("\n=== STEP 6: CONVERGENCE DIAGNOSTICS ===\n")

check_convergence <- function(fit, name) {
  rhat_max <- max(rhat(fit), na.rm = TRUE)
  ess_min  <- min(ess_bulk(fit), na.rm = TRUE)
  divs     <- sum(
    nuts_params(fit)$Value[nuts_params(fit)$Parameter == "divergent__"]
  )

  status <- ifelse(rhat_max < 1.01 & ess_min > 400 & divs == 0, "OK", "WARN")
  cat(sprintf("[%s] %s | Rhat_max=%.4f | ESS_min=%.0f | Divergences=%d\n",
              status, name, rhat_max, ess_min, divs))

  if (rhat_max > 1.01)
    warning(name, ": R-hat > 1.01 — check traceplots and increase iterations")
  if (ess_min < 400)
    warning(name, ": ESS < 400 — consider more iterations")
  if (divs > 0)
    warning(name, ": Divergences detected — increase adapt_delta")
}

check_convergence(fit_mace_re,  "MACE RE (primary)")
check_convergence(fit_mace_fe,  "MACE FE")
check_convergence(fit_bleed_re, "Bleeding RE")

# Trace plots — primary model key parameters
p_trace <- mcmc_trace(
  as_draws_array(fit_mace_re),
  pars = c("b_Intercept", "b_treatT", "b_treatP",
           "sd_study__treatT", "sd_study__treatP")
) +
  labs(title = "Trace plots: MACE RE model (primary)") +
  theme_bw(base_size = 9)
ggsave("outputs/figures/oct/02_traceplots.pdf", p_trace, width = 10, height = 7)

# Pairs plot
p_pairs <- mcmc_pairs(
  as_draws_array(fit_mace_re),
  pars = c("b_treatT", "b_treatP",
           "sd_study__treatT", "sd_study__treatP"),
  off_diag_args = list(size = 0.4, alpha = 0.25)
)
ggsave("outputs/figures/oct/03_pairs.pdf", p_pairs, width = 8, height = 8)

# =============================================================================
# STEP 7: POSTERIOR PREDICTIVE CHECKS
# =============================================================================
cat("\n=== STEP 7: POSTERIOR PREDICTIVE CHECKS ===\n")

ppc_plot <- function(fit, obs_data, title_str) {
  pp <- posterior_predict(fit)

  # Proportion of observed values inside 95% CrI
  lo <- apply(pp, 2, quantile, 0.025)
  hi <- apply(pp, 2, quantile, 0.975)
  coverage <- mean(obs_data$r_cc >= lo & obs_data$r_cc <= hi)
  cat(title_str, "— 95% CrI coverage:",
      round(coverage * 100, 1), "% (target ~95%)\n")

  pp_summary <- tibble(
    study     = obs_data$study,
    treat     = obs_data$treat,
    obs_rate  = obs_data$r_cc / obs_data$n_cc,
    pred_med  = apply(pp, 2, median) / obs_data$n_cc,
    pred_lo   = lo / obs_data$n_cc,
    pred_hi   = hi / obs_data$n_cc,
    label     = paste0(study, " (", treat, ")")
  )

  ggplot(pp_summary, aes(x = reorder(label, obs_rate))) +
    geom_point(aes(y = obs_rate), color = "black", size = 2) +
    geom_point(aes(y = pred_med), color = "#2166AC",
               size = 2, shape = 17) +
    geom_errorbar(aes(ymin = pred_lo, ymax = pred_hi),
                  color = "#2166AC", width = 0.3, alpha = 0.5) +
    scale_y_continuous(labels = percent) +
    coord_flip() +
    labs(
      title    = title_str,
      subtitle = "Black = observed; Blue = posterior median + 95% CrI",
      x = NULL, y = "Event rate"
    ) +
    theme_bw(base_size = 8)
}

p_ppc_mace  <- ppc_plot(fit_mace_re,  mace_data,  "PPC: MACE")
p_ppc_bleed <- ppc_plot(fit_bleed_re, bleed_data, "PPC: Major Bleeding")
ggsave("outputs/figures/oct/04_ppc_mace.pdf",  p_ppc_mace,  width = 8, height = 13)
ggsave("outputs/figures/oct/05_ppc_bleed.pdf", p_ppc_bleed, width = 8, height = 13)

# =============================================================================
# STEP 8: LOO-CV MODEL COMPARISON
# =============================================================================
cat("\n=== STEP 8: LOO-CV MODEL COMPARISON ===\n")

fit_mace_re <- add_criterion(fit_mace_re, "loo", moment_match = FALSE)
# Check Pareto-k diagnostics; rerun with reloo=TRUE if k > 0.7
pareto_k_mace <- loo(fit_mace_re)$diagnostics$pareto_k
cat("MACE RE: Pareto-k > 0.7:", sum(pareto_k_mace > 0.7), "observations\n")
cat("MACE RE: Pareto-k > 0.5:", sum(pareto_k_mace > 0.5), "observations\n")
# If any k > 0.7: uncomment and run exact LOO for those observations:
# fit_mace_re <- add_criterion(fit_mace_re, "loo", reloo = TRUE,
#                              moment_match = FALSE)
fit_mace_fe <- add_criterion(fit_mace_fe, "loo", moment_match = FALSE)
loo_comp    <- loo_compare(fit_mace_re, fit_mace_fe)

cat("\nLOO comparison (positive elpd_diff favours first model listed):\n")
print(loo_comp)
cat("Rule: |elpd_diff| > 4 SE indicates meaningful preference.\n")
write_csv(
  tibble(model     = rownames(loo_comp),
         elpd_diff = loo_comp[, "elpd_diff"],
         se_diff   = loo_comp[, "se_diff"]),
  "outputs/tables/oct/06_loo.csv"
)
cat("RE model expected to be preferred given observed heterogeneity.\n")

# =============================================================================
# STEP 9: POSTERIOR SUMMARIES 
# =============================================================================

# ROPE: HR [0.9, 1.11] — matches Kutcher 2023 and TC4 paper thresholds
ROPE_lo <- log(0.9)
ROPE_hi <- log(1.11)

extract_posterior <- function(fit, label) {
  d <- as_draws_df(fit)

  # All three pairwise comparisons
  log_hr <- list(
    "T vs C" = d$b_treatT,
    "P vs C" = d$b_treatP,
    "P vs T" = d$b_treatP - d$b_treatT
  )

  map_dfr(names(log_hr), function(comp) {
    lh <- log_hr[[comp]]
    hr <- exp(lh)
    tibble(
      outcome          = label,
      comparison       = comp,
      n_draws          = length(lh),
      median_HR        = median(hr),
      mean_HR          = mean(hr),
      lo95_HR          = quantile(hr, 0.025),
      hi95_HR          = quantile(hr, 0.975),
      lo90_HR          = quantile(hr, 0.05),
      hi90_HR          = quantile(hr, 0.95),
      pr_any_benefit   = mean(lh < 0),
      pr_clin_benefit  = mean(lh < ROPE_lo),   # HR < 0.9
      pr_rope          = mean(lh >= ROPE_lo & lh <= ROPE_hi),
      pr_clin_harm     = mean(lh > ROPE_hi)    # HR > 1.11
    )
  })
}

results_mace  <- extract_posterior(fit_mace_re,  "MACE (cloglog RE, primary)")
results_bleed <- extract_posterior(fit_bleed_re, "Major Bleeding (cloglog RE)")

cat("\n=== MAIN RESULTS ===\n")
print(
  bind_rows(results_mace, results_bleed) %>%
    mutate(across(c(median_HR, lo95_HR, hi95_HR), ~round(., 3)),
           across(starts_with("pr_"), ~sprintf("%.1f%%", . * 100))) %>%
    select(outcome, comparison, median_HR, lo95_HR, hi95_HR,
           pr_clin_benefit, pr_rope, pr_clin_harm),
  n = Inf
)

# =============================================================================
# STEP 10: SENSITIVITY ANALYSES
# =============================================================================
cat("\n=== STEP 10: SENSITIVITY ANALYSES ===\n")

# ── SA1: log-OR (logit link, no time offset) ──────────────────────────────────
fit_sa1_logor <- brm(
  formula = nma_formula_logit,
  data    = mace_data,
  family  = binomial(link = "logit"),
  prior   = priors_logor,
  chains  = 4, iter = 4000, warmup = 2000,
  control   = mcmc_ctrl, seed = 20260529,
  save_pars = save_pars(all = TRUE),
  file      = "outputs/models/oct/fit_sa1_logor"
)

# ── SA2: PCI-only trials ──────────────────────────────────────────────────────
fit_sa2_pci <- brm(
  formula = nma_formula_cloglog,
  data    = mace_pci,
  family  = binomial(link = "cloglog"),
  prior   = priors_primary,
  chains  = 4, iter = 4000, warmup = 2000,
  control   = mcmc_ctrl, seed = 20260529,
  save_pars = save_pars(all = TRUE),
  file      = "outputs/models/oct/fit_sa2_pci"
)

# ── SA3: Exclude TC4 (assess cluster-RCT design sensitivity) ─────────────────
fit_sa3_no_tc4 <- brm(
  formula = nma_formula_cloglog,
  data    = filter(mace_data, study != "TC4"),
  family  = binomial(link = "cloglog"),
  prior   = priors_primary,
  chains  = 4, iter = 4000, warmup = 2000,
  control   = mcmc_ctrl, seed = 20260529,
  save_pars = save_pars(all = TRUE),
  file      = "outputs/models/oct/fit_sa3_no_tc4"
)

# ── SA4: Exclude TUXEDO-2 (selected DM+MVD population) ───────────────────────
fit_sa4_no_tux <- brm(
  formula = nma_formula_cloglog,
  data    = filter(mace_data, study != "TUXEDO-2"),
  family  = binomial(link = "cloglog"),
  prior   = priors_primary,
  chains  = 4, iter = 4000, warmup = 2000,
  control   = mcmc_ctrl, seed = 20260529,
  save_pars = save_pars(all = TRUE),
  file      = "outputs/models/oct/fit_sa4_no_tux"
)

# ── SA5: Exclude PRAGUE-18 (39% crossover) ────────────────────────────────────
fit_sa5_no_pra <- brm(
  formula = nma_formula_cloglog,
  data    = filter(mace_data, study != "PRAGUE-18"),
  family  = binomial(link = "cloglog"),
  prior   = priors_primary,
  chains  = 4, iter = 4000, warmup = 2000,
  control   = mcmc_ctrl, seed = 20260529,
  save_pars = save_pars(all = TRUE),
  file      = "outputs/models/oct/fit_sa5_no_pra"
)

# ── SA6: Tight tau prior ──────────────────────────────────────────────────────
fit_sa6_tight <- brm(
  formula = nma_formula_cloglog,
  data    = mace_data,
  family  = binomial(link = "cloglog"),
  prior   = priors_tight,
  chains  = 4, iter = 4000, warmup = 2000,
  control   = mcmc_ctrl, seed = 20260529,
  save_pars = save_pars(all = TRUE),
  file      = "outputs/models/oct/fit_sa6_tight"
)

# ── SA7: Wide tau prior ────────────────────────────────────────────────────────
fit_sa7_wide <- brm(
  formula = nma_formula_cloglog,
  data    = mace_data,
  family  = binomial(link = "cloglog"),
  prior   = priors_wide,
  chains  = 4, iter = 4000, warmup = 2000,
  control   = mcmc_ctrl, seed = 20260529,
  save_pars = save_pars(all = TRUE),
  file      = "outputs/models/oct/fit_sa7_wide"
)

# ── SA8: Exclude SWITCH-SWEDEHEART (cluster stepped wedge design) ────────────────────────────────────
fit_sa8_no_ss <- brm(
  formula = nma_formula_cloglog,
  data    = filter(mace_data, study != "SWITCH-SWEDEHEART"),
  family  = binomial(link = "cloglog"),
  prior   = priors_primary,
  chains  = 4, iter = 4000, warmup = 2000,
  control   = mcmc_ctrl, seed = 20260529,
  save_pars = save_pars(all = TRUE),
  file      = "outputs/models/oct/fit_sa8_no_ss"
)

# ── SA9: SWITCH-SWEDEHEART crude (unadjusted) counts ─────────────────────────
# Primary analysis uses effective counts reproducing the trial's adjusted OR
# (see 1c-bis). SA9 replaces them with the crude arm-level counts.
fit_sa9_crude <- brm(
  formula = nma_formula_cloglog,
  data    = mace_data_crude,
  family  = binomial(link = "cloglog"),
  prior   = priors_primary,
  chains  = 4, iter = 4000, warmup = 2000,
  control   = mcmc_ctrl, seed = 20260529,
  save_pars = save_pars(all = TRUE),
  file      = "outputs/models/oct/fit_sa9_crude"
)

# SA9 for bleeding (reported separately; Table 3 is MACE only)
fit_sa9_bleed_crude <- brm(
  formula = nma_formula_cloglog,
  data    = bleed_data_crude,
  family  = binomial(link = "cloglog"),
  prior   = priors_primary,
  chains  = 4, iter = 4000, warmup = 2000,
  control   = mcmc_ctrl, seed = 20260529,
  save_pars = save_pars(all = TRUE),
  file      = "outputs/models/oct/fit_sa9_bleed_crude"
)

# Convergence checks for all sensitivity models
walk2(
  list(fit_sa1_logor, fit_sa2_pci, fit_sa3_no_tc4,
       fit_sa4_no_tux, fit_sa5_no_pra, fit_sa6_tight, fit_sa7_wide,
       fit_sa8_no_ss, fit_sa9_crude, fit_sa9_bleed_crude),
  c("SA1 logOR", "SA2 PCI-only", "SA3 excl TC4",
    "SA4 excl TUXEDO-2", "SA5 excl PRAGUE-18",
    "SA6 tight tau", "SA7 wide tau", "SA8 excl SWITCH-SWEDEHEART",
    "SA9 SS crude (MACE)", "SA9 SS crude (bleeding)"),
  check_convergence
)

# Collect all sensitivity results
sa_results <- bind_rows(
  extract_posterior(fit_mace_re,    "Primary (cloglog RE)"),
  extract_posterior(fit_sa1_logor,  "SA1: log-OR"),
  extract_posterior(fit_sa2_pci,    "SA2: PCI-only"),
  extract_posterior(fit_sa3_no_tc4, "SA3: excl. TC4"),
  extract_posterior(fit_sa4_no_tux, "SA4: excl. TUXEDO-2"),
  extract_posterior(fit_sa5_no_pra, "SA5: excl. PRAGUE-18"),
  extract_posterior(fit_sa6_tight,  "SA6: tight tau"),
  extract_posterior(fit_sa7_wide,   "SA7: wide tau"),
  extract_posterior(fit_sa8_no_ss,  "SA8: excl. SWITCH-SWEDEHEART"),
  extract_posterior(fit_sa9_crude,  "SA9: SWITCH-SWEDEHEART crude counts")
)

results_bleed_sa9 <- extract_posterior(fit_sa9_bleed_crude,
                                       "Bleeding SA9: SWITCH-SWEDEHEART crude")
cat("\nBleeding, SA9 (SWITCH-SWEDEHEART crude counts):\n")
print(results_bleed_sa9 %>%
        mutate(across(c(median_HR, lo95_HR, hi95_HR), ~round(., 3))) %>%
        select(comparison, median_HR, lo95_HR, hi95_HR,
               pr_clin_benefit, pr_rope, pr_clin_harm))
write_csv(results_bleed_sa9, "outputs/tables/oct/03b_sensitivity_bleed_sa9.csv")

# =============================================================================
# STEP 10b: NETWORK CONSISTENCY (node-splitting, P vs T comparison)
# =============================================================================
cat("\n=== STEP 10b: NETWORK CONSISTENCY ===\n")

# Direct P vs T trials: PRAGUE-18, ISAR-REACT 5, REDUCE-MVI, TUXEDO-2,
#   SWITCH-SWEDEHEART
# NOTE: "indirect" below is the full NETWORK estimate (direct + indirect),
#   not a pure indirect estimate. A pure indirect estimate requires refitting
#   the network without the 5 P vs T trials.
direct_PT <- mace_data %>%
  filter(study %in% c("PRAGUE-18", "ISAR-REACT 5",
                      "REDUCE-MVI", "TUXEDO-2", "SWITCH-SWEDEHEART")) %>%
  mutate(treat = factor(treat, levels = c("T", "P")))  # T as reference

fit_direct_PT <- brm(
  formula = r_cc | trials(n_cc) ~ treat + offset(log_time) + (treat | study),
  data    = direct_PT,
  family  = binomial(link = "cloglog"),
  prior   = c(prior(normal(-2.0, 1.0), class = "Intercept"),
              prior(normal(0, 0.5),    class = "b"),
              prior(normal(0, 0.5), class = "sd")),
  chains  = 4, iter = 4000, warmup = 2000,
  control   = mcmc_ctrl, seed = 20260529,
  save_pars = save_pars(all = TRUE),
  file      = "outputs/models/oct/fit_direct_PT"
)

# Indirect P vs T from network model
d_main <- as_draws_df(fit_mace_re)
indirect_PT_loghr <- d_main$b_treatP - d_main$b_treatT

d_direct <- as_draws_df(fit_direct_PT)
direct_PT_loghr  <- d_direct$b_treatP  # reference is T in this model

cat("\nNode-splitting: P vs T\n")
cat(sprintf("  Direct   HR: %.3f (95%%CrI %.3f–%.3f)\n",
    exp(median(direct_PT_loghr)),
    exp(quantile(direct_PT_loghr, 0.025)),
    exp(quantile(direct_PT_loghr, 0.975))))
cat(sprintf("  Indirect HR: %.3f (95%%CrI %.3f–%.3f)\n",
    exp(median(indirect_PT_loghr)),
    exp(quantile(indirect_PT_loghr, 0.025)),
    exp(quantile(indirect_PT_loghr, 0.975))))
cat(sprintf("  IF (log scale): %.3f\n",
    median(direct_PT_loghr) - median(indirect_PT_loghr)))

# NOTE: 'indirect_PT_loghr' above is the full NETWORK estimate. A pure
# indirect P vs T estimate uses only the P vs C and T vs C trials:
fit_indirect_PT <- brm(
  formula = nma_formula_cloglog,
  data    = filter(mace_data, !study %in% c("PRAGUE-18", "ISAR-REACT 5",
                   "REDUCE-MVI", "TUXEDO-2", "SWITCH-SWEDEHEART")),
  family  = binomial(link = "cloglog"),
  prior   = priors_primary,
  chains  = 4, iter = 4000, warmup = 2000,
  control   = mcmc_ctrl, seed = 20260529,
  save_pars = save_pars(all = TRUE),
  file      = "outputs/models/oct/fit_indirect_PT"
)
d_ind <- as_draws_df(fit_indirect_PT)
pure_indirect_PT_loghr <- d_ind$b_treatP - d_ind$b_treatT
n_min <- min(length(direct_PT_loghr), length(pure_indirect_PT_loghr))
incons_loghr <- direct_PT_loghr[seq_len(n_min)] -
                pure_indirect_PT_loghr[seq_len(n_min)]

summ_hr <- function(lh, label) tibble(
  estimate  = label,
  median_HR = exp(median(lh)),
  lo95_HR   = exp(quantile(lh, 0.025)),
  hi95_HR   = exp(quantile(lh, 0.975)),
  pr_clin_benefit = mean(lh < ROPE_lo)
)
nodesplit <- bind_rows(
  summ_hr(direct_PT_loghr,        "direct"),
  summ_hr(pure_indirect_PT_loghr, "indirect"),
  summ_hr(indirect_PT_loghr,      "network"),
  summ_hr(incons_loghr,           "ratio_direct_indirect")
)
print(nodesplit)
write_csv(nodesplit, "outputs/tables/oct/07_nodesplit_PT.csv")

# =============================================================================
# STEP 11: FIGURES
# =============================================================================
cat("\n=== STEP 11: FIGURES ===\n")

# Colour palette consistent throughout
col_TC <- "#E31A1C"
col_PC <- "#1F78B4"
col_PT <- "#33A02C"
comp_colours <- c("T vs C" = col_TC, "P vs C" = col_PC, "P vs T" = col_PT)

# ── Figure 1: Main forest plot ────────────────────────────────────────────────
make_forest <- function(results_df, title_str) {
  d <- results_df %>%
    mutate(
      comparison = factor(comparison, levels = c("T vs C","P vs C","P vs T")),
      sig_lo = lo95_HR < 1,
      sig_hi = hi95_HR > 1,
      sig    = !(sig_lo & sig_hi),  # significant if CI doesn't cross 1
      label  = sprintf("%.2f (%.2f–%.2f)", median_HR, lo95_HR, hi95_HR)
    )

  max_x <- max(d$hi95_HR) * 1.5

  ggplot(d, aes(x = median_HR, y = comparison)) +
    annotate("rect", xmin = 0.9, xmax = 1.11,
             ymin = -Inf, ymax = Inf, fill = "#4393C3", alpha = 0.10) +
    geom_vline(xintercept = 1,    linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(0.9, 1.11), linetype = "dotted",
               color = "#4393C3", linewidth = 0.6) +
    geom_errorbarh(aes(xmin = lo90_HR, xmax = hi90_HR, color = comparison),
                   height = 0, linewidth = 2.0, alpha = 0.55) +
    geom_errorbarh(aes(xmin = lo95_HR, xmax = hi95_HR, color = comparison),
                   height = 0.22, linewidth = 0.8) +
    geom_point(aes(color = comparison, shape = sig), size = 4) +
    geom_text(aes(label = label), x = max_x * 0.62,
              hjust = 0, size = 3.0, color = "grey20") +
    scale_color_manual(values = comp_colours, guide = "none") +
    scale_shape_manual(values = c("TRUE" = 18, "FALSE" = 23),
                       labels = c("TRUE" = "p<0.05", "FALSE" = "ns"),
                       name   = NULL) +
    scale_x_log10(
      limits = c(0.35, max_x),
      breaks = c(0.5, 0.7, 0.9, 1.0, 1.11, 1.5, 2.0),
      labels = c("0.50","0.70","0.90","1.00","1.11","1.50","2.00")
    ) +
    facet_wrap(~outcome, ncol = 1) +
    labs(
      title    = title_str,
      subtitle = "Blue band = ROPE [0.90–1.11]\nThick = 90% CrI; thin = 95% CrI\nDiamond = CrI does not cross 1.0",
      x = "Hazard Ratio (log scale, first treatment vs second)",
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          strip.text      = element_text(face = "bold", size = 9))
}

p_forest <- make_forest(
  bind_rows(results_mace, results_bleed),
  "Bayesian NMA: DAPT after ACS (20 trials, n ≈ 77,714)"
)
ggsave("outputs/figures/oct/06_forest_main.pdf", p_forest,
       width = 10, height = 7)

# ── Figure 2: Posterior density plots ─────────────────────────────────────────
make_density <- function(fit, title_str) {
  d <- as_draws_df(fit)
  tibble(
    `T vs C` = exp(d$b_treatT),
    `P vs C` = exp(d$b_treatP),
    `P vs T` = exp(d$b_treatP - d$b_treatT)
  ) %>%
    pivot_longer(everything(), names_to = "comparison", values_to = "HR") %>%
    ggplot(aes(x = HR, fill = comparison, color = comparison)) +
    annotate("rect", xmin = 0.9, xmax = 1.11,
             ymin = -Inf, ymax = Inf, fill = "#4393C3", alpha = 0.08) +
    geom_density(alpha = 0.30, linewidth = 0.8) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey30") +
    geom_vline(xintercept = c(0.9, 1.11), linetype = "dotted",
               color = "#4393C3", linewidth = 0.6) +
    scale_x_continuous(limits = c(0.4, 2.0), name = "Hazard Ratio") +
    scale_fill_manual(values  = comp_colours, name = "Comparison") +
    scale_color_manual(values = comp_colours, name = "Comparison") +
    labs(
      title    = title_str,
      subtitle = "Blue band = ROPE [0.90–1.11]",
      y        = "Posterior density"
    ) +
    theme_bw(base_size = 11)
}

p_dens <- make_density(fit_mace_re,  "MACE: posterior HR distributions") /
          make_density(fit_bleed_re, "Major Bleeding: posterior HR distributions")
ggsave("outputs/figures/oct/07_posterior_densities.pdf", p_dens,
       width = 9, height = 9)

# ── Figure 3: ROPE stacked bar ────────────────────────────────────────────────
p_rope <- bind_rows(results_mace, results_bleed) %>%
  select(outcome, comparison, pr_clin_benefit, pr_rope, pr_clin_harm) %>%
  pivot_longer(
    cols = c(pr_clin_benefit, pr_rope, pr_clin_harm),
    names_to = "region", values_to = "prob"
  ) %>%
  mutate(
    region = factor(region,
      levels = c("pr_clin_harm", "pr_rope", "pr_clin_benefit"),
      labels = c("Pr(HR>1.11)\nHarm", "Pr(0.9≤HR≤1.11)\nEquivalence",
                 "Pr(HR<0.9)\nBenefit")),
    comparison = factor(comparison, levels = c("T vs C","P vs C","P vs T"))
  ) %>%
  ggplot(aes(x = comparison, y = prob, fill = region)) +
  geom_col(position = "stack", width = 0.55) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "white",
             linewidth = 0.8) +
  scale_fill_manual(values = c(
    "Pr(HR>1.11)\nHarm"          = "#D73027",
    "Pr(0.9≤HR≤1.11)\nEquivalence" = "#FEE08B",
    "Pr(HR<0.9)\nBenefit"        = "#1A9850"
  ), name = NULL) +
  scale_y_continuous(labels = percent, name = "Posterior probability") +
  facet_wrap(~outcome) +
  labs(
    title    = "ROPE Analysis: Probability of Benefit, Equivalence, or Harm",
    subtitle = "Thresholds: HR<0.90 = clinical benefit; HR>1.11 = clinical harm",
    x = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        strip.text      = element_text(face = "bold", size = 9))

ggsave("outputs/figures/oct/08_rope.pdf", p_rope, width = 10, height = 5.5)

# ── Figure 4: Sensitivity forest plot ─────────────────────────────────────────
p_sa <- sa_results %>%
  mutate(
    outcome    = factor(outcome, levels = unique(outcome)),
    comparison = factor(comparison, levels = c("T vs C","P vs C","P vs T"))
  ) %>%
  ggplot(aes(x = median_HR, y = outcome, color = comparison,
             group = comparison)) +
  annotate("rect", xmin = 0.9, xmax = 1.11,
           ymin = -Inf, ymax = Inf, fill = "#4393C3", alpha = 0.08) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_errorbarh(
    aes(xmin = lo95_HR, xmax = hi95_HR),
    height = 0.35, linewidth = 0.6,
    position = position_dodge(width = 0.6)
  ) +
  geom_point(size = 2.5, shape = 18,
             position = position_dodge(width = 0.6)) +
  scale_color_manual(values = comp_colours, name = "Comparison") +
  scale_x_log10(
    breaks = c(0.6, 0.8, 0.9, 1.0, 1.11, 1.3, 1.6),
    labels = c("0.60","0.80","0.90","1.00","1.11","1.30","1.60")
  ) +
  labs(
    title    = "Sensitivity Analysis: MACE (all models)",
    subtitle = "Blue band = ROPE; bars = 95% CrI",
    x        = "Hazard Ratio (log scale)",
    y        = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        axis.text.y     = element_text(size = 8))

ggsave("outputs/figures/oct/09_sensitivity.pdf", p_sa, width = 10, height = 8)
ggsave("outputs/figures/oct/09_sensitivity.png", p_sa, width = 10, height = 8)
# ── Figures 10 & 11: Bayesian updating diagrams ───────────────────────────────

# Shared helper: Bayesian update under normal approximation
bayes_update <- function(prior_mean, prior_sd, like_mean, like_sd) {
  prec_prior <- 1 / prior_sd^2
  prec_like  <- 1 / like_sd^2
  post_prec  <- prec_prior + prec_like
  list(
    mean = (prior_mean * prec_prior + like_mean * prec_like) / post_prec,
    sd   = sqrt(1 / post_prec)
  )
}

# Shared plot function
update_plot <- function(prior_mean, prior_sd,
                        like_mean,  like_sd,
                        prior_label, like_label, post_label,
                        x_label, title_str, subtitle_str,
                        xlim = c(0.4, 1.75)) {
  post <- bayes_update(prior_mean, prior_sd, like_mean, like_sd)
  x_seq <- seq(log(xlim[1]), log(xlim[2]), length.out = 600)

  bind_rows(
    tibble(x = x_seq, y = dnorm(x_seq, prior_mean, prior_sd), dist = prior_label),
    tibble(x = x_seq, y = dnorm(x_seq, like_mean,  like_sd),  dist = like_label),
    tibble(x = x_seq, y = dnorm(x_seq, post$mean,  post$sd),  dist = post_label)
  ) %>%
    mutate(dist = factor(dist, levels = c(prior_label, like_label, post_label))) %>%
    ggplot(aes(x = exp(x), y = y, fill = dist, color = dist)) +
    annotate("rect", xmin = 0.9, xmax = 1.11,
             ymin = -Inf, ymax = Inf, fill = "#4393C3", alpha = 0.08) +
    geom_area(alpha = 0.25, position = "identity") +
    geom_line(linewidth = 1.0) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
    scale_x_continuous(name = x_label, limits = xlim,
                       breaks = c(0.5, 0.7, 0.9, 1.0, 1.11, 1.4, 1.75)) +
    scale_fill_manual(values  = c("#377EB8", "#4DAF4A", "#E41A1C"), name = NULL) +
    scale_color_manual(values = c("#377EB8", "#4DAF4A", "#E41A1C"), name = NULL) +
    labs(title    = title_str,
         subtitle = subtitle_str,
         y        = "Density") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          legend.text     = element_text(size = 9))
}

# ── Figure 10: TC4 updating (T vs C, MACE) ────────────────────────────────────
# Prior: BNMA 2023 T vs C MACE — HR 0.95, 95%CrI 0.81–1.14 [Table 2 RE]
# Likelihood: TC4 vague-prior result — RR 0.95, 95%CrI 0.67–1.35 [CMAJ Table 2]
p_fig10 <- update_plot(
  prior_mean   = log(0.95),
  prior_sd     = (log(1.14) - log(0.81)) / (2 * 1.96),
  like_mean    = log(0.95),
  like_sd      = (log(1.35) - log(0.67)) / (2 * 1.96),
  prior_label  = "BNMA 2023 prior\n(T vs C, n=57,814)",
  like_label   = "TC4 likelihood\n(T vs C, n=1,005)",
  post_label   = "Updated posterior\n(prior \u00d7 likelihood)",
  x_label      = "Hazard Ratio: T vs C (MACE)",
  title_str    = "Figure 10. Bayesian Updating: T vs C for MACE",
  subtitle_str = paste0(
    "Prior: Kutcher 2023 BNMA (HR 0.95, 95%CrI 0.81\u20131.14)\n",
    "Likelihood: TC4 trial (Kutcher 2025 CMAJ, HR 0.95, 95%CrI 0.67\u20131.35)\n",
    "Blue band = ROPE [0.90\u20131.11]"
  )
)

ggsave("outputs/figures/oct/10_bayesian_updating_TC4.pdf", p_fig10,
       width = 8, height = 5)

# ── Figure 11: TUXEDO-2 updating (P vs T, MACE + Bleeding) ────────────────────
# Prior MACE P vs T: BNMA 2023 RE model — HR 0.87, 95%CrI 0.74–1.06
# Derived as b_treatP - b_treatT; SD from CI width
# Prior bleeding P vs T: derived from BNMA P vs C and T vs C bleeding posteriors
#   P vs C bleed HR 1.24 (0.99-1.57); T vs C bleed HR 1.16 (0.98-1.48)
#   P vs T = log(1.24) - log(1.16) on log scale

# TUXEDO-2 data:
#   MACE P vs T: 77/899 P vs 90/901 T → HR = 0.857; SE from event counts
#   Bleeding P vs T: reported HR(T vs P) 1.31 (0.88-1.95) → P vs T = 0.763

p_fig11_mace <- update_plot(
  prior_mean   = -0.1393,   # log(0.87)
  prior_sd     =  0.0917,   # from (log(1.06)-log(0.74))/(2*1.96)
  like_mean    = -0.1538,   # log(77/899 / (90/901))
  like_sd      =  0.1552,   # sqrt(1/77 + 1/90)
  prior_label  = "BNMA 2023 prior\n(P vs T MACE, n=57,814)",
  like_label   = "TUXEDO-2 likelihood\n(P vs T MACE, n=1,800)",
  post_label   = "Updated posterior\n(prior \u00d7 likelihood)",
  x_label      = "Hazard Ratio: P vs T (MACE)",
  title_str    = "Figure 11a. Bayesian Updating: P vs T for MACE",
  subtitle_str = paste0(
    "Prior: Kutcher 2023 BNMA (HR 0.87, 95%CrI 0.74\u20131.06)\n",
    "Likelihood: TUXEDO-2 (Bangalore 2026, HR 0.86, 95%CrI 0.64\u20131.15)\n",
    "Posterior: HR 0.87, 95%CrI [0.74, 1.01] | Blue band = ROPE [0.90\u20131.11]"
  )
)

p_fig11_bleed <- update_plot(
  prior_mean   =  0.0667,   # log(1.07) — P vs T bleeding from BNMA
  prior_sd     =  0.1578,   # combined uncertainty P vs C and T vs C bleeding
  like_mean    = -0.2700,   # log(0.763) — TUXEDO-2 P vs T bleeding
  like_sd      =  0.2030,   # from (log(1.95)-log(0.88))/(2*1.96)
  prior_label  = "BNMA 2023 prior\n(P vs T bleeding, n=57,110)",
  like_label   = "TUXEDO-2 likelihood\n(P vs T bleeding, n=1,800)",
  post_label   = "Updated posterior\n(prior \u00d7 likelihood)",
  x_label      = "Hazard Ratio: P vs T (Major Bleeding)",
  title_str    = "Figure 11b. Bayesian Updating: P vs T for Major Bleeding",
  subtitle_str = paste0(
    "Prior: BNMA 2023 derived (HR 1.07, 95%CrI 0.76\u20131.51)\n",
    "Likelihood: TUXEDO-2 (Bangalore 2026, HR 0.76, 95%CrI 0.51\u20131.14)\n",
    "Posterior: HR 0.94, 95%CrI [0.74, 1.20] | Blue band = ROPE [0.90\u20131.11]"
  )
)

p_fig11 <- p_fig11_mace / p_fig11_bleed +
  plot_annotation(
    title   = "Figure 11. Bayesian Updating: P vs T After TUXEDO-2",
    caption = paste0(
      "TUXEDO-2 (Bangalore 2026): diabetes + multivessel PCI, India, n=1,800.\n",
      "Prior derived from Kutcher 2023 BNMA (17 RCTs, n=57,814).\n",
      "Note: TUXEDO-2 bleeding prior derived from P vs C and T vs C BNMA ",
      "posteriors; TUXEDO-2 population (100% diabetes, 85% triple-vessel) ",
      "may not generalise to the broader ACS population."
    )
  )

ggsave("outputs/figures/oct/11_bayesian_updating_TUXEDO2.pdf", p_fig11,
       width = 9, height = 9)

# =============================================================================
# STEP 12: DISCUSSION (printed to console and saved as text)
# =============================================================================
cat("\n=== DISCUSSION ===\n")

discussion_text <- "
DISCUSSION  (October 2026 revision: SWITCH-SWEDEHEART added)
==========================================================================

Comparison of Bayesian NMA (this analysis) with Maqsood et al. 2026 JAMA
Cardiology frequentist NMA, and with our June 2026 (19-trial) analysis.

--------------------------------------------------------------------------
1. HEADLINE RESULTS COMPARED
--------------------------------------------------------------------------

Maqsood 2026 (frequentist NMA, 15 trials, n=48,904):
  MACE:     Prasugrel vs clopidogrel OR 0.80 (95%CI 0.69-0.93) - significant
            Ticagrelor vs clopidogrel OR 0.97 (95%CI 0.84-1.12) - non-significant
            Prasugrel vs ticagrelor   OR 0.83 (95%CI 0.70-0.98) - significant
  Bleeding: Ticagrelor vs clopidogrel OR 1.24 (95%CI 1.01-1.52) - significant
  Conclusion: 'Prasugrel provided the optimal balance for efficacy and safety.'

This Bayesian NMA (20 trials, n=77,714; primary cloglog RE model,
SWITCH-SWEDEHEART entered at its adjusted estimate - see section 2):
  MACE:     P vs C  HR 0.87 (95%CrI 0.77-1.01), Pr(benefit) 71%, Pr(ROPE) 29%
            T vs C  HR 0.97 (95%CrI 0.83-1.14), Pr(benefit) 16%, Pr(ROPE) 80%
            P vs T  HR 0.90 (95%CrI 0.76-1.10), Pr(benefit) 52%, Pr(ROPE) 46%
  Bleeding: T vs C  HR 1.27 (95%CrI 1.04-1.61), Pr(harm) 90%
            P vs C  HR 1.10 (95%CrI 0.89-1.34), Pr(harm) 45%
            P vs T  HR 0.86 (95%CrI 0.69-1.07), Pr(benefit) 68%, Pr(ROPE) 31%
  (Pr(benefit) = Pr(HR < 0.90); Pr(harm) = Pr(HR > 1.11);
   ROPE = HR 0.90-1.11)

--------------------------------------------------------------------------
2. SWITCH-SWEDEHEART: WHAT IT ADDS AND HOW IT WAS INCORPORATED
--------------------------------------------------------------------------

SWITCH-SWEDEHEART (Omerovic 2026 NEJM) is a registry-based, open-label,
stepped-wedge cluster-randomised trial comparing a default prasugrel
with a default ticagrelor POLICY in 17,095 unselected Swedish ACS-PCI
patients (mean age 70; 37% >=75 yr; patients on oral anticoagulation or
with prior stroke retained). Seven regions grouped into 3 clusters
switched from ticagrelor to prasugrel in random order over four 9-month
periods. It is the largest P vs T trial, and adds 17,095 patients to the
network (60,619 -> 77,714).

Raw (crude) results:
  MACE:     1113/9444 (11.8%) T policy vs 848/7651 (11.1%) P policy
            RR 0.94 (95%CI 0.86-1.02); OR 0.93 (0.85-1.03)
  Bleeding: 413/9444 (4.4%) vs 321/7651 (4.2%); OR 0.96 (0.83-1.11)
Adjusted (trial's prespecified primary analysis):
  MACE:     OR 0.90 (95%CI 0.77-1.06)
  Bleeding: OR 0.80 (95%CI 0.64-0.99)

Why the adjusted estimate is used in the primary analysis.
  The crude comparison is not a valid randomised contrast in a
  stepped-wedge trial. By design, every cluster starts on ticagrelor and
  ends on prasugrel, so policy is confounded with calendar time: any
  secular change in case-mix, PCI technique or registry coding over
  2021-2025 appears as a treatment effect in a crude analysis.
  Randomisation in this design protects only the ORDER in which clusters
  switch, and the treatment effect is identified only after adjusting
  for time. In addition, patients within a region share practice
  patterns, so outcomes are correlated within clusters and a crude
  binomial analysis that treats 17,095 patients as independent
  overstates the precision (crude MACE log-OR SE 0.048 vs 0.082 implied
  by the reported CI). The trial's primary model - a generalised linear
  mixed model with a random intercept for cluster, a restricted cubic
  spline for calendar time, prespecified baseline covariates and
  multiple imputation - addresses both problems and is the estimate the
  investigators designed the trial to report. Using crude counts would
  both bias the estimate and give this trial roughly 2.8 times its
  proper weight for MACE (1.9 times for bleeding).

How the adjustment was carried into the NMA.
  Our model is arm-based (binomial counts per arm), so the adjusted OR
  cannot be entered directly. We therefore converted it into integer
  'effective counts', a generalisation of the Cochrane Handbook
  effective-sample-size method for cluster-randomised trials
  (section 23.1.4): (i) the observed ticagrelor-arm risk (11.8% MACE,
  4.4% bleeding) was kept; (ii) the prasugrel-arm risk was set so that
  the odds ratio equals the adjusted OR; (iii) both arms were scaled by
  a common factor, preserving the 9444:7651 allocation, so that the
  binomial standard error of the log-OR equals that implied by the
  published 95% CI. This gives MACE 399/3384 (T) vs 294/2742 (P),
  reproducing OR 0.898 (0.766-1.054), and bleeding 212/4858 vs
  139/3936, reproducing OR 0.802 (0.645-0.998). The implied design
  effects are 2.8 (MACE) and 1.9 (bleeding). Because the arm risks are
  around 11% and 4%, the cloglog HRs these counts imply (0.90 and 0.81)
  are essentially the same as the ORs.

Sensitivity to this choice (SA9, crude counts):
  MACE P vs T:  HR 0.90 (0.76-1.10) adjusted -> 0.90 (0.77-1.09) crude;
                Pr(benefit) 52% -> 48%
  Bleed P vs T: HR 0.86 (0.69-1.07) adjusted -> 0.91 (0.71-1.11) crude;
                Pr(benefit) 68% -> 47%
  Bleed T vs C: HR 1.27 (1.04-1.61) adjusted -> 1.23 (1.01-1.56) crude
  MACE T vs C and P vs C were essentially unchanged.
  The choice matters little for MACE but materially for bleeding,
  where the crude counts both dilute SWITCH-SWEDEHEART's prasugrel
  advantage and over-weight it.

Caveats (informed by the Supplementary Appendix, Tables S7-S9).
  a) Robustness of the adjusted estimate. The trial's prespecified
     sensitivity analyses (Table S9) give ORs of 0.82 to 0.94:
       additional cluster-period random effect  0.90 (0.77-1.06)
       alternative covariate set (adds OAC)     0.90 (0.77-1.06)
       log link (risk ratio)                    0.91 (0.78-1.05)
       time-varying policy effect (post hoc)    0.87 (0.73-1.03)
       complete case                            0.94 (0.77-1.15)
       fixed-effects (cluster) model            0.94 (0.85-1.05)
       categorical time                         0.82 (0.70-0.96)
       GEE (population-averaged)                0.88 (0.85-0.92)
     The primary GLMM estimate sits in the middle of this range. With
     only 3 clusters, the GEE interval is unreliable (the investigators
     warn that sandwich variance estimates are too narrow) and the
     categorical-time model cannot use periods 1 and 4, where policy and
     period are fully aligned. Neither is a better choice than the
     primary model. The crude counts used in SA9 (OR 0.93) are close to
     the least favourable sensitivity estimates (0.94), so SA9 also
     serves as a lower bound on the prasugrel effect from this trial.
  b) Estimand. The adjusted OR is conditional on baseline covariates
     and, being non-collapsible, could in principle differ from the
     unadjusted (marginal) estimates of the other trials. The
     covariate-adjusted risk ratio, which is collapsible, is
     essentially the same (RR 0.91 vs OR 0.90), so this distinction
     is negligible here.
  c) Policy dilution. The estimand is a treatment POLICY effect. At
     discharge, 63% of prasugrel-policy patients on a P2Y12 inhibitor
     received prasugrel, 10% ticagrelor and 27% clopidogrel; in the
     ticagrelor-policy group, 77% received ticagrelor and 22%
     clopidogrel (Table S7). Because roughly a quarter of both groups
     received clopidogrel and a tenth of the prasugrel group received
     ticagrelor, the P vs T contrast is biased toward the null as an
     estimate of the drug-versus-drug effect.
  d) Outcomes are registry-ascertained (ICD-10) and not adjudicated,
     and the broad ICD-10 bleeding definition differs from the
     BARC/TIMI/PLATO definitions used in the other trials. The trial
     reports no sensitivity analyses for bleeding, so the adjusted
     bleeding OR (0.80) cannot be checked in the same way.

Change from the June 19-trial analysis (June -> October, primary):
  MACE P vs T:   HR 0.89 (0.73-1.12) -> 0.90 (0.76-1.10)
                 Pr(benefit) 54% -> 52%; Pr(ROPE) 43% -> 46%
  MACE P vs C:   HR 0.87 (0.76-1.03) -> 0.87 (0.77-1.01); Pr(benefit) 70% -> 71%
  MACE T vs C:   HR 0.975 (0.82-1.18) -> 0.97 (0.83-1.14); Pr(ROPE) 76% -> 80%
  Bleed P vs T:  HR 0.88 (0.67-1.14) -> 0.86 (0.69-1.07); Pr(benefit) 58% -> 68%
  Bleed T vs C:  HR 1.26 (1.01-1.63) -> 1.27 (1.04-1.61); Pr(harm) 86% -> 90%
SWITCH-SWEDEHEART narrowed the credible intervals without moving the
MACE estimates; its main contribution is to strengthen the evidence
that prasugrel causes less bleeding than ticagrelor. Excluding it (SA8)
reproduces the June results.

--------------------------------------------------------------------------
3. METHODOLOGICAL ADVANTAGES OF THE BAYESIAN APPROACH
--------------------------------------------------------------------------

A. Uncertainty is quantified rather than dichotomised.
   The frequentist approach forces a binary significant/non-significant
   conclusion. Maqsood 2026 reports prasugrel vs clopidogrel as
   'significant' (OR 0.80) and ticagrelor vs clopidogrel as
   'non-significant' (OR 0.97), inviting the reader to treat these as
   qualitatively different findings. The posterior distributions make
   explicit that both estimates carry substantial uncertainty: the 71%
   probability of clinically meaningful MACE benefit with prasugrel is
   notable but far from certain, and the 16% probability with
   ticagrelor, while lower, is not zero. Clinical decisions involve
   weighing probabilities, not accepting binary verdicts.

B. The ROPE framework acknowledges practical equivalence.
   For MACE, 80% of the posterior for T vs C and 46% of the posterior
   for P vs T lie within the ROPE [0.90, 1.11]. These probability
   statements are more informative than declaring 'prasugrel is
   superior' on the basis of p < 0.05.

C. Between-study heterogeneity is modelled explicitly.
   The between-study SD for T vs C is tau = 0.14 (95%CrI 0.02-0.37) for
   MACE and 0.19 (0.04-0.43) for bleeding, with Pr(tau > 0.10) of 72%
   and 85%. Heterogeneity of this size is clinically meaningful and is
   consistent with the divergence between the PLATO global result
   (HR 0.84), the PLATO North American subgroup (HR 1.25) and TC4
   (RR 0.95). Ignoring it matters: a fixed-effect model gives T vs C
   HR 0.92 (0.87-0.98), P vs C 0.85 (0.80-0.90) and P vs T 0.92
   (0.86-0.99), all spuriously precise compared with the random-effects
   estimates above.

D. The time offset is correctly modelled.
   The cloglog link with a log(time) offset places all trials on a
   common hazard scale, accounting for TRILOGY ACS (30 months), TRITON
   (15 months) and the 6-month trials. The JAMA 2026 analysis uses odds
   ratios that ignore variation in follow-up duration. (SA1, a logit
   model without offset, gives similar results: P vs C 0.85, P vs T
   0.87.)

E. Prior information is made explicit rather than hidden.
   The weakly informative priors (normal(0, 0.5) on log-HR) encode only
   that extreme effects (HR < 0.37 or > 2.7) are implausible. The TC4
   trial (Kutcher 2025) used our 2023 BNMA as an explicit summary
   prior, and the present analysis is itself a sequential update:
   2023 BNMA -> June 2026 (+TC4, TUXEDO-2) -> October 2026
   (+SWITCH-SWEDEHEART), each step shifting the conclusions only
   modestly.

F. Robustness.
   Across SA1-SA9 the posterior medians are stable: T vs C 0.96-0.97,
   P vs C 0.83-0.88, P vs T 0.87-0.90. No single trial or modelling
   choice drives the MACE results. The largest shift comes from SA2
   (PCI-only), where P vs C Pr(benefit) rises to 84%. For bleeding, the
   handling of SWITCH-SWEDEHEART (SA9) is the most influential choice.

G. Network consistency (P vs T).
   Direct evidence (5 P vs T trials) gives HR 0.82 (95%CrI 0.66-1.03);
   the indirect estimate via clopidogrel (P vs C and T vs C trials only)
   gives HR 0.97 (0.73-1.34). The ratio of direct to indirect HRs is
   approximately 0.85 (0.58-1.24): the direct trials favour prasugrel
   somewhat more than the indirect route, but the difference is well
   within uncertainty and gives no evidence of important inconsistency.

--------------------------------------------------------------------------
4. WHAT THE DATA SHOW AFTER ~78,000 PATIENTS
--------------------------------------------------------------------------

The best summary of 20 RCTs and 77,714 patients is one of modest,
uncertain and context-dependent treatment differences, not the clear
superiority of prasugrel claimed by Maqsood 2026.

Ticagrelor vs clopidogrel (MACE):
  The posterior is centred just below 1.0 (HR 0.97) with 80% of its mass
  in the ROPE, a 16% probability of clinically meaningful benefit and a
  5% probability of harm. Ticagrelor and clopidogrel produce similar
  MACE outcomes in the average trial population. The North American
  evidence (PLATO NA subgroup, TC4) shows no benefit, raising the
  possibility that the pooled result reflects geographic heterogeneity
  that may not apply to North American practice.

Ticagrelor vs clopidogrel (bleeding):
  This remains the clearest signal in the dataset: a 90% posterior
  probability of clinically meaningful bleeding harm with ticagrelor
  (HR 1.27, 95%CrI 1.04-1.61). The intracranial haemorrhage signal
  reported by Maqsood 2026 (OR 1.89, 95%CI 1.08-3.33) is consistent.
  Whatever the uncertain MACE benefit of ticagrelor, it comes at a real
  cost in bleeding.

Prasugrel vs clopidogrel (MACE):
  A 71% probability of clinically meaningful MACE benefit
  (HR 0.87, 95%CrI 0.77-1.01) is the strongest efficacy signal in the
  network. The direct evidence comes mainly from TRITON-TIMI 38, with
  further indirect support from the P vs T trials. Nevertheless 29% of
  the posterior lies in the ROPE, so practical equivalence with
  clopidogrel cannot be excluded, and the 45% probability of clinically
  meaningful excess bleeding with prasugrel is a non-trivial safety
  concern.

Prasugrel vs ticagrelor (net clinical balance):
  Both outcomes lean toward prasugrel: a 52% probability of clinically
  meaningful MACE benefit and a 68% probability of clinically
  meaningful bleeding reduction, with at most a 2% probability of
  clinically meaningful harm on either. SWITCH-SWEDEHEART, the largest
  P vs T trial, left the MACE comparison essentially unchanged and
  strengthened the bleeding advantage. 'Probably no worse and possibly
  modestly better' is the appropriate characterisation; 'definitively
  superior' is not.

--------------------------------------------------------------------------
5. OVERALL CONCLUSION
--------------------------------------------------------------------------

After 20 RCTs and 77,714 randomised patients, including the large
pragmatic SWITCH-SWEDEHEART trial, the Bayesian summary is that there
is little to choose between the three P2Y12 inhibitors in terms of a
definitive treatment hierarchy.

Ticagrelor is the weakest performer: modest and uncertain MACE efficacy
relative to clopidogrel combined with a clear bleeding penalty. It is
difficult to recommend over clopidogrel outside specific circumstances,
and it is probably no better than prasugrel for efficacy and probably
worse for bleeding.

Prasugrel has the strongest efficacy signal but a possible bleeding
cost relative to clopidogrel. In patients at low bleeding risk
(younger, no prior stroke, not on anticoagulation), prasugrel likely
provides the best net clinical benefit. Its MACE advantage over
ticagrelor is smaller and less certain than ISAR-REACT 5 alone
suggested, but SWITCH-SWEDEHEART adds evidence that it causes less
bleeding.

Clopidogrel performs at least as well as ticagrelor on MACE with less
bleeding, and its efficacy gap with prasugrel remains uncertain (the
95%CrI for P vs C MACE crosses 1.0). In high bleeding-risk and North
American patients in particular, where ticagrelor shows no consistent
benefit, clopidogrel remains a reasonable first-line choice that
current guidelines undervalue.

The claim by Maqsood 2026 that 'prasugrel provided the optimal balance
for efficacy and safety' is not wrong, but it overstates the certainty
of the evidence. A more accurate summary: prasugrel is probably the
most efficacious option for post-PCI patients without bleeding risk
factors, ticagrelor offers uncertain and possibly geographically
variable efficacy with a real bleeding cost, and clopidogrel remains a
reasonable first-line choice, especially in older patients, those at
high bleeding risk, and North American populations where the PLATO
results have not been replicated.
==========================================================================
"

cat(discussion_text)
writeLines(discussion_text, "outputs/tables/oct/05_discussion.txt")
cat("Discussion saved to outputs/tables/oct/05_discussion.txt\n")

# =============================================================================
# STEP 12: TABLES  (gt package — publication-quality PDF output)
# =============================================================================
cat("\n=== STEP 12: TABLES ===\n")

# gt PDF export requires webshot2 + Chrome/Chromium.
# Install once with: install.packages("webshot2"); webshot2::install_chromium()
# Fallback: all tables also saved as .html (open in browser, print to PDF)
# Raw data always saved as .csv for sharing regardless of rendering.

# ── Shared gt theme ───────────────────────────────────────────────────────────
gt_theme_dapt <- function(gt_obj) {
  gt_obj %>%
    tab_options(
      table.font.size                    = px(11),
      table.width                        = pct(100),
      column_labels.font.weight          = "bold",
      column_labels.background.color     = "#2C3E50",
      row_group.font.weight              = "bold",
      row_group.background.color         = "#ECF0F1",
      stub.font.weight                   = "bold",
      heading.title.font.size            = px(13),
      heading.subtitle.font.size         = px(10),
      table.border.top.color             = "#2C3E50",
      table.border.bottom.color          = "#2C3E50",
      column_labels.border.bottom.color  = "#2C3E50",
      data_row.padding                   = px(4)
    ) %>%
    # column_labels.font.color removed from tab_options in newer gt versions
    # set white text on dark header via tab_style instead
    tab_style(
      style     = cell_text(color = "white", weight = "bold"),
      locations = cells_column_labels()
    ) %>%
    tab_style(
      style     = cell_text(color = "white", weight = "bold"),
      locations = cells_row_groups()
    )
}

# ── Save helper: html always; pdf via webshot2 if available ───────────────────
save_gt <- function(gt_obj, stem) {
  html_path <- paste0("outputs/tables/oct/", stem, ".html")
  pdf_path  <- paste0("outputs/tables/oct/", stem, ".pdf")
  gtsave(gt_obj, html_path)
  tryCatch(
    gtsave(gt_obj, pdf_path),
    error = function(e) {
      message("PDF save failed for ", stem,
              " — webshot2/Chrome may not be installed.\n",
              "HTML saved at: ", html_path, "\n",
              "Install with: install.packages('webshot2')")
    }
  )
  invisible(gt_obj)
}

# ── Table 1: Trial reconciliation ─────────────────────────────────────────────
recon <- tribble(
  ~Trial,           ~Comparison, ~N,     ~Followup, ~Kutcher23,  ~JAMA26,     ~Main,       ~PCIonly, ~Note,
  "TREAT",           "T vs C",  3799,   "12m",     "Yes",       "No",        "Yes",       "No",     "Non-PCI (fibrinolysis)",
  "Gasecka",         "T vs C",  55,     "6m",      "Yes",       "No",        "Yes",       "Yes",    "0 MACE events",
  "POPular AGE",     "T vs C",  1002,   "12m",     "Yes",       "Yes",       "Yes",       "Yes",    "",
  "PHILO",           "T vs C",  801,    "12m",     "Yes",       "Yes",       "Yes",       "Yes",    "",
  "He",              "T vs C",  266,    "6m",      "Yes",       "No",        "Yes",       "Yes",    "",
  "Mohareb",         "T vs C",  943,    "12m",     "Yes",       "No",        "Yes",       "Yes",    "",
  "TICAKOREA",       "T vs C",  800,    "12m",     "Yes",       "Yes",       "Yes",       "Yes",    "",
  "PLATO",           "T vs C",  18624,  "12m",     "Yes",       "Yes",       "Yes",       "Yes",    "NA subgroup: no benefit",
  "Wu",              "T vs C",  350,    "12m",     "Yes",       "No",        "Yes",       "Yes",    "0 MACE events",
  "TC4",             "T vs C",  1005,   "12m",     "No",        "No",        "Yes \u2605","Yes",    "North American; cluster RCT; no COI",
  "TRILOGY ACS",     "P vs C",  9326,   "30m",     "Yes",       "No",        "Yes",       "No",     "0% PCI; excluded PCI-only SA",
  "PRASFIT-ACS",     "P vs C",  1363,   "12m",     "Yes",       "Yes",       "Yes",       "Yes",    "",
  "Elderly ACS II",  "P vs C",  1443,   "12m",     "Yes",       "Yes",       "Yes",       "Yes",    "",
  "TRITON-TIMI 38",  "P vs C",  13608,  "15m",     "Yes",       "Yes",       "Yes",       "Yes",    "",
  "Yabe",            "P vs C",  76,     "12m",     "Yes",       "No",        "Yes",       "Yes",    "",
  "PRAGUE-18",       "P vs T",  1230,   "12m",     "Yes (1yr)", "Yes (1mo)", "Yes (1yr)", "Yes",    "39% crossover to C \u2014 SA5",
  "ISAR-REACT 5",    "P vs T",  4018,   "12m",     "Yes",       "Yes",       "Yes",       "Yes",    "",
  "REDUCE-MVI",      "P vs T",  110,    "12m",     "Yes",       "No",        "Yes",       "Yes",    "Small STEMI trial",
  "TUXEDO-2",        "P vs T",  1800,   "12m",     "No",        "Yes",       "Yes \u2605","Yes",    "DM+MVD; SA4",
  "SWITCH-SWEDEHEART","P vs T", 17095,  "12m",     "No",        "No",        "Yes \u2605","Yes",    "Stepped-wedge cluster RCT; adjusted OR as effective counts; SA8, SA9",
  "Tang 2016",       "T vs C",  400,    "6m",      "No (RoB)",  "Yes",       "No",        "No",     "Excluded: high RoB",
  "Wang 2019",       "T vs C",  298,    "6m",      "No (RoB)",  "Yes",       "No",        "No",     "Excluded: high RoB",
  "ALPHEUS",         "T vs C",  1883,   "1m",      "No (<6m)",  "Yes",       "No",        "No",     "Excluded: 1-month follow-up"
)

write_csv(recon, "outputs/tables/oct/01_reconciliation.csv")

gt_recon <- recon[, c(-9)] %>%
  mutate(N = formatC(N, format = "d", big.mark = ",")) %>%
  gt(groupname_col = "Comparison") %>%
  tab_header(
    title    = "Table 1. Trial Reconciliation",
    subtitle = paste0("Comparison of included trials across Kutcher 2023 BNMA, ",
                      "Maqsood 2026 JAMA NMA, and current analysis")
  ) %>%
  cols_label(
    Trial     = "Trial",
    N         = "N",
    Followup  = "Follow-up",
    Kutcher23 = "Kutcher 2023",
    JAMA26    = "Maqsood 2026",
    Main      = "This analysis",
    PCIonly   = "PCI-only SA"
 #   Note      = "Notes"
  ) %>%
  cols_width(
    Trial     ~ px(110),
    N         ~ px(65),
    Followup  ~ px(60),
    Kutcher23 ~ px(75),
    JAMA26    ~ px(75),
    Main      ~ px(80),
    PCIonly   ~ px(65)
  #  Note      ~ px(200)
  ) %>%


  tab_footnote(
    footnote  = "\u2605 New trials added to this analysis (not in Kutcher 2023).",
    locations = cells_column_labels(columns = Main)
  ) %>%
  gt_theme_dapt()

save_gt(gt_recon, "01_reconciliation")

# ── Table 2: Main results ──────────────────────────────────────────────────────
write_csv(
  bind_rows(results_mace, results_bleed),
  "outputs/tables/oct/02_main_results_raw.csv"
)

main_gt_data <- bind_rows(results_mace, results_bleed) %>%
  mutate(
    Outcome     = outcome,
    Comparison  = comparison,
    hr_str      = sprintf("%.2f (%.2f\u2013%.2f)", median_HR, lo95_HR, hi95_HR),
    benefit_str = sprintf("%.1f%%", pr_clin_benefit * 100),
    rope_str    = sprintf("%.1f%%", pr_rope          * 100),
    harm_str    = sprintf("%.1f%%", pr_clin_harm     * 100)
  )

gt_main <- main_gt_data %>%
  select(Outcome, Comparison, hr_str, benefit_str, rope_str, harm_str) %>%
  gt(groupname_col = "Outcome") %>%
  tab_header(
    title    = "Table 2. Main Results: Bayesian NMA",
    subtitle = paste0("20 trials, n \u2248 77,714 | ",
                      "Primary model: binomial cloglog, log(time) offset, ",
                      "random effects")
  ) %>%
  cols_label(
    Comparison  = "Comparison",
    hr_str      = "HR (95% CrI)",
    benefit_str = "Pr(HR < 0.90) Benefit",
    rope_str    = "Pr(0.90\u20131.11) Equivalence",
    harm_str    = "Pr(HR > 1.11) Harm"
  ) %>%
  cols_width(
    Comparison  ~ px(90),
    hr_str      ~ px(135),
    benefit_str ~ px(115),
    rope_str    ~ px(125),
    harm_str    ~ px(105)
  ) %>%
  cols_align(align = "center",
             columns = c(hr_str, benefit_str, rope_str, harm_str)) %>%
  tab_footnote(
    footnote  = paste0("HR: hazard ratio from cloglog model with log(time) ",
                       "offset. CrI: Bayesian credible interval. ",
                       "ROPE: region of practical equivalence [0.90, 1.11]."),
    locations = cells_column_labels(columns = hr_str)
  ) %>%
  tab_footnote(
    footnote  = paste0("Pr(HR < 0.90): probability of clinically meaningful ",
                       "benefit (>10% relative reduction). ",
                       "Pr(HR > 1.11): probability of clinically meaningful ",
                       "harm (>11% relative increase). ",
                       "Thresholds match Kutcher 2023 and TC4 (Kutcher 2025)."),
    locations = cells_column_labels(columns = benefit_str)
  ) %>%
  gt_theme_dapt()

save_gt(gt_main, "02_main_results")
gtsave(gt_main, "outputs/tables/oct/02_main_results.png")

# ── Table 3: Sensitivity analysis ─────────────────────────────────────────────
write_csv(sa_results, "outputs/tables/oct/03_sensitivity_raw.csv")

sa_gt_data <- sa_results %>%
  mutate(
    Model       = outcome,
    Comparison  = factor(comparison,
                         levels = c("T vs C", "P vs C", "P vs T")),
    hr_str      = sprintf("%.2f (%.2f\u2013%.2f)", median_HR, lo95_HR, hi95_HR),
    benefit_str = sprintf("%.1f%%", pr_clin_benefit * 100),
    rope_str    = sprintf("%.1f%%", pr_rope          * 100),
    harm_str    = sprintf("%.1f%%", pr_clin_harm     * 100)
  )

gt_sa <- sa_gt_data %>%
  select(Model, Comparison, hr_str, benefit_str, rope_str, harm_str,
         pr_clin_benefit, pr_clin_harm) %>%
  gt(groupname_col = "Model") %>%
  tab_header(
    title    = "Table 3. Sensitivity Analyses: MACE Outcomes",
    subtitle = "All models use binomial cloglog link with random effects unless noted"
  ) %>%
  cols_label(
    Comparison  = "Comparison",
    hr_str      = "HR (95% CrI)",
    benefit_str = "Pr(benefit)",
    rope_str    = "Pr(equivalence)",
    harm_str    = "Pr(harm)"
  ) %>%
  cols_width(
    Comparison  ~ px(90),
    hr_str      ~ px(135),
    benefit_str ~ px(95),
    rope_str    ~ px(105),
    harm_str    ~ px(90)
  ) %>%
  cols_align(align = "center",
             columns = c(hr_str, benefit_str, rope_str, harm_str)) %>%
  cols_hide(columns = c(pr_clin_benefit, pr_clin_harm)) %>%
  tab_footnote(
    footnote  = paste0(
      "SA1: logistic (logOR) model, no time offset. ",
      "SA2: PCI-only trials (excludes TREAT and TRILOGY ACS). ",
      "SA3: excludes TC4 (cluster-RCT design). ",
      "SA4: excludes TUXEDO-2 (selected DM+MVD population). ",
      "SA5: excludes PRAGUE-18 (39% crossover to clopidogrel at 1 year). ",
      "SA6: tighter between-study SD prior [half-Normal(0, 0.25)]. ",
      "SA7: wider between-study SD prior [half-Normal(0, 1.0)]. ",
      "SA8: excludes SWITCH-SWEDEHEART (cluster stepped wedge RCT design). ",
      "SA9: SWITCH-SWEDEHEART entered as crude arm-level counts instead of ",
      "effective counts reproducing its cluster- and calendar-time-adjusted OR."
    ),
    locations = cells_title(groups = "title")
  ) %>%
  gt_theme_dapt()

save_gt(gt_sa, "03_sensitivity")
gtsave(gt_sa, "outputs/tables/oct/03_sensitivity.png")

# ── Table 4: Between-study heterogeneity ──────────────────────────────────────
tau_tbl <- map_dfr(
  list(fit_mace_re, fit_bleed_re),
  ~{
    d        <- as_draws_df(.x)
    tau_cols <- grep("^sd_study__treat", names(d), value = TRUE)
    map_dfr(tau_cols, function(col) {
      vals <- d[[col]]
      tibble(
        Parameter = col %>%
                      gsub("sd_study__treat", "", .) %>%
                      gsub("Intercept", "Baseline (C)", .),
        Median    = round(median(vals), 3),
        lo95      = round(quantile(vals, 0.025), 3),
        hi95      = round(quantile(vals, 0.975), 3),
        pr_gt_01  = round(mean(vals > 0.1), 3),
        pr_gt_03  = round(mean(vals > 0.3), 3)
      )
    })
  },
  .id = "Outcome"
) %>%
  mutate(
    Outcome   = ifelse(Outcome == "1", "MACE", "Major Bleeding"),
    tau_str   = sprintf("%.3f (%.3f\u2013%.3f)", Median, lo95, hi95),
    gt_01_str = sprintf("%.1f%%", pr_gt_01 * 100),
    gt_03_str = sprintf("%.1f%%", pr_gt_03 * 100)
  )

write_csv(tau_tbl, "outputs/tables/oct/04_heterogeneity.csv")

gt_tau <- tau_tbl %>%
  select(Outcome, Parameter, tau_str, gt_01_str, gt_03_str) %>%
  gt(groupname_col = "Outcome") %>%
  tab_header(
    title    = "Table 4. Between-Study Heterogeneity (\u03c4)",
    subtitle = "Posterior distribution of between-study SD on log-hazard ratio scale"
  ) %>%
  cols_label(
    Parameter = "Parameter",
    tau_str   = "\u03c4 median (95% CrI)",
    gt_01_str = "Pr(\u03c4 > 0.10)",
    gt_03_str = "Pr(\u03c4 > 0.30)"
  ) %>%
  cols_width(
    Parameter ~ px(145),
    tau_str   ~ px(165),
    gt_01_str ~ px(105),
    gt_03_str ~ px(105)
  ) %>%
  cols_align(align = "center",
             columns = c(tau_str, gt_01_str, gt_03_str)) %>%
  tab_footnote(
    footnote  = paste0(
      "\u03c4: between-study SD on log-HR scale. ",
      "Pr(\u03c4 > 0.10): posterior probability of meaningful heterogeneity. ",
      "Pr(\u03c4 > 0.30): posterior probability of substantial heterogeneity. ",
      "Prior: half-Normal(0, 0.5) on SD class (equivalent to normal(0, 0.5) ",
      "with lower bound 0 in brms)."
    ),
    locations = cells_column_labels(columns = tau_str)
  ) %>%
  gt_theme_dapt()

save_gt(gt_tau, "04_heterogeneity")

cat("Tables saved to outputs/tables/ as .html and .pdf\n")
cat("If PDF export failed, open .html in browser and print to PDF\n")

# =============================================================================
# STEP 12c: EXPORTS FOR THE MANUSCRIPT (CJC/oct/manuscript_oct2026.qmd)
# =============================================================================
# Trial-level input data (crude counts as reported; SWITCH-SWEDEHEART primary
# analysis uses the effective counts saved separately below)
write_csv(mace_raw_crude,  "outputs/tables/oct/00_input_mace.csv")
write_csv(bleed_raw_crude, "outputs/tables/oct/00_input_bleed.csv")
write_csv(
  bind_rows(
    as_tibble(ss_mace_adj)  %>% mutate(outcome = "MACE",
                                       or = 0.90, lo = 0.77, hi = 1.06),
    as_tibble(ss_bleed_adj) %>% mutate(outcome = "Bleeding",
                                       or = 0.80, lo = 0.64, hi = 0.99)
  ),
  "outputs/tables/oct/00_switch_effective_counts.csv"
)
# Fixed-effect MACE results (contrast with the random-effects primary model)
write_csv(extract_posterior(fit_mace_fe, "MACE (cloglog FE)"),
          "outputs/tables/oct/02b_fe_results.csv")
# Convergence summary for the main models
conv_summary <- map_dfr(
  list("MACE RE" = fit_mace_re, "MACE FE" = fit_mace_fe,
       "Bleeding RE" = fit_bleed_re),
  ~ tibble(rhat_max = max(rhat(.x), na.rm = TRUE),
           ess_min  = min(ess_bulk(.x), na.rm = TRUE),
           divergences = sum(nuts_params(.x)$Value[
             nuts_params(.x)$Parameter == "divergent__"])),
  .id = "model"
)
write_csv(conv_summary, "outputs/tables/oct/08_convergence.csv")

# =============================================================================
# STEP 13: SESSION INFO
# =============================================================================
cat("\n=== ANALYSIS COMPLETE ===\n")
writeLines(
  c(
    paste("Date:", Sys.time()),
    paste("R:", R.version.string),
    paste("brms:", packageVersion("brms")),
    paste("cmdstanr:", packageVersion("cmdstanr")),
    paste("CmdStan:", cmdstanr::cmdstan_version()),
    paste("Seed: 20260529"),
    paste("Trials (main):", n_distinct(mace_data$study)),
    paste("Patients (MACE):", sum(mace_raw_crude$n))
  ),
  "outputs/tables/oct/00_session.txt"
)

cat("\nKey results (MACE primary model):\n")
print(results_mace %>%
  mutate(across(c(median_HR, lo95_HR, hi95_HR), ~round(., 3)),
         across(starts_with("pr_"), ~sprintf("%.1f%%", . * 100))) %>%
  select(comparison, median_HR, lo95_HR, hi95_HR,
         pr_clin_benefit, pr_rope, pr_clin_harm))

cat("\nAll outputs written to outputs/figures/ and outputs/tables/\n")
cat("Models cached in outputs/models/ (rerun will load from cache)\n")

