# =============================================================================
# RoB 2 traffic-light plot for the Supplementary Material
# 20 trials included in the NMA + 2 trials excluded for risk of bias
# Run from CJC/oct:  source("rob2_traffic_light.R")
# Requires: install.packages("robvis")
# =============================================================================
library(robvis)
library(ggplot2)

# RoB 2 domains:
#   D1 randomization process
#   D2 deviations from intended interventions
#   D3 missing outcome data
#   D4 measurement of the outcome
#   D5 selection of the reported result
# Values must be "Low", "Some concerns" or "High".

L <- "Low"; S <- "Some concerns"; H <- "High"

included <- data.frame(
  Study = c(
    # Ticagrelor vs clopidogrel
    "Wallentin 2009 (PLATO)",
    "Goto 2015 (PHILO)",
    "Berwanger 2019 (TREAT)",
    "Park 2019 (TICAKOREA)",
    "Gasecka 2020",
    "Gimbel 2020 (POPular AGE)",
    "Mohareb 2020",
    "He 2021",
    "Wu 2021",
    "Kutcher 2025 (TC4)",
    # Prasugrel vs clopidogrel
    "Wiviott 2007 (TRITON-TIMI 38)",
    "Roe 2012 (TRILOGY ACS)",
    "Saito 2014 (PRASFIT-ACS)",
    "Savonitto 2018 (Elderly ACS II)",
    "Yabe 2022",
    # Prasugrel vs ticagrelor
    "Motovska 2018 (PRAGUE-18)",
    "Schupke 2019 (ISAR-REACT 5)",
    "van der Hoeven 2020 (REDUCE-MVI)",
    "Bangalore 2026 (TUXEDO-2)",
    "Omerovic 2026 (SWITCH-SWEDEHEART)"
  ),
  D1 = L, D2 = L, D3 = L, D4 = L, D5 = L, Overall = L,
  stringsAsFactors = FALSE
)

# SWITCH-SWEDEHEART: stepped-wedge cluster RCT with only three randomized
# clusters (randomized switching order); policy confounded with calendar time
# by design -> some concerns for the randomization process (D1). No concerns in
# the other domains. Retained in the primary analysis (see footnote / Methods).
# TC4: cluster (time-period) randomization with a randomly generated cluster
# sequence -> low risk in all domains.
sw <- included$Study == "Omerovic 2026 (SWITCH-SWEDEHEART)"
included$D1[sw]      <- S
included$Overall[sw] <- S

# Trials excluded from the primary analysis for risk of bias.
# Judgments based on the published reports:
#
# Tang 2016 (J Cardiovasc Pharmacol 68:115-20), STEMI PPCI, n = 420 randomized
#   D1 Some concerns: "envelope method"; sequence generation and concealment
#      (sealed/opaque/sequential) not described; no baseline imbalance reported.
#   D2 High: 20 randomized patients (10 per arm) excluded after randomization
#      because of angiographic findings (stent thrombosis, multivessel disease
#      requiring CABG, no culprit lesion) and not followed; analysis was not by
#      intention to treat. Blinding of patients/carers unclear.
#   D3 Some concerns: outcomes missing for 21/420 (5%) randomized patients,
#      relevant against only 19 MACCE events in total.
#   D4 Some concerns: no blinded or adjudicated outcome assessment reported;
#      composite includes unplanned revascularization.
#   D5 High: no trial registration or protocol; several "primary" end points
#      listed, with significance reported only for selected composites.
#   Overall High.
#
# Wang 2019 (Coron Artery Dis 30:317-22), STEMI PPCI, n = 298 randomized
#   D1 Some concerns: "assigned randomly"; sequence generation and allocation
#      concealment not described; post-angiography exclusions before
#      randomization; baseline imbalances in procedural variables (thrombus
#      aspiration 38% vs 48%, direct stenting 12% vs 20%, single-vessel disease
#      47% vs 37%).
#   D2 Some concerns: described as "double-blind" but no placebo/double-dummy
#      or other blinding procedure described; adherence and analysis population
#      not reported.
#   D3 Low: outcomes reported for all 298 randomized patients.
#   D4 Some concerns: imaging/ECG readers blinded, but no blinded assessment or
#      adjudication of clinical events; MACE includes urgent target vessel
#      revascularization.
#   D5 High: no trial registration or protocol; surrogate primary end point
#      (ST resolution); clinical outcomes reported at selected time points
#      (1 and 6 months) as a non-standard composite (death, MI, TVR; stroke
#      not reported).
#   Overall High (one domain High plus some concerns in three others).
excluded <- data.frame(
  Study   = c("Tang 2016 [excluded]", "Wang 2019 [excluded]"),
  D1      = c(S, S),
  D2      = c(H, S),
  D3      = c(S, L),
  D4      = c(S, S),
  D5      = c(H, H),
  Overall = c(H, H),
  stringsAsFactors = FALSE
)

rob_df <- rbind(included, excluded)
write.csv(rob_df, "rob2_judgments.csv", row.names = FALSE)

# Traffic-light plot
rob_footnote <- paste(strwrap(paste(
  "SWITCH-SWEDEHEART: randomization concerns were noted because of cluster",
  "randomization (three clusters), but given no other concerns, its contemporary",
  "nature and the sensitivity analysis results (SA8, SA9), the trial was retained",
  "in the primary analysis. TC4 also used cluster randomization, but the cluster",
  "sequence was random and it was considered at low risk of bias. Tang 2016 and",
  "Wang 2019 were excluded for high risk of bias."), width = 110), collapse = "\n")

p_tl <- rob_traffic_light(rob_df, tool = "ROB2", psize = 8) +
  patchwork::plot_annotation(
    caption = rob_footnote,
    theme   = ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0, size = 9))
  )

# Summary (weighted bar) plot of the included trials only
p_sum <- rob_summary(included, tool = "ROB2", overall = TRUE, weighted = FALSE)

dir.create("figures", showWarnings = FALSE)
ggsave("figures/RoB2_traffic_light.png", p_tl, width = 8, height = 11.8, dpi = 300, bg = "white")
ggsave("figures/RoB2_traffic_light.pdf", p_tl, width = 8, height = 11.8)
ggsave("figures/RoB2_summary_included.png", p_sum, width = 8, height = 2.6, dpi = 300, bg = "white")

p_tl
