# Bayesian-NMA

Code and data for:

> Kutcher SA, Dendukuri N, Dandona S, Brophy JM. **Revisiting P2Y12 Inhibitors After Acute Coronary Syndromes - Additional Insights from a Bayesian Network Meta-Analysis.** *Can J Cardiol* 2026 (accepted, CJC-D-26-02373R1).

This is a Bayesian network meta-analysis of 20 randomized trials (77,714 patients). The trials compare clopidogrel, ticagrelor and prasugrel, each given with aspirin, after acute coronary syndrome. It updates our 2023 analysis ([medRxiv 2023.08.12.23294021](https://doi.org/10.1101/2023.08.12.23294021)) with three trials: TC4, TUXEDO-2 and SWITCH-SWEDEHEART.

Every number, table and figure in the manuscript and Supplementary Material can be regenerated from this repository.

## Repository layout

```
dapt_bayesian_nma_oct2026.R        Main analysis (trial data, models, diagnostics, sensitivity analyses, exports)
outputs/
  tables/oct/                      CSV results written by the main script and read by the manuscript
  figures/oct/                     Diagnostic and supplementary figures (prior/posterior predictive checks, traceplots, ...)
  models/oct/                      Fitted brms models (not tracked; created when the script is run)
CJC/oct/
  manuscript_oct2026.qmd           Manuscript (Quarto); all numbers are computed from outputs/tables/oct
  supplementary_material_oct2026.qmd  Supplementary Material (Supplemental Appendices S1-S6)
  central_illustration.R           Central Illustration
  prisma_flow_update.R             PRISMA 2020 flow diagram (Supplemental Figure S1)
  rob2_traffic_light.R             Cochrane RoB 2 judgments and traffic-light plot (Supplemental Figure S2)
  search_*.R / search_*.csv / .md  Updated literature search: strategies, de-duplication and screening log
  figures/                         Figures used by the manuscript and supplement
  references_oct.bib, vancouver.csl, reference_cjc*.docx   Bibliography, citation style and Word templates
archive/dapt_bayesian_nma_v1.R     Earlier version of the analysis (19 trials, before SWITCH-SWEDEHEART)
```

## Reproducing the results

The pipeline has two stages.

**1. Run the analysis** from the repository root. This is the slow step: it fits all the Stan models and rewrites everything in `outputs/`.

```r
source("dapt_bayesian_nma_oct2026.R")
```

**2. Build the figures and documents.** These steps read only the CSV files in `outputs/tables/oct/`, which are included in the repository. They work without step 1 and need no Stan installation.

```r
source("CJC/oct/central_illustration.R")   # Central Illustration
source("CJC/oct/prisma_flow_update.R")     # PRISMA flow diagram
setwd("CJC/oct"); source("rob2_traffic_light.R"); setwd("../..")   # RoB 2 figure
```

```bash
quarto render CJC/oct/manuscript_oct2026.qmd               # html, pdf and docx
quarto render CJC/oct/supplementary_material_oct2026.qmd
```

Figures 1 and 2 are drawn inside `manuscript_oct2026.qmd`. Supplemental Figure S3 (prior predictive check) is written by the main script to `outputs/figures/oct/`.

### Where the results come from

| Manuscript item | Source |
|---|---|
| Trial-level data (Tables 2-3) | Hard-coded in Step 1 of `dapt_bayesian_nma_oct2026.R`; exported to `outputs/tables/oct/00_input_*.csv` |
| SWITCH-SWEDEHEART effective counts (Supplemental Table S1) | Section 1c-bis; `00_switch_effective_counts.csv` |
| Main results, Table 4, Figure 2, Central Illustration | Steps 5 and 9; `02_main_results_raw.csv` |
| Fixed-effect model and LOO-CV (Supplemental Table S5) | Steps 5 and 8; `02b_fe_results.csv`, `06_loo.csv` |
| Sensitivity analyses SA1-SA9 (Supplemental Tables S2-S3, Figure S4) | Step 10; `03_sensitivity_raw.csv`, `03b_sensitivity_bleed_sa9.csv` |
| Heterogeneity (Supplemental Table S4) | `04_heterogeneity.csv` |
| Node-splitting, prasugrel vs ticagrelor (Supplemental Table S6) | Step 10b; `07_nodesplit_PT.csv` |
| Convergence diagnostics | Step 6; `08_convergence.csv` |
| Search and screening counts (Supplemental Figure S1) | `CJC/oct/search_update_all_sources_log.csv` |
| Risk-of-bias judgments (Supplemental Figure S2) | `CJC/oct/rob2_judgments.csv` |

## Software

The published results were produced with R 4.6.1, brms 2.23.0, cmdstanr 0.9.0 and CmdStan 2.39.0, using seed 20260529 (`outputs/tables/oct/00_session.txt`).

R packages: brms, cmdstanr, tidyverse, posterior, bayesplot, loo, patchwork, scales, gt, webshot2, here, glue, knitr, PRISMA2020, robvis and readxl. Install CmdStan once with `cmdstanr::install_cmdstan()`. Exporting gt tables to PDF through webshot2 also needs a Chrome or Chromium browser. Rendering the documents requires [Quarto](https://quarto.org), and the PDF output needs a LaTeX installation with xelatex.

MCMC results are reproducible to Monte Carlo error. Exact draws depend on the CmdStan version and platform.

## Notes on the literature search

The update searched PubMed, Embase (Ovid) and CENTRAL from August 2022 to October 2026, plus a hand search. The strategies are in Supplemental Appendix S1 and `CJC/oct/search_strategies_embase_central.md`. The screening decisions for every record are in `CJC/oct/search_update_*.csv`.

The raw Embase and CENTRAL export files are not included, because of database licence terms. `search_merge_dedup.R` documents how they were merged and de-duplicated. `prisma_flow_update.R` builds the flow diagram from the merged log, which is included.

## Statistical workflow

The analysis follows the Bayesian workflow of Gelman et al. (2020, [arXiv:2011.01808](https://arxiv.org/abs/2011.01808)):

1. data and reconciliation;
2. prior specification;
3. prior predictive checks;
4. fake-data simulation and parameter recovery;
5. model fitting;
6. convergence diagnostics;
7. posterior predictive checks;
8. LOO-CV model comparison;
9. posterior summaries;
10. sensitivity analyses and node-splitting;
11. figures;
12. tables and exports.
