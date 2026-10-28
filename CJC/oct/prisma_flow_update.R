# =============================================================================
# PRISMA 2020 flow diagram — updated Bayesian NMA (search 2022-09-01 to 2026-10-20)
# Requires: install.packages("PRISMA2020")
# Run from the CJC/oct folder:  source("prisma_flow_update.R")
# Counts come from search_update_all_sources_log.csv (PubMed + Embase + CENTRAL,
# built by search_merge_dedup.R) or, if absent, the PubMed-only screening log
# =============================================================================
library(PRISMA2020)

# ---- locate CJC/oct (works from CJC/oct or from the dapt_ma project root) ------
.old_wd <- getwd()
for (d in c(".", "CJC/oct", "../CJC/oct", "../../CJC/oct")) {
  if (file.exists(file.path(d, "search_update_screening_log.csv"))) { setwd(d); break }
}
if (!file.exists("search_update_screening_log.csv"))
  stop("Cannot find CJC/oct; setwd() to that folder first. Current directory: ", getwd())
message("Working in: ", getwd())

# Start from the package's own template so all required rows are present
csv  <- read.csv(system.file("extdata", "PRISMA.csv", package = "PRISMA2020"),
                 stringsAsFactors = FALSE)
set_n <- function(df, field, value) { df$n[which(df$data == field)] <- value; df }
set_box <- function(df, field, text) { df$boxtext[which(df$data == field)] <- text; df }

# ---- counts ------------------------------------------------------------------
# If search_merge_dedup.R has been run (PubMed + Embase + CENTRAL), counts are
# computed from search_update_all_sources_log.csv; otherwise the PubMed-only
# counts below are used.
ft_cat <- function(r) {
  r <- sub("^Full text: ", "", r)
  ifelse(grepl("pharmacodynamic|mechanistic", r, TRUE), "Pharmacodynamic or mechanistic end point only",
  ifelse(grepl("follow-up", r, TRUE), "Follow-up <6 months",
  ifelse(grepl("MACE not reported", r, TRUE), "MACE not reported", sub(" \\(.*$", "", r))))
}
# Screening (title/abstract) exclusions grouped into major categories
screen_cat <- function(r) {
  r <- tolower(r)
  ifelse(grepl("secondary|post hoc|design paper|extended follow-up|substudy", r), "Secondary report of an existing trial",
  ifelse(grepl("not an acs-pci|not coronary|post-cabg|elective|chronic coronary|stable patients", r), "Population not ACS with PCI",
  ifelse(grepl("pharmacodynamic|mechanistic|in vitro|angiographic", r), "Pharmacodynamic or mechanistic only",
  ifelse(grepl("not a randomized", r), "Not a randomized comparison",
  ifelse(grepl("review|meta-analysis|guidance", r), "Review or meta-analysis",
  ifelse(grepl("not a comparison|not a p2y12|de-escalation|monotherapy|guided|stent randomization|other drug", r),
         "Not a P2Y12 comparison within DAPT",
         "Other"))))))
}
sx <- NULL
if (file.exists("search_update_all_sources_log.csv")) {
  lg <- read.csv("search_update_all_sources_log.csv", stringsAsFactors = FALSE)
  if (any(lg$stage_decision == "TO SCREEN"))
    stop("search_update_all_sources_log.csv still has records marked TO SCREEN")
  src <- table(factor(lg$source, c("PubMed", "Embase", "CENTRAL")))
  uq  <- lg[lg$stage_decision != "DUPLICATE", ]
  ft  <- uq[uq$stage_decision == "FT-EXCLUDED", ]
  fr  <- sort(table(ft_cat(ft$reason)), decreasing = TRUE)
  sx  <- table(screen_cat(uq$reason[uq$stage_decision == "SCREEN-EXCLUDED"]))
  sx  <- c(sort(sx[names(sx) != "Other"], decreasing = TRUE),
           sx[names(sx) == "Other"])
  print(sx)
  db  <- list(
    database_results          = sum(src),
    database_specific_results = paste0(c("PubMed/MEDLINE", "Embase", "CENTRAL"), ", ", src, collapse = "; "),
    duplicates                = sum(lg$stage_decision == "DUPLICATE"),
    records_screened          = nrow(uq),
    records_excluded          = sum(uq$stage_decision == "SCREEN-EXCLUDED"),
    dbr_sought_reports        = sum(uq$stage_decision != "SCREEN-EXCLUDED"),
    dbr_assessed              = sum(uq$stage_decision != "SCREEN-EXCLUDED"),
    dbr_excluded              = paste0(names(fr), ", ", fr, collapse = "; "),
    db_reports_included       = sum(grepl("^INCLUDED", uq$stage_decision)))
} else {
  db <- list(
    database_results = 71, database_specific_results = "PubMed/MEDLINE, 71",
    duplicates = 0, records_screened = 71, records_excluded = 59,
    dbr_sought_reports = 12, dbr_assessed = 12,
    dbr_excluded = "Pharmacodynamic or mechanistic end point only, 7; Follow-up <6 months, 2; MACE not reported, 1",
    db_reports_included = 2)
}
hand_reports <- 1            # SWITCH-SWEDEHEART (hand search)
new_studies  <- 3            # TUXEDO-2, TC4, SWITCH-SWEDEHEART (update if a new trial is found)

counts <- list(
  # Previous version of the review (Kutcher 2023 BNMA)
  previous_studies   = 17,
  previous_reports   = 17,     # [AUTHOR TO CONFIRM number of reports]
  # Identification: databases
  database_results   = db$database_results,
  database_specific_results = db$database_specific_results,
  register_results   = 0,
  duplicates         = db$duplicates,
  excluded_automatic = 0,
  excluded_other     = 0,
  # Screening
  records_screened   = db$records_screened,
  records_excluded   = db$records_excluded,
  dbr_sought_reports = db$dbr_sought_reports,
  dbr_notretrieved_reports = 0,
  dbr_assessed       = db$dbr_assessed,
  dbr_excluded       = db$dbr_excluded,
  # Identification: other methods (hand search)
  website_results      = 0,
  organisation_results = 0,
  citations_results    = hand_reports,   # relabelled below as "Hand searching"
  other_sought_reports = hand_reports,
  other_notretrieved_reports = 0,
  other_assessed       = hand_reports,
  other_excluded       = 0,
  # Included
  new_studies   = new_studies,
  new_reports   = db$db_reports_included + hand_reports,
  total_studies = 17 + new_studies,
  total_reports = 17 + db$db_reports_included + hand_reports  # [AUTHOR TO CONFIRM previous_reports]
)
for (f in names(counts)) csv <- set_n(csv, f, counts[[f]])
csv <- set_box(csv, "citations_results", "Hand searching")

prisma_data <- PRISMA_data(csv)
# Add the screening-exclusion categories to the "Records excluded" box.
# PRISMA2020 builds that label as: <text> "\n(n = " <number> ")", so the
# breakdown is appended to the number and the package supplies the final ")".
if (!is.null(sx))
  prisma_data$records_excluded <- paste0(
    prisma_data$records_excluded, "):\n",
    paste0(names(sx), " (n = ", sx, collapse = ")\n"))

plot <- PRISMA_flowdiagram(
  prisma_data,
  interactive = FALSE,
  previous    = TRUE,    # show "previous studies" column (updated review)
  other       = TRUE,    # show "other methods" (hand search) column
  detail_databases = TRUE, # list records per database
  fontsize    = 10,
  font        = "Helvetica"
)

dir.create("figures", showWarnings = FALSE)
PRISMA_save(plot, filename = "figures/PRISMA_flow_update.png", filetype = "PNG", overwrite = TRUE)
PRISMA_save(plot, filename = "figures/PRISMA_flow_update.pdf", filetype = "PDF", overwrite = TRUE)
PRISMA_save(plot, filename = "figures/PRISMA_flow_update.svg", filetype = "SVG", overwrite = TRUE)
setwd(.old_wd)   # restore the original working directory
plot
