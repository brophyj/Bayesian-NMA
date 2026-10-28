# =============================================================================
# Merge and de-duplicate the updated searches (PubMed + Embase + CENTRAL)
# Run from CJC/oct:  source("search_merge_dedup.R")
# Inputs : search_update_screening_log.csv  (71 PubMed records + hand search, screened)
#          embase.csv                        (Ovid Embase export, CSV; 1st row = copyright; embase.xls also accepted)
#          cochrane.csv                      (Cochrane Library CENTRAL export, CSV)
#          search_update_embase_central_decisions.csv
#                                            (screening decisions for unique non-PubMed records)
# Outputs: search_update_all_sources_log.csv (every record: source, duplicate status, decision)
#          search_update_new_to_screen.csv   (unique non-PubMed records still without a decision)
# Base R only (readxl needed only if reading embase.xls)
# =============================================================================


# ---- locate CJC/oct (works from CJC/oct or from the dapt_ma project root) ------
.old_wd <- getwd()
for (d in c(".", "CJC/oct", "../CJC/oct", "../../CJC/oct")) {
  if (file.exists(file.path(d, "search_update_screening_log.csv"))) { setwd(d); break }
}
if (!file.exists("search_update_screening_log.csv"))
  stop("Cannot find CJC/oct; setwd() to that folder first. Current directory: ", getwd())
message("Working in: ", getwd())

norm_doi   <- function(d) { d <- tolower(trimws(d)); d <- sub("^https?://(dx\\.)?doi\\.org/", "", d); sub("[.]$", "", d) }
norm_title <- function(t) {
  t <- tolower(iconv(t, to = "ASCII//TRANSLIT", sub = ""))
  t <- sub("\\[[a-z]+\\]\\s*$", "", t)                    # trailing [Chinese], [Russian] ...
  t <- gsub("[^a-z0-9]", "", t); substr(t, 1, 70)
}
blank <- function(x) is.na(x) | trimws(x) == ""

# ---- PubMed (already screened) -----------------------------------------------
pm <- read.csv("search_update_screening_log.csv", stringsAsFactors = FALSE)
pm_db <- pm[pm$no != "hand", ]
cit <- sub("^[0-9]+: ", "", pm_db$citation)
# title = text after the author list ("... ; <group>. Title. Journal ...")
pm_title <- vapply(cit, function(s) {
  s2 <- sub("^.*?[A-Z]{1,3}( Jr)?(; [^.]*)?\\. ", "", s)  # drop authors (initials end with ". ")
  sub("\\. [A-Z][^.]*\\. (19|20)[0-9]{2}.*$", "", s2)      # drop journal/date onward
}, "", USE.NAMES = FALSE)
pubmed <- data.frame(
  source = "PubMed", id = as.character(pm_db$pmid), title = pm_title, journal = NA, year = NA,
  type = NA, doi = vapply(cit, function(s) {
    m <- regmatches(s, regexpr("doi: 10\\.[^ ]+", s))
    if (length(m)) sub("[.]$", "", sub("^doi: ", "", m)) else NA_character_ }, "", USE.NAMES = FALSE),
  pmid = as.character(pm_db$pmid), registry = FALSE,
  stage_decision = pm_db$stage_decision, reason = pm_db$reason, stringsAsFactors = FALSE)

# ---- Embase (Ovid export: CSV preferred, Excel as fallback; 1st row = copyright) --
e <- if (file.exists("embase.csv")) {
  read.csv("embase.csv", skip = 1, stringsAsFactors = FALSE, check.names = FALSE,
           encoding = "UTF-8")
} else {
  if (!requireNamespace("readxl", quietly = TRUE)) stop("embase.csv not found; install readxl to read embase.xls")
  as.data.frame(readxl::read_excel("embase.xls", skip = 1))
}
e <- e[!blank(e$TI), ]
embase <- data.frame(
  source = "Embase", id = paste0("EMB", e$ORN), title = e$TI, journal = e$SO,
  year = regmatches(e$SO, regexpr("(19|20)[0-9]{2}", e$SO)), type = e$PT,
  doi = NA_character_, pmid = NA_character_, registry = FALSE,
  stage_decision = NA_character_, reason = NA_character_, stringsAsFactors = FALSE)

# ---- CENTRAL (Cochrane Library CSV export) -------------------------------------
cc <- read.csv("cochrane.csv", stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8")
central <- data.frame(
  source = "CENTRAL", id = cc[["CENTRAL ID"]], title = cc$Title, journal = cc$Source,
  year = as.character(cc$Year), type = cc[["Publication Type"]],
  doi = ifelse(blank(cc$DOI), NA, cc$DOI),
  pmid = ifelse(blank(cc[["PubMed ID"]]), NA, as.character(cc[["PubMed ID"]])),
  registry = !blank(cc[["CTgov ID"]]) | !blank(cc[["ICTRP ID"]]) |
             grepl("clinicaltrials\\.gov|trialsearch\\.who|ICTRP|CTRI|ChiCTR|ISRCTN|DRKS|JPRN|IRCT",
                   paste(cc$Source, cc$URL), ignore.case = TRUE),
  stage_decision = NA_character_, reason = NA_character_, stringsAsFactors = FALSE)

all <- rbind(pubmed, embase, central)

# ---- de-duplicate (PubMed first, then Embase, then CENTRAL) --------------------
all$k_doi   <- ifelse(blank(all$doi),  NA, norm_doi(all$doi))
all$k_pmid  <- ifelse(blank(all$pmid), NA, all$pmid)
all$k_title <- ifelse(blank(all$title) | nchar(norm_title(all$title)) < 25, NA, norm_title(all$title))
all$dup_of  <- NA_character_
for (i in seq_len(nrow(all))[-1]) {
  prev <- seq_len(i - 1)
  prev <- prev[is.na(all$dup_of[prev])]                    # match against retained records only
  hit <- prev[(!is.na(all$k_doi[i])   & all$k_doi[prev]   %in% all$k_doi[i])   |
              (!is.na(all$k_pmid[i])  & all$k_pmid[prev]  %in% all$k_pmid[i])  |
              (!is.na(all$k_title[i]) & all$k_title[prev] %in% all$k_title[i])]
  if (length(hit) && all$source[i] != "PubMed") {
    all$dup_of[i] <- paste(all$source[hit[1]], all$id[hit[1]])
    all$stage_decision[i] <- "DUPLICATE"
    all$reason[i] <- paste("Duplicate of", all$dup_of[i])
  }
}
# second pass: title (first 45 normalized characters) found inside a full PubMed citation
pm_cit_norm <- gsub("[^a-z0-9]", "", tolower(iconv(cit, to = "ASCII//TRANSLIT", sub = "")))
for (i in which(all$source != "PubMed" & is.na(all$dup_of) & !is.na(all$k_title))) {
  h <- which(grepl(substr(all$k_title[i], 1, 45), pm_cit_norm, fixed = TRUE))
  if (length(h)) {
    all$dup_of[i] <- paste("PubMed", pubmed$id[h[1]])
    all$stage_decision[i] <- "DUPLICATE"
    all$reason[i] <- paste("Duplicate of", all$dup_of[i])
  }
}

# ---- screening decisions for unique Embase/CENTRAL records ----------------------
if (file.exists("search_update_embase_central_decisions.csv")) {
  dec <- read.csv("search_update_embase_central_decisions.csv", stringsAsFactors = FALSE)
  j <- match(all$id, dec$id)
  use <- !is.na(j) & is.na(all$stage_decision)
  all$stage_decision[use] <- dec$stage_decision[j[use]]
  all$reason[use]         <- dec$reason[j[use]]
}
all$stage_decision[is.na(all$stage_decision)] <- "TO SCREEN"

# ---- outputs -----------------------------------------------------------------
keep <- c("source", "id", "title", "journal", "year", "type", "doi", "pmid", "registry",
          "dup_of", "stage_decision", "reason")
write.csv(all[, keep], "search_update_all_sources_log.csv", row.names = FALSE)
write.csv(all[all$stage_decision == "TO SCREEN", keep], "search_update_new_to_screen.csv", row.names = FALSE)

tab <- table(factor(all$source, c("PubMed", "Embase", "CENTRAL")))
cat("\nRecords identified:\n"); print(tab)
cat("Duplicates removed:", sum(all$stage_decision == "DUPLICATE"), "\n")
cat("Unique records screened:", sum(all$stage_decision != "DUPLICATE"), "\n")
print(table(all$source, all$stage_decision))
cat("Still to screen:", sum(all$stage_decision == "TO SCREEN"), "\n")

setwd(.old_wd)   # restore the original working directory
