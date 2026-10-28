# Title-level screening decisions for unique Embase and CENTRAL records
a <- read.csv("search_update_all_sources_log.csv", stringsAsFactors = FALSE)
todo <- a[a$stage_decision == "TO SCREEN", c("source", "id", "title")]
R <- c(
 A = "Not an ACS-PCI population (neurovascular, stroke, GI/endoscopy, surgery, peripheral or other)",
 B = "Not a randomized comparison (observational, registry, survey, case report or abstract compilation)",
 C = "Review, meta-analysis, guidance or commentary",
 D = "Not a comparison of P2Y12 inhibitors as standard DAPT (de-escalation, guided therapy, monotherapy or other antithrombotic)",
 E = "Pharmacodynamic, mechanistic, angiographic or in vitro end point only",
 F = "Elective PCI or chronic coronary syndrome",
 G = "Erratum or correction",
 H = "Secondary analysis, design paper or earlier report of a trial already included or previously assessed",
 J = "Not a P2Y12 inhibitor comparison (stent, procedural or other intervention)")
emb <- list(
 A = c(149,5,6,7,9,11,14,21,22,24,25,31,34,36,44,45,49,51,53,56,59,60,64,65,66,68,69,71,78,82,83,86,89,91,93,94,98,102,103,107,112,114,115,118,129,132,133,134,143,144,146,151,152,160,171,177,189,196,197,198,201,203),
 B = c(4,10,15,16,17,18,23,26,27,28,29,30,32,33,35,38,41,46,47,48,50,54,57,61,62,63,67,72,73,76,79,80,81,84,88,95,96,97,99,100,101,104,105,109,113,121,123,124,127,128,135,136,138,139,141,147,148,150,154,155,157,159,161,163,164,165,167,170,172,173,174,179,182,183,185,188,192,194,195,199,200,207),
 C = c(2,3,8,19,20,74,75,110,111,120,145,158,175,180,187,190,193,205),
 D = c(1,13,42,85,92,106,125),
 E = c(40,55,58,70,119,131,142,168,169,178,181,186,202),
 F = c(87),
 G = c(52,116))
cen <- list(
 A = c("CN-02269217","CN-01469079","CN-01703403"),
 B = c("CN-02243043","CN-01088389","CN-01068615","CN-01719653"),
 D = c("CN-02233635","CN-02499776","CN-00905986","CN-02117978","CN-01798616","CN-01810103","CN-01036985"),
 J = c("CN-01420382","CN-01339323","CN-01857288","CN-02036024","CN-01975641","CN-02579392"),
 E = c("CN-03028116"),
 F = c("CN-01101370","CN-00842378"),
 H = c("CN-02101150","CN-00803002","CN-00682219","CN-00733635","CN-02523480"))
ids <- c(unlist(lapply(names(emb), function(k) setNames(rep(k, length(emb[[k]])), paste0("EMB", emb[[k]])))),
         unlist(lapply(names(cen), function(k) setNames(rep(k, length(cen[[k]])), cen[[k]]))))
stopifnot(!any(duplicated(names(ids))))
miss <- setdiff(todo$id, names(ids)); extra <- setdiff(names(ids), todo$id)
cat("Unassigned:", miss, "\nAssigned but not TO SCREEN (duplicates):", extra, "\n")
dec <- data.frame(id = todo$id, source = todo$source, title = todo$title,
                  stage_decision = "SCREEN-EXCLUDED", reason = unname(R[ids[todo$id]]))
write.csv(dec, "search_update_embase_central_decisions.csv", row.names = FALSE)
print(table(dec$source, substr(dec$reason, 1, 40)))
