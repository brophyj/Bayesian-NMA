# =============================================================================
# Central Illustration (draft) — Bayesian NMA of P2Y12 inhibitors after ACS
# Run from CJC/oct or the dapt_ma project root:  source("CJC/oct/central_illustration.R")
# Reads outputs/tables/oct/02_main_results_raw.csv and 00_input_mace.csv
# Writes CJC/oct/figures/Central_Illustration.{png,tiff,pdf}
# Requires: ggplot2, dplyr, patchwork
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(patchwork) })

.old_wd <- getwd()
for (d in c(".", "CJC/oct", "../CJC/oct", "../../CJC/oct")) {
  if (file.exists(file.path(d, "manuscript_oct2026.qmd"))) { setwd(d); break }
}
tdir <- file.path("..", "..", "outputs", "tables", "oct")
res  <- read.csv(file.path(tdir, "02_main_results_raw.csv"), stringsAsFactors = FALSE)
inp  <- read.csv(file.path(tdir, "00_input_mace.csv"), stringsAsFactors = FALSE)

# ---- colours: diverging (benefit = cool, equivalence = neutral, harm = warm) ----
col <- c(benefit = "#2B6CB0", rope = "#C9C9C9", harm = "#C2502A")
ink <- "#1F2328"; ink2 <- "#57606A"
node_col <- c(C = "#6E7781", T = "#8250DF", P = "#0F7B6C")

# ---- network numbers ---------------------------------------------------------
K     <- n_distinct(inp$study)
N_pat <- sum(inp$n)
pairs <- inp |> group_by(study) |> summarise(pair = paste(sort(unique(treat)), collapse = "-"), .groups = "drop") |>
  count(pair)
n_by  <- inp |> group_by(treat) |> summarise(n = sum(n))
k_of  <- function(p) pairs$n[pairs$pair == p]

# ---- Panel A: the evidence (network) -------------------------------------------
nodes <- data.frame(treat = c("C", "T", "P"), x = c(0, -1, 1), y = c(0, 1.55, 1.55),
                    lab = c("Clopidogrel", "Ticagrelor", "Prasugrel")) |>
  left_join(n_by, by = "treat")
edges <- data.frame(from = c("C", "C", "P"), to = c("T", "P", "T"),
                    k = c(k_of("C-T"), k_of("C-P"), k_of("P-T"))) |>
  left_join(nodes[, c("treat", "x", "y")], by = c(from = "treat")) |>
  left_join(nodes[, c("treat", "x", "y")], by = c(to = "treat"), suffix = c("", "end"))
pA <- ggplot() +
  geom_segment(data = edges, aes(x, y, xend = xend, yend = yend, linewidth = k), colour = "#D0D7DE",
               lineend = "round") +
  geom_label(data = edges, aes((x + xend) / 2, (y + yend) / 2, label = paste(k, "RCTs")),
             size = 3.1, colour = ink2, label.size = 0, fill = "white") +
  geom_point(data = nodes, aes(x, y, size = n, colour = treat)) +
  geom_text(data = nodes, aes(x, y - ifelse(treat == "C", 0.42, -0.40), label = lab),
            size = 3.6, fontface = "bold", colour = ink) +
  scale_linewidth(range = c(2, 6), guide = "none") + scale_size(range = c(9, 17), guide = "none") +
  scale_colour_manual(values = node_col, guide = "none") +
  coord_equal(xlim = c(-1.6, 1.6), ylim = c(-0.6, 2.1), clip = "off") +
  labs(title = "The evidence",
       subtitle = paste0(K, " randomized trials, ", format(N_pat, big.mark = ","), " patients\nACS, all regimens with aspirin")) +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", colour = ink, size = 13),
        plot.subtitle = element_text(colour = ink2, size = 9.5, margin = margin(b = 6)))

# ---- Panel B: posterior probabilities -----------------------------------------
lab <- c("T vs C" = "Ticagrelor vs clopidogrel", "P vs C" = "Prasugrel vs clopidogrel",
         "P vs T" = "Prasugrel vs ticagrelor")
pb <- res |>
  mutate(Outcome = ifelse(grepl("MACE", outcome), "MACE", "Major bleeding"),
         Comparison = factor(lab[comparison], levels = rev(lab))) |>
  select(Outcome, Comparison, benefit = pr_clin_benefit, rope = pr_rope, harm = pr_clin_harm) |>
  tidyr::pivot_longer(c(benefit, rope, harm), names_to = "cat", values_to = "p") |>
  mutate(cat = factor(cat, levels = c("harm", "rope", "benefit"))) |>
  group_by(Outcome, Comparison) |> arrange(desc(cat), .by_group = TRUE) |>
  mutate(xmax = cumsum(p), xmin = xmax - p, mid = (xmin + xmax) / 2, yn = as.numeric(Comparison)) |> ungroup()
pB <- ggplot(pb) +
  geom_rect(aes(xmin = xmin + 0.004, xmax = xmax - 0.004, ymin = yn - 0.32, ymax = yn + 0.32, fill = cat)) +
  geom_text(data = filter(pb, p >= 0.07), aes(x = mid, y = yn, label = paste0(round(100 * p), "%"),
            colour = cat == "rope"), size = 3.3, fontface = "bold") +
  scale_fill_manual(values = col, breaks = c("benefit", "rope", "harm"),
                    labels = c("Meaningful benefit (HR < 0.90)", "Practical equivalence (HR 0.90 to 1.11)",
                               "Meaningful harm (HR > 1.11)"), name = NULL) +
  scale_colour_manual(values = c(`TRUE` = ink, `FALSE` = "white"), guide = "none") +
  scale_x_continuous(expand = c(0, 0), labels = NULL, breaks = NULL) +
  scale_y_continuous(breaks = seq_along(lab), labels = rev(unname(lab)), expand = c(0.02, 0.02)) +
  facet_wrap(~ Outcome, ncol = 1) +
  labs(title = "What the Bayesian analysis shows",
       subtitle = "Posterior probability that the first-named drug is meaningfully better,\nequivalent or meaningfully worse", x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", colour = ink, size = 13),
        plot.subtitle = element_text(colour = ink2, size = 9.5),
        strip.text = element_text(face = "bold", hjust = 0, colour = ink, size = 11),
        axis.text.y = element_text(colour = ink, size = 10), panel.grid = element_blank(),
        legend.position = "bottom", legend.text = element_text(size = 8.5, colour = ink2),
        legend.key.size = unit(0.4, "cm"))

# ---- Panel C: clinical message -------------------------------------------------
P <- function(cmp, o) res[grepl(o, res$outcome) & res$comparison == cmp, ]
msg <- c(
  sprintf("No definitive hierarchy: the probability\nthat prasugrel is meaningfully better\nthan ticagrelor for MACE is %d%%,\nroughly a coin toss", round(100 * P("P vs T", "MACE")$pr_clin_benefit)),
  sprintf("Ticagrelor and clopidogrel are probably\nequivalent for MACE (%d%%), but ticagrelor\nprobably causes more bleeding (%d%%)",
          round(100 * P("T vs C", "MACE")$pr_rope), round(100 * P("T vs C", "Bleeding")$pr_clin_harm)),
  "Clopidogrel remains a reasonable\nchoice, particularly for older and\nhigh-bleeding-risk patients")
msgdf <- data.frame(y = c(3, 2, 1), txt = msg)
pC <- ggplot(msgdf) +
  geom_point(aes(0, y), shape = 21, size = 7, fill = "#F6F8FA", colour = "#D0D7DE", stroke = 0.6) +
  geom_text(aes(0, y, label = 4 - y), fontface = "bold", colour = ink, size = 3.6) +
  geom_text(aes(0.16, y, label = txt), hjust = 0, vjust = 0.5, lineheight = 0.95, size = 3.3, colour = ink) +
  coord_cartesian(xlim = c(-0.1, 1.25), ylim = c(0.45, 3.55), clip = "off") +
  labs(title = "Clinical message") + theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", colour = ink, size = 13))

p <- (pA | pB | pC) + plot_layout(widths = c(0.9, 1.45, 1.05)) +
  plot_annotation(
    title = "Revisiting P2Y12 inhibitors after acute coronary syndromes: a Bayesian network meta-analysis",
    caption = "MACE, death, myocardial infarction or stroke. HR, hazard ratio. Random-effects network meta-analysis with a complementary log-log link and follow-up offset.",
    theme = theme(plot.title = element_text(face = "bold", size = 14.5, colour = ink),
                  plot.caption = element_text(hjust = 0, size = 8.5, colour = ink2),
                  plot.background = element_rect(fill = "white", colour = NA)))

dir.create("figures", showWarnings = FALSE)
ggsave("figures/Central_Illustration.png", p, width = 13, height = 6, dpi = 300, bg = "white")
ggsave("figures/Central_Illustration.tiff", p, width = 13, height = 6, dpi = 300, bg = "white", compression = "lzw")
ggsave("figures/Central_Illustration.pdf", p, width = 13, height = 6, bg = "white", device = cairo_pdf)
setwd(.old_wd)
p
