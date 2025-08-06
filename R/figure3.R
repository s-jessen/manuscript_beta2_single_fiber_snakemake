# ---- Setup ----
source("renv/activate.R")
source(snakemake@input[["functions"]])
source(snakemake@input[["settings"]])

# Load libraries
library(tidyverse)
library(ggvenn)
library(ggrepel)
library(patchwork)

# Load data --------------------------------------------------------------
df_long <- readRDS(snakemake@input[["df_long"]])
df_long_l2fc <- readRDS(snakemake@input[["df_long_l2fc"]])

# GSEA's
gsea_i_and_ii_interaction <- readRDS(snakemake@input[["gsea_i_and_ii_interaction"]])
gsea_i_interaction <- readRDS(snakemake@input[["gsea_i_interaction"]])
gsea_ii_interaction <- readRDS(snakemake@input[["gsea_ii_interaction"]])


# Load results  ----------------------------------------------------------
# Interaction between intervention groups
results_i_and_ii_interaction <- readRDS(snakemake@input[["results_i_and_ii_interaction"]])
results_i_interaction <- readRDS(snakemake@input[["results_i_interaction"]])
results_ii_interaction <- readRDS(snakemake@input[["results_ii_interaction"]])


# Volcanoes --------------------------------------------------------------
cat("\033[1;33mCreating volcano plots (Fig. 3)\033[0m\n")
# Independent, intervention x time interaction
labels_i_and_ii <- list(
  "MUSTN1",
  "VIM",
  "THY1",
  "MYBPH",
  "UFL1",
  "RPL35A",
  "MYLK2",
  "API5"
)

volcano_i_and_ii_interaction <- volcano_plot(results_i_and_ii_interaction, labels_i_and_ii) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "effect size (l2fc RES vs. l2fc B2A)",
    title = "Two-way interaction effect (fiber type independent)\n(intervention × time)"
  )

# Type I, intervention x time interaction
labels_i <- list(
  "MYH8",
  "ATP2B1",
  "COQ7",
  "NDUFA4",
  "CREG1",
  "TUBB4A",
  "ACTN1",
  "ACTN4",
  "SH3BGRL3",
  "DSP"
)

volcano_i_interaction <- volcano_plot(results_i_interaction, labels_i) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "effect size (l2fc RES vs. l2fc B2A)",
    title = "Two-way interaction effect in type I fibers\n(intervention × time)"
  )

# Type II, intervention x time interaction
labels_ii <- list(
  "NUDT5",
  "API5",
  "PRRC1",
  "MYLK2",
  "MYH6",
  "MYH4",
  "COL1A2",
  "ACTC1",
  "MYL2",
  "MYBPH",
  "MUSTN1",
  "MYH3",
  "TNNC1",
  "TNNI1",
  "CSRP3",
  "VIM",
  "RPS5",
  "CRYAB"
)

volcano_ii_interaction <- volcano_plot(results_ii_interaction, labels_ii) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "effect size (l2fc RES vs. l2fc B2A)",
    title = "Two-way interaction effect in type II fibers\n(intervention × time)"
  )

# Collect volcanoes into single layout
plot_all_volcanoes <- (volcano_i_and_ii_interaction | volcano_i_interaction | volcano_ii_interaction)

ggsave(snakemake@output[["fig_volcanoes"]], width = 180, height = 75, units = "mm", plot = plot_all_volcanoes)


# GSEA's -----------------------------------------------------------------
cat("\033[1;33mCreating GSEA plots (Fig. 3)\033[0m\n")

# Select terms of interest, and terms to highlight for each comparison
terms_i_and_ii_interaction <- c(
  "extracellular matrix",
  "actin cytoskeleton organization",
  "cell-cell adhesion",
  "cytoskeleton organization",
  "cellular developmental process",
  "intermediate filament",
  "ADP metabolic process",
  "pyruvate metabolic process"
)

highlight_i_and_ii_interaction <- c(
  "extracellular matrix",
  "cytoskeleton organization"
)

terms_i_interaction <- c(
  "cell adhesion",
  "extracellular matrix",
  "ribosomal small subunit biogenesis",
  "cytoskeleton organization",
  "cellular developmental process",
  "intermediate filament",
  "cellular respiration",
  "mitochondrial envelope"
)

highlight_i_interaction <- c(
  "mitochondrial envelope",
  "ribosomal small subunit biogenesis",
  "extracellular matrix"
)

terms_ii_interaction <- c(
  "muscle contraction",
  "muscle system process",
  "muscle structure development",
  "mitochondrial matrix",
  "protein refolding",
  "actin filament-based process",
  "cytosolic ribosome",
  "rRNA metabolic process",
  "extracellular matrix disassembly"
)

highlight_ii_interaction <- c(
  "muscle structure development",
  "mitochondrial matrix",
  "muscle contraction",
  "cytosolic ribosome"
)


# Collect to iterate over
gseas_list <- list(
  "i_and_ii" = gsea_i_and_ii_interaction,
  "i" = gsea_i_interaction,
  "ii" = gsea_ii_interaction
)

terms_list <- list(
  "i_and_ii" = terms_i_and_ii_interaction,
  "i" = terms_i_interaction,
  "ii" = terms_ii_interaction
)

highlights_list <- list(
  "i_and_ii" = highlight_i_and_ii_interaction,
  "i" = highlight_i_interaction,
  "ii" = highlight_ii_interaction
)

# list to store plots
gsea_plots <- list()

for (comparison in names(gseas_list)) {
  df <- gseas_list[[comparison]]
  terms <- terms_list[[comparison]]
  highlight <- highlights_list[[comparison]]

  df_plot <- df %>%
    as.data.frame() %>%
    dplyr::mutate(color = ifelse(Description %in% highlight, "highlight", "no_color")) %>%
    dplyr::filter(Description %in% terms) %>%
    ggplot(aes(x = -log10(pvalue), y = reorder(Description, -pvalue), fill = color)) +
    geom_col() +
    theme(
      legend.position = "none"
    ) +
    scale_fill_manual(values = c("highlight" = "#a51e22", "no_color" = "#d3d3d3")) +
    labs(
      y = NULL,
      x = "-log10(p)",
      title = "GSEA"
    )

  gsea_plots[[comparison]] <- df_plot
}

# Combine plots
plot_all_gseas <- (gsea_plots[["i_and_ii"]] | gsea_plots[["i"]] | gsea_plots[["ii"]])

# Save plot
ggplot2::ggsave(snakemake@output[["fig_gsea"]], height = 50, width = 200, units = "mm", plot = plot_all_gseas)


# Cytosolic ribosomal proteins -------------------------------------------
cat("\033[1;33mCreating Ribosomal cytosolic proteins plots (Fig. 3)\033[0m\n")

# Terms of interest
terms <- c("GO:0006413", "GO:0006414", "GO:0022625", "GO:0022627")
term_labels <- c(
  "GO:0006413" = "Init.\nfactors",
  "GO:0006414" = "Elong.\nfactors",
  "GO:0022627" = "Small ribo.\nsubunits",
  "GO:0022625" = "Large ribo.\nsubunits"
)

# Set order
term_order <- c(
  "GO:0006413",
  "GO:0006414",
  "GO:0022627",
  "GO:0022625"
)

# Retrieve plot and results
ribosomal <- plot_terms(terms, term_labels, term_order)

plot_ribo <- ribosomal[["plot"]] +
  coord_cartesian(ylim = c(-0.15, 0.6))

# Save plot
ggplot2::ggsave(snakemake@output[["fig_ribosomal_proteins"]], width = 200, height = 50, units = "mm", plot = plot_ribo)


# Mitochondrial proteins -------------------------------------------------
cat("\033[1;33mCreating Mitochondrial complex proteins plots (Fig. 3)\033[0m\n")

# Complexes to iterate over
complexes <- c("CI", "CII", "CIII", "CIV", "CV")

# Retrieve plot and results
mito <- plot_complexes(complexes)

plot_mito <- mito[["plot"]] +
  coord_cartesian(ylim = c(-0.6, 0.3))

# Save plot
ggplot2::ggsave(snakemake@output[["fig_mitochondrial_proteins"]], width = 200, height = 50, units = "mm", plot = plot_mito)


# Structural GO-terms ----------------------------------------------------
cat("\033[1;33mCreating structural GO-term plots (Fig. 3)\033[0m\n")

# Choose GO:terms to analyze
terms <- c("GO:0031012", "GO:0006936", "GO:0030239", "GO:0030036")
term_labels <- c(
  "GO:0031012" = "ECM",
  "GO:0006936" = "Muscle\ncontraction",
  "GO:0030239" = "Myofibril\nassembly",
  "GO:0030036" = "Actin\ncytosk.org."
)

# Set order
term_order <- c(
  "GO:0030036",
  "GO:0006936",
  "GO:0031012",
  "GO:0030239"
)

# Retrieve plot and results
structural <- plot_terms(terms, term_labels, term_order)

plot_structural <- structural[["plot"]] +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  ))

# Save plot
ggplot2::ggsave(snakemake@output[["fig_structural_go_terms"]], width = 100, height = 50, units = "mm", plot = plot_structural)


# Metabolic GO-terms -----------------------------------------------------
cat("\033[1;33mCreating metabolic GO-term plots (Fig. 3)\033[0m\n")

# Choose GO:terms to analyze
terms <- c("GO:0006099", "GO:0006635", "GO:0061621")
term_labels <- c(
  "GO:0006099" = "TCA",
  "GO:0006635" = "FA ox.",
  "GO:0061621" = "Gly-\ncolysis"
)

# Set order
term_order <- c(
  "GO:0006635",
  "GO:0061621",
  "GO:0006099"
)

# Retrieve plot and results
metabolic <- plot_terms(terms, term_labels, term_order)

plot_metabolic <- metabolic[["plot"]] +
  coord_cartesian(ylim = c(-0.35, 0.3))

# Save plot
ggsave(snakemake@output[["fig_metabolic_go_terms"]], width = 100, height = 50, units = "mm", plot = plot_metabolic)
