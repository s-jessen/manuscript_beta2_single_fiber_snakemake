# ---- Setup ----
source("renv/activate.R")
source(snakemake@input[["functions"]])
source(snakemake@input[["settings"]])

# Load libraries
library(tidyverse)
library(ggrepel)
library(patchwork)

# Load data --------------------------------------------------------------
df_long <- readRDS(snakemake@input[["df_long"]])
df_long_l2fc <- readRDS(snakemake@input[["df_long_l2fc"]])
df_long_l2fc_mean <- readRDS(snakemake@input[["df_long_l2fc_mean"]])
gsea_all <- readRDS(snakemake@input[["gsea_all"]])


# Load results  ----------------------------------------------------------
# Intervention group independent
results_main <- readRDS(snakemake@input[["results_main"]])
results_i <- readRDS(snakemake@input[["results_i"]])
results_ii <- readRDS(snakemake@input[["results_ii"]])
results_interaction <- readRDS(snakemake@input[["results_interaction"]])

# Terbutaline
results_ter_main <- readRDS(snakemake@input[["results_ter_main"]])
results_ter_i <- readRDS(snakemake@input[["results_ter_i"]])
results_ter_ii <- readRDS(snakemake@input[["results_ter_ii"]])
results_ter_interaction <- readRDS(snakemake@input[["results_ter_interaction"]])

# Resistance training
results_res_main <- readRDS(snakemake@input[["results_res_main"]])
results_res_i <- readRDS(snakemake@input[["results_res_i"]])
results_res_ii <- readRDS(snakemake@input[["results_res_ii"]])
results_res_interaction <- readRDS(snakemake@input[["results_res_interaction"]])


# Volcano plots ----------------------------------------------------------
cat("\033[1;33mCreating volcano plots (Fig. 2)\033[0m\n")

# Independent, main
labels_main <- c(
  "S100A13",
  "MYBPH",
  "THY1",
  "RRAD",
  "AKR1C3",
  "ANKRD2",
  "VIM",
  "YBX1",
  "SERPINH1",
  "ART3",
  "FAM185A",
  "DHRSC7"
)

volcano_main <- volcano_plot(results_main, labels_main) +
  labs(title = "Main effect")

# Independent, type 1
labels_i <- c(
  "S100A13",
  "THY1",
  "TUBB2B",
  "VIM",
  "SULT1A3",
  "CYC1",
  "DNAJC11"
)

volcano_i <- volcano_plot(results_i, labels_i) +
  labs(title = "Type I")

# Independent, type 2
labels_ii <- c(
  "MYBPH",
  "S100A13",
  "ACTC1",
  "MYL2",
  "THY1",
  "AKR1C3",
  "AKR1C2",
  "SERPINH1",
  "YBX1",
  "SERPINB1",
  "ALDH1A2",
  "ART3"
)

volcano_ii <- volcano_plot(results_ii, labels_ii) +
  labs(title = "Type II")

# Independent, interaction
labels_interaction <- c(
  "MYH6",
  "MYL2",
  "ATP2A2",
  "PPP1R3D",
  "CREG1",
  "NIBAN1",
  "SERPINB3",
  "MAPRE3"
)

volcano_interaction <- volcano_plot(results_interaction, labels_interaction) +
  labs(title = "Interaction effect (Type I vs. Type II)")

# Res, main
labels_res_main <- c(
  "MUSTN1",
  "ANKRD2",
  "MYBPH",
  "THY1",
  "S100A13",
  "TMSB4X",
  "CSRP3",
  "KLHL40",
  "VIM"
)

volcano_res_main <- volcano_plot(results_res_main, labels_res_main) +
  labs(title = "Main effect")

# Res, type 1
labels_res_i <- c(
  "MUSTN1",
  "ANKRD2",
  "MYBPH",
  "S100A13"
)

volcano_res_i <- volcano_plot(results_res_i, labels_res_i) +
  labs(title = "Type I")

# Res, type 2
labels_res_ii <- c(
  "MYH6",
  "KLHL40",
  "DNAJA4",
  "MUSTN1",
  "S100A13",
  "THY1",
  "CSRP3",
  "KLHL41"
)

volcano_res_ii <- volcano_plot(results_res_ii, labels_res_ii) +
  labs(title = "Type II")

# Res, interaction
labels_res_interaction <- c(
  "SERPINB3",
  "PGM5",
  "MYH6",
  "MYH4",
  "TNNC1",
  "TNNI1",
  "MYL3",
  "MYL2"
)

volcano_res_interaction <- volcano_plot(results_res_interaction, labels_res_interaction) +
  labs(title = "Interaction effect (Type I vs. Type II)")

# Ter main
labels_ter_main <- c(
  "S100A13",
  "FTL",
  "PRKAA1",
  "AAMDC",
  "SORB51"
)

volcano_ter_main <- volcano_plot(results_ter_main, labels_ter_main) +
  labs(title = "Main effect")

# Ter type I
labels_ter_i <- c(
  "S100A13",
  "UNC45B",
  "FTL",
  "ATP2B1",
  "RPS19",
  "RPL3L",
  "DSP"
)

volcano_ter_i <- volcano_plot(results_ter_i, labels_ter_i) +
  labs(title = "Type I")

# Ter type II
labels_ter_ii <- c(
  "S100A13",
  "LCP1",
  "SOD3",
  "HIBADH",
  "PRKAA1",
  "OSTF1"
)

volcano_ter_ii <- volcano_plot(results_ter_ii, labels_ter_ii) +
  labs(title = "Type II")

# Ter interaction
labels_ter_interaction <- c(
  "HOMER2",
  "COLA1A2",
  "RPL3L"
)

volcano_ter_interaction <- volcano_plot(results_ter_interaction, labels_ter_interaction) +
  labs(title = "Interaction effect (Type I vs. Type II)")

# Collect volcanoes into single layout
all_volcanoes <- (volcano_main | volcano_i | volcano_ii | volcano_interaction) /
  (volcano_res_main | volcano_res_i | volcano_res_ii | volcano_res_interaction) /
  (volcano_ter_main | volcano_ter_i | volcano_ter_ii | volcano_ter_interaction)

ggplot2::ggsave(snakemake@output[["fig_volcanoes"]], width = 180, height = 150, units = "mm", plot = all_volcanoes)



# Venn diagram -----------------------------------------------------------
cat("\033[1;33mCreating Venn diagram plot (Fig. 2)\033[0m\n")

plot_venn <- ggvenn::ggvenn(
  list(
    "ter_i" = dplyr::filter(results_ter_i, xiao < 0.05)$protein,
    "ter_ii" = dplyr::filter(results_ter_ii, xiao < 0.05)$protein,
    "res_i" = dplyr::filter(results_res_i, xiao < 0.05)$protein,
    "res_ii" = dplyr::filter(results_res_ii, xiao < 0.05)$protein
  ),
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4, show_elements = F
)

ggplot2::ggsave(snakemake@output[["fig_venn"]], width = 100, height = 100, units = "mm", plot = plot_venn)


# Top 25 regulated proteins ----------------------------------------------
# Identify the top 25 most upregulated proteins
cat("\033[1;33mCreating Top25 upregulated proteins plot (Fig. 2)\033[0m\n")

sum_most_up <- df_long_l2fc_mean %>%
  dplyr::filter(n > 4) %>%
  dplyr::group_by(protein) %>%
  dplyr::summarize(sum = sum(mean_l2fc, na.rm = T)) %>%
  dplyr::slice_max(order_by = sum, n = 25) %>%
  dplyr::pull(protein)

# Dataframe for plot
df_top25 <- df_long_l2fc_mean %>%
  dplyr::filter(protein %in% sum_most_up) %>%
  dplyr::mutate(color = paste(intervention, fiber_type, sep = "_"))

# Save csv
readr::write_csv(df_top25, snakemake@output[["top25_upregulated"]])

# Plot
plot_top25 <- ggplot2::ggplot(df_top25, aes(x = reorder(protein, mean_l2fc), y = mean_l2fc, fill = color)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, color = "black", size = 7)
  ) +
  scale_fill_manual(
    values = c(
      resistance_I = res_i_color,
      resistance_II = res_ii_color,
      terbutaline_I = ter_i_color,
      terbutaline_II = ter_ii_color
    ),
    labels = c(
      "resistance_I" = "Resistance Type I",
      "resistance_II" = "Resistance Type II",
      "terbutaline_I" = "B2A, Type I",
      "terbutaline_II" = "B2a, Type II"
    )
  ) +
  labs(y = expression(paste(Sigma, " Log2fold change")), x = "")

ggplot2::ggsave(snakemake@output[["fig_top25_upregulated"]], width = 120, height = 50, units = "mm", plot = plot_top25)


# MYH6 plot --------------------------------------------------------------
cat("\033[1;33mCreating MYH6 plot (Fig. 2)\033[0m\n")

myh6_abundances <- df_long %>%
  dplyr::filter(protein == "MYH6")

readr::write_csv(myh6_abundances, snakemake@output[["myh6_abundances"]])

plot_myh6 <- myh6_abundances %>%
  dplyr::mutate(fiber_intervention = interaction(intervention, fiber_type, sep = "_")) %>%
  ggplot(aes(x = factor(time, levels = c("pre", "post")), y = abundance, color = fiber_intervention)) +
  geom_point(size = 1) +
  geom_path(aes(group = interaction(id, fiber_type))) +
  facet_grid(~intervention) +
  scale_color_manual(
    values = c(
      resistance_I = res_i_color,
      resistance_II = res_ii_color,
      terbutaline_I = ter_i_color,
      terbutaline_II = ter_ii_color
    ),
    labels = c(
      "resistance_I" = "Resistance Type I",
      "resistance_II" = "Resistance Type II",
      "terbutaline_I" = "B2A, Type I",
      "terbutaline_II" = "B2a, Type II"
    )
  ) +
  scale_x_discrete(labels = c(pre = "Pre", post = "Post")) +
  theme(
  ) +
  labs(x = "", y = "log2 (abundance)")

ggplot2::ggsave(snakemake@output[["fig_myh6"]], width = 35, height = 40, units = "mm", plot = plot_myh6)


# Correlations - type II vs. type I --------------------------------------
cat("\033[1;33mCreating correlation plot (type II vs. type I) (Fig. 2)\033[0m\n")
# Create plotting data frame
df <- df_long_l2fc_mean %>%
  dplyr::select(protein, intervention, fiber_type, mean_l2fc) %>%
  tidyr::pivot_wider(
    names_from = fiber_type,
    values_from = mean_l2fc
  )

readr::write_csv(df, snakemake@output[["cor_i_vs_ii"]])

# Prepare labels
labels <- c(
  "MYH6",
  "MYL2",
  "MYL3",
  "TNNC1",
  "TNNI1",
  "ANXA1",
  "CTSG",
  "MYH11",
  "MYH8"
)

# Plot
plot_cor_i_vs_ii <- df %>%
  dplyr::mutate(label_col = ifelse(protein %in% labels, protein, NA)) %>%
  ggplot(aes(x = I, y = II)) +
  geom_point(alpha = 0.1, shape = 16, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.1) +
  labs(
    x = "Log2fold change, Type 1",
    y = "Log2fold change, Type 2",
    title = ""
  ) +
  geom_smooth(
    method = "lm",
    linewidth = 0.25,
    color = "black",
    se = F
  ) +
  geom_text_repel(aes(label = label_col),
    point.size = 1,
    size = 1,
    min.segment.length = 0.1,
    force = 0.3,
    segment.size = 0.1
  ) +
  ggpubr::stat_cor(
    method = "pearson",
    label.x = -Inf,
    label.y = Inf,
    hjust = -0.1,
    vjust = 1.1,
    size = 2.5
  ) +
  facet_grid(~intervention,
    labeller = as_labeller(c(resistance = "RES", terbutaline = "B2A"))
  )

# Save plot
ggplot2::ggsave(snakemake@output[["fig_cor_i_vs_ii"]], width = 100, height = 50, units = "mm", plot = plot_cor_i_vs_ii)


# Correlations - RES vs. B2A ---------------------------------------------
cat("\033[1;33mCreating correlation plot (RES vs. B2a) (Fig. 2)\033[0m\n")

# Create plotting data frame
df <- df_long_l2fc_mean %>%
  dplyr::select(protein, intervention, fiber_type, mean_l2fc) %>%
  tidyr::pivot_wider(
    names_from = intervention,
    values_from = mean_l2fc
  )

readr::write_csv(df, snakemake@output[["cor_res_vs_b2a"]])

# Prepare labels
labels <- c(
  "MYH6",
  "MYL2",
  "MYL3",
  "TNNC1",
  "TNNI1",
  "ANXA1",
  "CTSG",
  "MYH11",
  "MYH8"
)

# Plot
plot_cor_res_vs_b2a <- df %>%
  dplyr::mutate(label_col = ifelse(protein %in% labels, protein, NA)) %>%
  ggplot(aes(x = resistance, y = terbutaline)) +
  geom_point(alpha = 0.1, shape = 16, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.1) +
  labs(
    x = "Log2fold change, RES",
    y = "Log2fold change, B2A",
    title = ""
  ) +
  geom_smooth(
    method = "lm",
    linewidth = 0.25,
    color = "black",
    se = F
  ) +
  geom_text_repel(aes(label = label_col),
    point.size = 1,
    size = 1,
    min.segment.length = 0.1,
    force = 0.3,
    segment.size = 0.1
  ) +
  ggpubr::stat_cor(
    method = "pearson",
    label.x = -Inf,
    label.y = Inf,
    hjust = -0.1,
    vjust = 1.1,
    size = 2.5
  ) +
  facet_grid(~fiber_type,
    labeller = as_labeller(c(I = "Type I", II = "Type II"))
  )

# Save plot
ggplot2::ggsave(snakemake@output[["fig_cor_res_vs_b2a"]], width = 100, height = 50, units = "mm", plot = plot_cor_res_vs_b2a)


# GSEA results -----------------------------------------------------------
cat("\033[1;33mCreating GSEA plot (Fig. 2)\033[0m\n")

# Convert to data frame
df_gsea <- gsea_all@compareClusterResult

# Mark terms of interest
terms <- c(
  "cellular respiration",
  "mitochondrion organization",
  "cytoplasmic translation",
  "cytoskeleton organization",
  "muscle contraction",
  "muscle cell development",
  "striated muscle contration",
  "fatty acid beta-oxidation",
  "tricarboxylic acid cycle",
  "gluconeogenesis",
  "cell differentiation",
  "cell development",
  "translation",
  "cell population proliferation",
  "ribosome assembly",
  "actin filament bundle assembly",
  "striated muscle cell differentation",
  "ribosome biogenesis",
  "actin cytoskeleton organization",
  "muscle structure development",
  "G protein-coupled receptor signaling pathway",
  "striated muscle adaptation",
  "myofibril assembly",
  "actin filament organization"
)

# Clean up data frame
df_gsea_plot <- df_gsea %>%
  dplyr::filter(Description %in% terms) %>%
  dplyr::mutate(direction = ifelse(NES > 0, "Enriched", "Depleted")) %>%
  dplyr::mutate(gene_ratio = (str_count(core_enrichment, "/") + 1) / setSize)

# Determine plotting order
order <- df_gsea_plot %>%
  dplyr::group_by(Description) %>%
  dplyr::summarize(mean = mean(NES)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(abs(mean))

# Set order
df_gsea_plot <- df_gsea_plot %>%
  dplyr::mutate(
    Description = factor(Description, levels = (order$Description)),
    Cluster = factor(Cluster, levels = c(
      "ter_i",
      "ter_ii",
      "res_i",
      "res_ii"
    ))
  )

# Plot
plot_gsea <- ggplot(df_gsea_plot, aes(x = Cluster, y = Description, size = gene_ratio)) +
  geom_point(aes(fill = p.adjust), shape = 21, color = "black", stroke = 0.25) +
  scale_size(
    range = c(2, 10),
    name = "Gene Ratio"
  ) +
  scale_fill_gradient(
    low = "#e26664",
    high = "#327ebd",
    name = "p.adjust"
  ) +
  scale_x_discrete(labels = c(
    ter_i = "B2A\nType I",
    ter_ii = "B2A\nType II",
    res_i = "RES\nType I",
    res_ii = "RES\nType II"
  )) +
  facet_grid(~ factor(direction, levels = c("Enriched", "Depleted"))) +
  theme(
    legend.title = element_text(),
    axis.text.x = element_text(size = 7)
  ) +
  labs(y = NULL, x = NULL)

# Save plot
ggplot2::ggsave(snakemake@output[["fig_gsea"]], height = 95, width = 145, units = "mm", plot = plot_gsea)
