
#' used for creating a line of xiao values on volcano
#'
#' @param logFC
#' @param p
#'
#' @return
#' @export
#'
#' @examples
xiao_line <- function(logFC, p_threshold = 0.05) {
    -log10(10^(log10(p_threshold) / abs(logFC)))
}

#' Volcano plot
#'
#' @param results_sheet
#' @param labels
#'
#' @return
#' @export
#'
#' @examples
volcano_plot <- function(results_sheet, labels, xlim = NULL, ylim = NULL){

  #If no limits are supplied, these are standards
  if (is.null(xlim)){
    xlim <- c(-2.5, 2.5)
  }

  if (is.null(ylim)){
    ylim <- c(0, 10)
  }

  #Count upregulated proteins
  upregulated_count <- results_sheet %>%
    dplyr::filter(xiao < 0.05 & logFC > 0) %>%
    nrow()

  #Count downregulated proteins
  downregulated_count <- results_sheet %>%
    dplyr::filter(xiao < 0.05 & logFC < 0) %>%
    nrow()

  #Plot
  results_sheet %>%
    dplyr::mutate(label = case_when(
      protein %in% labels ~protein
    )) %>%
    dplyr::mutate(color = case_when(
      logFC >= 0 & xiao <= 0.05 ~ "Upregulated",
      logFC <= 0 & xiao <= 0.05 ~ "Downregulated",
      TRUE ~ "Unchanged")) %>%
    ggplot(aes(x=logFC, y=-log10(P.Value)))+
    geom_point(aes(color = color, alpha=color), size = 1, shape = 16)+
    #Add dashed Xiao line
    geom_function(fun = xiao_line, linetype = "dotted",
                  linewidth = 0.25, alpha = 0.5) +
    geom_text_repel(aes(label = label), point.size=1, size=2,
                    min.segment.length = 0.1, force=0.3,
                    segment.size = 0.1, na.rm = T, max.overlaps = 20)+
    #Add protein counts in corners
    annotate("text", x = xlim[2]-0.5, y = ylim[2], label = paste("Up:", upregulated_count),
             size = 2, hjust = 1, vjust = 1, color = upregulated_color) +
    annotate("text", x = xlim[1]+0.5, y = ylim[2], label = paste("Down:", downregulated_count),
             size = 2, hjust = 0, vjust = 1, color = downregulated_color) +
    theme(
      legend.position = "none",
      aspect.ratio = 1,
      plot.title = element_text(hjust = 0.5)
    )+
    scale_color_manual(breaks = c("Upregulated", "Downregulated", "Unchanged"),
                       values=c(upregulated_color, downregulated_color, "gray50"))+
    scale_alpha_manual(breaks = c("Upregulated", "Downregulated", "Unchanged"),
                       values=c(1, 1, 0.5))+
    labs(
      x = "l2fc (post-pre)",
      y = "-log10 (p-value)"
    )+
    coord_cartesian(y = ylim,
                    x = xlim)

}


#' Plot and summarize GO term-specific changes in proteomics data
#'
#' This function performs statistical testing and visualization of
#' log2-fold changes and abundance changes for user-specified GO terms
#' (biological process, cellular component, or molecular function).
#' It requires the objects "df_long" and "df_long_l2fc" to be loaded in
#' the environment.
#'
#' @param term_list A character vector of GO terms (e.g., "GO:0031012").
#' Each term is matched against the "gocc", "gobp", or "gomf" columns
#'  in the dataset to subset proteins.
#' @param term_labels An optional named character vector
#' mapping GO terms to readable labels for plotting.
#' @param term_labels An optional character vector setting the order along
#' the x-axis of chosen GO-terms
#'
#' @details The function performs two types of statistical tests:
#' Within-group comparison: paired t-tests comparing pre vs. post intervention,
#' stratified by intervention group and muscle fiber type.
#'
#' Between-group comparison: unpaired t-tests comparing log2-fold
#' changes between interventions.
#' The results are visualized as bar plots with error bars (95% CI),
#' grouped by GO term and fiber type.
#'
#' @examples
#'
#' @examples
#' \dontrun{
#' term_vec <- c("GO:0031012", "GO:0006936")
#' term_names <- c("GO:0031012" = "ECM", "GO:0006936" = "Muscle contraction")
#' results <- plot_terms(term_vec, term_labels = term_names)
#' results$plot  # display the ggplot
#' }
plot_terms <- function(term_list, term_labels = NULL, term_order = NULL){

  fiber_types <- c("I", "II", "all")

  #Lists to store results
  term_results_interaction <- list()
  term_results_within <- list()

  for (term in term_list){
    for (fibertype in fiber_types){

      #Create key to store unique result names
      key <- paste(term, fibertype, sep = "_")

      #Interaction data
      interaction <- df_long_l2fc %>%
        dplyr::filter(grepl(term, gocc) | grepl(term, gobp) | grepl(term, gomf)) %>%
        dplyr::mutate(goterm = term)

      if (fibertype != "all"){
        interaction <- interaction %>%
          dplyr::filter(fiber_type == fibertype) %>%
          rstatix::t_test(l2fc ~ intervention, paired=F, detailed=TRUE) %>%
          dplyr::select(estimate, p)
      }
      else {
        interaction <- interaction %>%
          rstatix::t_test(l2fc ~ intervention, paired=F, detailed=TRUE) %>%
          dplyr::select(estimate, p)
      }

      #Within group data
      within <- df_long %>%
        dplyr::filter(grepl(term, gocc) | grepl(term, gobp) | grepl(term, gomf)) %>%
        dplyr::mutate(goterm = term)

      if (fibertype != "all"){
        within <- within %>%
          dplyr::group_by(intervention) %>%
          dplyr::filter(fiber_type == fibertype) %>%
          rstatix::t_test(abundance ~ time, paired=TRUE, detailed=TRUE, ref.group = "post") %>%
          dplyr::mutate(goterm = term) %>%
          dplyr::mutate(fiber_type = fibertype)
      }
      else {
        within <- within %>%
          dplyr::group_by(intervention) %>%
          rstatix::t_test(abundance ~ time, paired=TRUE, detailed=TRUE, ref.group = "post") %>%
          dplyr::mutate(goterm = term) %>%
          dplyr::mutate(fiber_type = fibertype)
      }

      #Add to results lists
      term_results_interaction[[key]] <- interaction
      term_results_within[[key]] <- within
  }
}

  #Combine results
  within_results <- dplyr::bind_rows(term_results_within)
  interaction_results <- dplyr::bind_rows(term_results_interaction)

  #Set order for plotting
  if (!is.null(term_order)){
    within_results <- within_results %>%
      dplyr::mutate(goterm = factor(goterm, levels = term_order))
  }

  #Plot results
  plot <- ggplot(within_results, aes(x = goterm, y = estimate,
                         fill = intervention,
                         ymin = conf.low,
                         ymax = conf.high))+
    stat_summary(fun = mean, geom = "col", position = position_dodge(width = 0.75), width = 0.7)+
    geom_errorbar(width = 0, position = position_dodge(0.75), color = "gray")+
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25)+
    theme(
      legend.position = "inside"
    )+
    scale_fill_manual(values = c(terbutaline = terbutaline_color,
                                 resistance = resistance_color),
                      label = c(terbutaline = "B2A", resistance = "RES")) +
    labs(y = "Log2fold change", x = "")+
    facet_grid(~fiber_type, labeller = as_labeller(c(I = "Type I",
                                                     II = "Type II",
                                                     all = "Fiber type independent")))

  if (!is.null(term_labels)){
    plot <- plot + scale_x_discrete(labels = term_labels)
  }

  return (list(plot = plot,
               within_results = within_results,
               interaction_results = interaction_results))

}

#' Plot and summarize mitochondrial changes in proteomics data
#'
#' This function performs statistical testing and visualization of
#' log2-fold changes and abundance changes for user-specified mitochondrial
#' proteins. It requires the objects "df_long" and "df_long_l2fc" to be
#' loaded in the environment.
#'
#' @param term_list A character vector of mitochondrial categories (from mitocarta).
#' Each term is matched against the "mito" column
#' in the dataset to subset proteins.
#' @param term_labels An optional named character vector
#' mapping mitochondrial terms to readable labels for plotting.
#' @param term_labels An optional character vector setting the order along
#' the x-axis of chosen mitochondrial terms
#'
#' @details The function performs two types of statistical tests:
#' Within-group comparison: paired t-tests comparing pre vs. post intervention,
#' stratified by intervention group and muscle fiber type.
#'
#' Between-group comparison: unpaired t-tests comparing log2-fold
#' changes between interventions.
#' The results are visualized as bar plots with error bars (95% CI),
#' grouped by GO term and fiber type.
#'
#' @examples
#'
#' @examples
#' \dontrun{
#' term_vec <- c("GO:0031012", "GO:0006936")
#' term_names <- c("GO:0031012" = "ECM", "GO:0006936" = "Muscle contraction")
#' results <- plot_terms(term_vec, term_labels = term_names)
#' results$plot  # display the ggplot
#' }
plot_complexes <- function(complex_list, complex_labels = NULL, complex_order = NULL){

  fiber_types <- c("I", "II", "all")

  #Lists to store results
  complex_results_interaction <- list()
  complex_results_within <- list()

  for (complex in complex_list){
    for (fibertype in fiber_types){

      #Create key to store unique result names
      key <- paste(complex, fibertype, sep = "_")

      #Create filter string to filter df_long_l2fc
      filter_string <- paste(complex, "subunits", sep = " ")

      #Interaction data
      interaction <- df_long_l2fc %>%
        dplyr::filter(grepl(filter_string, mito)) %>%
        dplyr::mutate(complex = complex)

      if (fibertype != "all"){
        interaction <- interaction %>%
          dplyr::filter(fiber_type == fibertype) %>%
          rstatix::t_test(l2fc ~ intervention, paired=F, detailed=TRUE) %>%
          dplyr::select(estimate, p)

      }
      else {
        interaction <- interaction %>%
          rstatix::t_test(l2fc ~ intervention, paired=F, detailed=TRUE) %>%
          dplyr::select(estimate, p)
      }

      #Within group data
      within <- df_long %>%
        dplyr::filter(grepl(filter_string, mito)) %>%
        dplyr::mutate(complex = complex)

      if (fibertype != "all"){
        within <- within %>%
          dplyr::group_by(intervention) %>%
          dplyr::filter(fiber_type == fibertype) %>%
          rstatix::t_test(abundance ~ time, paired=TRUE, detailed=TRUE, ref.group = "post") %>%
          dplyr::mutate(complex = complex) %>%
          dplyr::mutate(fiber_type = fibertype)
      }
      else {
        within <- within %>%
          dplyr::group_by(intervention) %>%
          rstatix::t_test(abundance ~ time, paired=TRUE, detailed=TRUE, ref.group = "post") %>%
          dplyr::mutate(complex = complex) %>%
          dplyr::mutate(fiber_type = fibertype)
      }

      #Add to results lists
      complex_results_interaction[[key]] <- interaction
      complex_results_within[[key]] <- within
    }
  }

  #Combine results
  within_results <- dplyr::bind_rows(complex_results_within)
  interaction_results <- dplyr::bind_rows(complex_results_interaction)

  #Set order for plotting
  if (!is.null(complex_order)){
    within_results <- within_results %>%
      dplyr::mutate(complex = factor(complex, levels = complex_order))
  }

  #Plot results
  plot <- ggplot(within_results, aes(x = complex, y = estimate,
                                     fill = intervention,
                                     ymin = conf.low,
                                     ymax = conf.high))+
    stat_summary(fun = mean, geom = "col", position = position_dodge(width = 0.75), width = 0.7)+
    geom_errorbar(width = 0, position = position_dodge(0.75), color = "gray")+
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25)+
    theme(
      legend.position = "inside"
    )+
    scale_fill_manual(values = c(terbutaline = terbutaline_color,
                                 resistance = resistance_color),
                      label = c(terbutaline = "B2A", resistance = "RES")) +
    labs(y = "Log2fold change", x = "")+
    facet_grid(~fiber_type, labeller = as_labeller(c(I = "Type I",
                                                     II = "Type II",
                                                     all = "Fiber type independent")))

  if (!is.null(complex_labels)){
    plot <- plot + scale_x_discrete(labels = complex_labels)
  }

  return (list(plot = plot,
               within_results = within_results,
               interaction_results = interaction_results))

}
