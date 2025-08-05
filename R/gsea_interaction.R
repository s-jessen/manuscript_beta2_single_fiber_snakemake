source("renv/activate.R")
source(snakemake@input[["settings"]])

#Load libraries
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)

#GSEAs on intervention x time comparison
#For fiber type independent, type 1, and type 2 comparisons
gsea_results <- list()
gsea_results_filtered <- list()

results_list <- list(
  "i_and_ii_interaction" = readRDS(snakemake@input[["results_i_and_ii_interaction"]]),
  "i_interaction" = readRDS(snakemake@input[["results_i_interaction"]]),
  "ii_interaction" = readRDS(snakemake@input[["results_ii_interaction"]])
)

for (comparison in names(results_list)){
  
  #Extract result
  res <- results_list[[comparison]]
  
  #Preparation of ranked protein list
  gsea_list <- as.numeric(res$logFC)
  names(gsea_list) = as.character(res$protein)
  gsea_list <- gsea_list[!is.na(gsea_list)]
  
  #GSEA analysis (GO:BP)
  gsea <- clusterProfiler::gseGO(
    geneList = gsea_list,
    OrgDb = org.Hs.eg.db,
    ont = "ALL",
    pvalueCutoff = 0.05,
    keyType = "SYMBOL",
    eps=0,
    minGSSize=10
    )
  
  #Save to gsea_results
  gsea_results[[comparison]] <- gsea
}

#Extract and save to data folder
gsea_i_and_ii_interaction <- gsea_results[["i_and_ii_interaction"]]
gsea_i_interaction <- gsea_results[["i_interaction"]]
gsea_ii_interaction <- gsea_results[["ii_interaction"]]

saveRDS(gsea_i_and_ii_interaction, snakemake@output[["gsea_i_and_ii_interaction"]])
saveRDS(gsea_i_interaction, snakemake@output[["gsea_i_interaction"]])
saveRDS(gsea_ii_interaction, snakemake@output[["gsea_ii_interaction"]])

for (res in names(gsea_results)){

  gsea = gsea_results[[res]]

  readr::write_csv(as.data.frame(gsea), snakemake@output[[paste0("df_gsea_", res)]])

}