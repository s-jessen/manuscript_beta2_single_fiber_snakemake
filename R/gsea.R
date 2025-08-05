source("renv/activate.R")
source(snakemake@input[["settings"]])

#Load libraries
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)

## Enrichment analysis for each fiber type and intervention
gsea_results <- list()
gsea_results_filtered <- list()

results_list <- list("res_i" = readRDS(snakemake@input[["results_res_i"]]),
                     "res_ii" = readRDS(snakemake@input[["results_res_ii"]]),
                    "ter_i" = readRDS(snakemake@input[["results_ter_i"]]),
                     "ter_ii" = readRDS(snakemake@input[["results_ter_ii"]])
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
    ont = "BP",
    pvalueCutoff = 0.05,
    keyType = "SYMBOL",
    eps=0,
    minGSSize=10
    )
  
  #Save to gsea_results
  gsea_results[[comparison]] <- gsea
}

#Extract and save to data folder
gsea_res_i <- gsea_results[["res_i"]]
gsea_res_ii <- gsea_results[["res_ii"]]
gsea_ter_i <- gsea_results[["ter_i"]]
gsea_ter_ii <- gsea_results[["ter_ii"]]

saveRDS(gsea_res_i, snakemake@output[["gsea_res_i"]])
saveRDS(gsea_res_ii, snakemake@output[["gsea_res_ii"]])
saveRDS(gsea_ter_i, snakemake@output[["gsea_ter_i"]])
saveRDS(gsea_ter_ii, snakemake@output[["gsea_ter_ii"]])

for (res in names(gsea_results)){

  gsea = gsea_results[[res]]

  readr::write_csv(as.data.frame(gsea), snakemake@output[[paste0("df_gsea_", res)]])

}

#Merge results
gsea_all <- clusterProfiler::merge_result(list("res_i" = gsea_res_i,
                              "res_ii" = gsea_res_ii,
                              "ter_i" = gsea_ter_i,
                              "ter_ii" = gsea_ter_ii))

#Save merged results
saveRDS(gsea_all, snakemake@output[["gsea_all"]])
readr::write_csv(as.data.frame(gsea_all), snakemake@output[["df_gsea_all"]])