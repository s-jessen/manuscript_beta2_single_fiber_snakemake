rule all:
  input:
    #Results
    'data/results.rds',
    'data/results_i_and_ii_interaction.rds',
    'data/results_i_interaction.rds',
    'data/results_ii_interaction.rds',
    'data/results_myh_ii_vs_i.rds',

    #GSEA outputs
    'data/gsea_res_i.rds',
    'data/gsea_res_ii.rds',
    'data/gsea_ter_i.rds',
    'data/gsea_ter_ii.rds',
    'data/gsea_all.rds',
    'data/gsea_i_and_ii_interaction.rds',
    'data/gsea_i_interaction.rds',
    'data/gsea_ii_interaction.rds',

    #GSEA csv
    'data/gsea/gsea_results_res_i.csv',
    'data/gsea/gsea_results_res_ii.csv',
    'data/gsea/gsea_results_ter_i.csv',
    'data/gsea/gsea_results_ter_ii.csv',
    'data/gsea/gsea_results_all.csv',
    'data/gsea/gsea_results_i_and_ii_interaction.csv',
    'data/gsea/gsea_results_i_interaction.csv',
    'data/gsea/gsea_results_ii_interaction.csv',

    #Figure 1 outputs
    'data/myh_distribution.csv',
    'data/figures/figure1/myh_distribution.svg',
    'data/figures/figure1/venn_myh.svg',
    'data/figures/figure1/median_log10_intensity.svg',
    'data/figures/figure1/myh_volcano.svg',
    'data/figures/figure1/PCA.svg',
    'data/figures/figure1/loadings.svg',

    #Figure 2 output
    'data/top25_upregulated.csv',
    'data/myh6_abundances.csv',
    'data/correlations_i_vs_ii.csv',
    'data/correlations_res_vs_b2a.csv',
    'data/figures/figure2/volcanoes.svg',
    'data/figures/figure2/venn.svg',
    'data/figures/figure2/top25_upregulated.svg',
    'data/figures/figure2/myh6.svg',
    'data/figures/figure2/cor_i_vs_ii.svg',
    'data/figures/figure2/cor_res_vs_b2a.svg',
    'data/figures/figure2/gsea.svg',

    #Figure 3 output
    'data/figures/figure3/volcanoes.svg',
    'data/figures/figure3/gsea.svg',
    'data/figures/figure3/ribosomal_proteins.svg',
    'data/figures/figure3/mitochondrial_proteins.svg',
    'data/figures/figure3/structural_go_terms.svg',
    'data/figures/figure3/metabolic_go_terms.svg'

rule data_prep:
  input:
    settings = 'R/settings.R',
    raw_data = 'data-raw/proteomic_data.csv',
    design = 'data-raw/design.xlsx',
    keywords = 'data-raw/keywords.xlsx',
    mitocarta = 'data-raw/mitocarta.xls'
    
  output:
    metadata = 'data/metadata.rds',
    se_raw = 'data/se_raw.rds',
    se = 'data/se.rds',
    df_long = 'data/df_long.rds',
    df_long_l2fc = 'data/df_long_l2fc.rds',
    df_long_l2fc_mean = 'data/df_long_l2fc_mean.rds'
    
  script:
    "R/data_preparation.R"

rule limma:
  input: 
    se = 'data/se.rds'

  output:
    #Intervention group independent results
    results_main = 'data/results_main.rds',
    results_i = 'data/results_i.rds',
    results_ii = 'data/results_ii.rds',
    results_interaction = 'data/results_interaction.rds',
    
    #Terbutaline results
    results_ter_main = 'data/results_ter_main.rds',
    results_ter_i = 'data/results_ter_i.rds',
    results_ter_ii = 'data/results_ter_ii.rds',
    results_ter_interaction = 'data/results_ter_interaction.rds',
    
    #Resistance training results
    results_res_main = 'data/results_res_main.rds',
    results_res_i = 'data/results_res_i.rds',
    results_res_ii = 'data/results_res_ii.rds',
    results_res_interaction = 'data/results_res_interaction.rds',

    #All the above results
    results = 'data/results.rds',
    
    #Interaction between intervention groups
    results_i_and_ii_interaction = 'data/results_i_and_ii_interaction.rds',  #Independent of fiber type
    results_i_interaction = 'data/results_i_interaction.rds',                #For type I fibers
    results_ii_interaction = 'data/results_ii_interaction.rds',              #For type II fibers
    
    #Baseline differences in fiber type pools
    results_myh_ii_vs_i = 'data/results_myh_ii_vs_i.rds'
    
  script:
    "R/limma.R"

rule gsea:
  input:
    #Resistance training results
    results_res_main = 'data/results_res_main.rds',
    results_res_i = 'data/results_res_i.rds',
    results_res_ii = 'data/results_res_ii.rds',
    results_res_interaction = 'data/results_res_interaction.rds',
    
    #Terbutaline results
    results_ter_main = 'data/results_ter_main.rds',
    results_ter_i = 'data/results_ter_i.rds',
    results_ter_ii = 'data/results_ter_ii.rds',
    results_ter_interaction = 'data/results_ter_interaction.rds'

  output:
    gsea_res_i = 'data/gsea_res_i.rds',
    gsea_res_ii = 'data/gsea_res_ii.rds',
    gsea_ter_i = 'data/gsea_ter_i.rds',
    gsea_ter_ii = 'data/gsea_ter_ii.rds',
    gsea_all = 'data/gsea_all.rds',

    df_gsea_res_i = 'data/gsea/gsea_results_res_i.csv',
    df_gsea_res_ii = 'data/gsea/gsea_results_res_ii.csv',
    df_gsea_ter_i = 'data/gsea/gsea_results_ter_i.csv',
    df_gsea_ter_ii = 'data/gsea/gsea_results_ter_ii.csv',
    df_gsea_all = 'data/gsea/gsea_results_all.csv'

  script:
    "R/gsea.R"

rule gsea_interaction:
  input:
    #Interaction results
    results_i_and_ii_interaction = 'data/results_i_and_ii_interaction.rds',
    results_i_interaction = 'data/results_i_interaction.rds',
    results_ii_interaction = 'data/results_ii_interaction.rds'

  output:
    #Files
    gsea_i_and_ii_interaction = 'data/gsea_i_and_ii_interaction.rds',
    gsea_i_interaction = 'data/gsea_i_interaction.rds',
    gsea_ii_interaction = 'data/gsea_ii_interaction.rds',

    #Data frames
    df_gsea_i_and_ii_interaction = 'data/gsea/gsea_results_i_and_ii_interaction.csv',
    df_gsea_i_interaction = 'data/gsea/gsea_results_i_interaction.csv',
    df_gsea_ii_interaction = 'data/gsea/gsea_results_ii_interaction.csv'

  script:
    "R/gsea_interaction.R"

rule figure1:
  input:
    functions = 'R/functions.R',
    settings = 'R/settings.R',
    df_long = 'data/df_long.rds',
    results_myh_ii_vs_i = 'data/results_myh_ii_vs_i.rds',
    se_raw = 'data/se_raw.rds',
    se = 'data/se.rds',
    metadata = 'data/metadata.rds'
  
  output:
    #Data frames
    df_long_myh = 'data/myh_distribution.csv', #MYH distribution dataset

    #Figures
    fig_myh_distribution = 'data/figures/figure1/myh_distribution.svg',
    fig_venn = 'data/figures/figure1/venn_myh.svg',
    fig_ranked_intensities = 'data/figures/figure1/median_log10_intensity.svg',
    fig_myh_volcano = 'data/figures/figure1/myh_volcano.svg',
    fig_pca = 'data/figures/figure1/PCA.svg',
    fig_pca_loadings = 'data/figures/figure1/loadings.svg'

  script:
    'R/figure1.R'

rule figure2:
  input:
    functions = 'R/functions.R',
    settings = 'R/settings.R',
    df_long = 'data/df_long.rds',
    df_long_l2fc = 'data/df_long_l2fc.rds',
    df_long_l2fc_mean = 'data/df_long_l2fc_mean.rds',

    #Intervention group independent results
    results_main = 'data/results_main.rds',
    results_i = 'data/results_i.rds',
    results_ii = 'data/results_ii.rds',
    results_interaction = 'data/results_interaction.rds',
    
    #Terbutaline results
    results_ter_main = 'data/results_ter_main.rds',
    results_ter_i = 'data/results_ter_i.rds',
    results_ter_ii = 'data/results_ter_ii.rds',
    results_ter_interaction = 'data/results_ter_interaction.rds',
    
    #Resistance training results
    results_res_main = 'data/results_res_main.rds',
    results_res_i = 'data/results_res_i.rds',
    results_res_ii = 'data/results_res_ii.rds',
    results_res_interaction = 'data/results_res_interaction.rds',

    #GSEA results
    gsea_all = 'data/gsea_all.rds'

  output:
    #Data frames
    top25_upregulated = 'data/top25_upregulated.csv',
    myh6_abundances = 'data/myh6_abundances.csv',
    cor_i_vs_ii = 'data/correlations_i_vs_ii.csv',
    cor_res_vs_b2a = 'data/correlations_res_vs_b2a.csv',

    #Figures
    fig_volcanoes = 'data/figures/figure2/volcanoes.svg',
    fig_venn = 'data/figures/figure2/venn.svg',
    fig_top25_upregulated = 'data/figures/figure2/top25_upregulated.svg',
    fig_myh6 = 'data/figures/figure2/myh6.svg',
    fig_cor_i_vs_ii = 'data/figures/figure2/cor_i_vs_ii.svg',
    fig_cor_res_vs_b2a = 'data/figures/figure2/cor_res_vs_b2a.svg',
    fig_gsea = 'data/figures/figure2/gsea.svg'

  script:
    'R/figure2.R'

rule figure3:
  input:
    functions = 'R/functions.R',
    settings = 'R/settings.R',
    df_long = 'data/df_long.rds',
    df_long_l2fc = 'data/df_long_l2fc.rds',

    #Interaction between intervention groups
    results_i_and_ii_interaction = 'data/results_i_and_ii_interaction.rds',  
    results_i_interaction = 'data/results_i_interaction.rds',                
    results_ii_interaction = 'data/results_ii_interaction.rds',

    #GSEA interaction results
    gsea_i_and_ii_interaction = 'data/gsea_i_and_ii_interaction.rds',
    gsea_i_interaction = 'data/gsea_i_interaction.rds',
    gsea_ii_interaction = 'data/gsea_ii_interaction.rds'

  output:
    #Figures
    fig_volcanoes = 'data/figures/figure3/volcanoes.svg',
    fig_gsea = 'data/figures/figure3/gsea.svg',
    fig_ribosomal_proteins = 'data/figures/figure3/ribosomal_proteins.svg',
    fig_mitochondrial_proteins = 'data/figures/figure3/mitochondrial_proteins.svg',
    fig_structural_go_terms = 'data/figures/figure3/structural_go_terms.svg',
    fig_metabolic_go_terms = 'data/figures/figure3/metabolic_go_terms.svg'

  script:
    'R/figure3.R'