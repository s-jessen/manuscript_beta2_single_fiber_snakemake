rule data_prep:
  input:
    raw_data = 'data-raw/proteomic_data.csv',
    design = 'data-raw/design.xlsx',
    keywords = 'data-raw/keywords.xlsx',
    mitocarta = 'data-raw/mitocarta.xls'
    
  output:
    metadata = 'data/metadata.rda',
    se_raw = 'data/se_raw.rda',
    se = 'data/se.rda',
    df_long = 'data/df_long.rda',
    df_long_l2fc = 'data/df_long_l2fc.rda',
    df_long_l2fc_mean = 'data/df_long_l2fc_mean.rda'
    
  script:
    "R/data_preparation.R"

rule limma:
  input: 
    se = 'data/se.rda'

  output:
    #Intervention group independent results
    results_main = 'data/results_main.rda',
    results_i = 'data/results_i.rda',
    results_ii = 'data/results_ii.rda',
    results_interaction = 'data/results_interaction.rda',
    
    #Terbutaline results
    results_ter_main = 'data/results_ter_main.rda',
    results_ter_i = 'data/results_ter_i.rda',
    results_ter_ii = 'data/results_ter_ii.rda',
    results_ter_interaction = 'data/results_ter_interaction.rda',
    
    #Resistance training results
    results_res_main = 'data/results_res_main.rda',
    results_res_i = 'data/results_res_i.rda',
    results_res_ii = 'data/results_res_ii.rda',
    results_res_interaction = 'data/results_res_interaction.rda',

    #All the above results
    results = 'data/results.rda',
    
    #Interaction between intervention groups
    results_i_and_ii_interaction = 'data/results_i_and_ii_interaction.rda',  #Independent of fiber type
    results_i_interaction = 'data/results_i_interaction.rda',                #For type I fibers
    results_ii_interaction = 'data/results_ii_interaction.rda',              #For type II fibers
    
    #Baseline differences in fiber type pools
    results_myh_ii_vs_i = 'data/results_myh_ii_vs_i.rda'
    
  script:
    "R/limma.R"

include: "rules/renv.smk"