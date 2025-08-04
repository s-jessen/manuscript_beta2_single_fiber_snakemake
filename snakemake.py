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
    "R/data_prep.R"
