## Proteomics Analysis: β2-Agonist and Resistance Training (single skeletal muscle fibers)

This repository contains the code and data to reproduce the analysis from an unpublished study on β2-agonist treatment combined with resistance training and their effects on the skeletal muscle proteome in
type I and IIa muscle fibers, respectively.

## Project Setup

This project uses [`renv`](https://rstudio.github.io/renv/) to manage R package dependencies in a fully self-contained way.

##Installation Steps

1. **Clone the repository**
```bash
git clone https://github.com/s-jessen/manuscript_beta2_single_fiber.git
cd manuscript_beta2_single_fiber
```

2. **Install packages and set up environment**

You can either run the setup file in a terminal:
```bash
Rscript setup.R
```

Or open setup.R in RStudio and run it interactively.
  
## 📁 Project Folder Structure

This project follows a reproducible research structure using the `prodigenr` framework. Below is an overview of the folder layout and purpose of key files:

```
├── data                                    #Stores output of scripts (for easy retrieval without need to rerun code)
│   ├── figures                             #Folder to store figures
│   └── ...                                 #Output files from analyses (GSEA results, statistical outputs, etc.)
├── data-raw                                #Raw data
│   ├── design.xlsx                         #Metadata file
│   ├── keywords.xlsx                       #Uniprot keywords file (for gene annotation)
│   ├── mitocarta.xls                       #MitoCarta 3.0 (for gene annotation)
│   └── proteomic_data.csv                  #Proteomic data
├── DESCRIPTION                             
├── docs                                    #Manuscript figures
├── manuscript_beta2_single_fiber.Rproj
├── R                                       #Main analysis scripts
│   ├── 1_data_preparation.qmd              #Cleans up raw data and stores SummarizedExperiment
│   ├── 2_gsea.qmd                          #Runs GSEA.
│   ├── 2_limma.qmd                         #Runs limma analyses.
│   ├── figure1.qmd                         #Code associated with Figure 1
│   ├── figure2.qmd                         #Code associated with Figure 2
│   ├── figure3.qmd                         #Code associated with Figure 3
│   ├── functions.R
│   ├── gsea_results.R                      #Loads GSEA results from /data
│   ├── libraries.R                         #Loads required libraries
│   ├── results.R                           #Loads limma results from /data
│   └── settings.R                          #Contains color coding, ggplot setup, etc.
├── README.md                               #This file
├── renv                                    #Renv for reproducible environment
├── renv.lock                               #Renv file for storing version information of packages
└── setup.R                                 #Script file to install required packages for reproducible analysis
```

The R project contains code to reproduce the following highlighted figures:

![Figure 1](docs/Figure1.png)
![Figure 2](docs/Figure2.png)
![Figure 3](docs/Figure3.png)
![Figure 4](docs/Figure4.png)
![Figure 5](docs/Figure5.png)
![Figure 6](docs/Figure6.png)
