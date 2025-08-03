#Set mirrors
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install renv if needed
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

# Ensure BiocManager is available
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Restore all project packages
renv::restore()
