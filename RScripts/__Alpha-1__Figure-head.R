################################################################################
# Rosenberger et al. #
# The proteomic landscape of proteotoxic stress in a fibrogenic liver disease #
################################################################################

# This code allows to reproduce the figures of the paper. For a custom report
# file about a particular protein or pathway of interest, please use
# Alpha-1-report-generator.R

## -- Prepare Workspace
cat("\014")
rm(list=ls())
set.seed(123)

## -- Package setup
packages <- c("rstudioapi",
              "tidyverse",
              "limma",
              "viridis",
              "ggrepel",
              "pheatmap")

for (p in packages) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p)
    library(p, character.only = TRUE)
  }
}

## -- Set Working Directory
setwd(dirname(getActiveDocumentContext()$path))

lapply(c("../output/Figures", "../output/Variables/", "../output/Tables"), function(paths){
  if (!dir.exists(paths)) {
    # create the folder if it doesn't exist
    dir.create(paths)
    cat("Folder created: ", paths, "\n")
  } else {
    cat("Folder already exists: ", paths, "\n")
  }
})

invisible(file.remove(list.files("../output/Figures", full.names = TRUE)))
invisible(file.remove(list.files("../output/Variables", full.names = TRUE)))
invisible(file.remove(list.files("../output/Tables", full.names = TRUE)))

source("_Alpha-1__Figure-biopsies.R")
source("_Alpha-1__Figure-scDVP.R")
source("_Alpha-1__Figure-mgDVP.R")
