####################################
## Morphology-guided DVP (Fig. 4) ##
####################################

## -- Prepare Workspace
cat("\014")
rm(list=ls())
set.seed(123)

setwd(dirname(getActiveDocumentContext()$path))

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

## --  Read data
pg <- read_tsv("../data/protein_quantification/mgDVP_pg.tsv") %>%
  select(contains(c("Protein", "Gene")))

d <- read_tsv("../data/protein_quantification/mgDVP_report.tsv")

meta <- read_tsv("../data/meta/meta_mgDVP_well-assigment.txt") %>%
  filter(!slide %in% c("Explant_DVP_02", "Explant_DVP_04"))

## -- Call scripts
source("Figure_S8A.R")
source("Figure_4D.R")
