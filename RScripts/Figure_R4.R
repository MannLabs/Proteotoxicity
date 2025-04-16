## -- Figure R2

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

## -- Read limma from bulk, and CORUM

d <- read_tsv("../output/Tables/Stats_biopsies.tsv")
corum <- read_tsv("../data/meta/protein_complexes/corum_humanComplexes.txt") %>%
  separate_rows(subunits_uniprot_id, sep = ";")

d %>%
  mutate(in_complex = Protein.Group %in% corum$subunits_uniprot_id) %>%
  mutate(label = abs(logFC) > 2) -> d_complex

ggplot()+
  geom_boxplot(data = d_complex, aes(x = in_complex, y = logFC))+
  geom_text_repel(data = d_complex %>% filter(label == TRUE), aes(x = in_complex, y = logFC, label = Genes)) +
  theme_bw()

ggsave(file = "../output/Figures/Figure_R2-CORUM.pdf", width = 3, height = 5)

## -- Use Hu.MAP 2.0

humap <- read_csv("../data/meta/protein_complexes/humap2_complexes_20200809.txt") %>%
  separate_rows(Uniprot_ACCs, sep = " ")

d %>%
  mutate(in_complex = Protein.Group %in% humap$Uniprot_ACCs) %>%
  mutate(label = abs(logFC) > 2) -> d_complex

ggplot()+
  geom_boxplot(data = d_complex, aes(x = in_complex, y = logFC))+
  geom_text_repel(data = d_complex %>% filter(label == TRUE), aes(x = in_complex, y = logFC, label = Genes))+
  theme_bw()

ggsave(file = "../output/Figures/Figure_R2-HUMAP.pdf", width = 3, height = 5)
