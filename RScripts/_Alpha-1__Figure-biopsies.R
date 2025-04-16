##################################
## Biopsies (Fig. 1 and Fig. 2) ##
##################################

# Note that the code for clinical summary statistics is not documented here
# due to data protection reasons.

cat("\014")
rm(list=ls())
set.seed(123)

setwd(dirname(getActiveDocumentContext()$path))

## -- Load data
d_biopsies_noexcl <- read_tsv("../data/protein_quantification/Biopsies_PiZ_report.tsv") %>%
  select(contains(c( "Protein.Group", "Evo12")))

meta_biopsies <- read_tsv("../data/meta/meta_biopsies.txt")

meta_pg <- read_tsv("../data/protein_quantification/Biopsies_PiZ_report.tsv") %>%
  select(!contains("Evo12"))

## -- Data filtering
d_biopsies_noexcl %>%
  gather(ms_id, int, contains("Evo12")) %>%
  mutate(int = as.numeric(int)) %>%
  filter(int != 0) -> d_long_noexcl

d_long_noexcl %>%
  group_by(ms_id) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  summarise(mean = mean(n),  sd = sd(n)) -> stats

stats_lower_n <- stats$mean - 1.5 * stats$sd

too_low_n <- d_long_noexcl %>%
  group_by(ms_id) %>%
  summarise(n = n()) %>%
  filter(n > stats_lower_n) %>%
  mutate(well_id = str_replace(str_replace(ms_id, ".*__", ""), "_.*", "")) %>%
  filter(well_id %in% meta_biopsies$well_id) %>%
  left_join(meta_biopsies) %>%
  pull(ms_id)

SA_incl <- d_long_noexcl %>%
  group_by(ms_id) %>%
  summarise(n = n()) %>%
  filter(n > stats_lower_n) %>%
  mutate(well_id = str_replace(str_replace(ms_id, ".*__", ""), "_.*", "")) %>%
  filter(well_id %in% meta_biopsies$well_id) %>%
  left_join(meta_biopsies) %>%
  filter(tech_rep %in% c(1,NA)) %>%
  filter(include == TRUE) %>%
  pull(ms_id)

d_long_noexcl %>%
  filter(ms_id %in% SA_incl) %>%
  mutate(well_id = str_replace(str_replace(ms_id, ".*__", ""), "_.*", "")) -> d_long

# Fraction of excluded samples
meta_biopsies %>%
  filter(tech_rep %in% c(1,NA)) %>%
  filter(include == TRUE) %>%
  mutate(QC_pass = well_id %in% d_long$well_id) %>%
  group_by(QC_pass) %>%
  summarise(n = n())

meta_biopsies %>%
  left_join(d_long %>% distinct (well_id, ms_id)) %>%
  filter(ms_id %in% SA_incl) -> meta_biopsies

## -- Call figures

# Figure 1 and associated extended data figures, tables
source("Figure_S1B.R")
source("Figure_S1C.R")
source("Figure_S1D.R")
source("Figure_S1E.R")
source("Figure_S1F.R")
source("Figure_S1G.R")
#source("Figure_S1H.R")
source("Figure_1C.R")
source("Figure_1D.R")
source("Figure_1E.R")
source("Figure_S2D.R")
source("Figure_S2E.R")
source("Figure_R2.R")
source("Figure_S3A.R")
source("Figure_S3B.R")
source("Table_S1-binned.R")

# Figure 2 and associated extended data figures
source("Figure_2A.R")
source("Figure_2B.R")
source("Figure_2C.R")
source("Figure_2D.R")
source("Figure_2E.R")
source("Figure_2F.R")
source("Figure_2G.R")
source("Figure_S3C.R")
source("Figure_S3DEF.R")
source("Figure_S4.R")
source("Figure_S5A.R")
source("Figure_S5B.R")
source("Figure_S5C.R")
source("Table_S2")






