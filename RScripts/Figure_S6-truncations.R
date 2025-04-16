## -- Figure S6 truncations

## -- Load necessary libraries
library()
library(stringr)
library(tidyverse)

# In silico digest
tryptic_digest <- function(sequence) {
  peptides <- unlist(strsplit(sequence, "(?<=[KR])", perl = TRUE))
  peptides <- peptides[peptides != ""]
  return(peptides)
}

# Read FASTA file into a named list of sequences
fasta_file <- readAAStringSet("../data/meta/HsSERPINA1.fasta")
protein_sequences <- as.character(fasta_file)
protein_ids <- names(fasta_file)

# Initialize an empty list to store results
digested_proteins <- list()

# Loop through each protein sequence to perform tryptic digestion
for (i in seq_along(protein_sequences)) {
  # Perform tryptic digest
  peptides <- tryptic_digest(protein_sequences[i])
  
  # Create a data frame with numbered peptides and calculate peptide length
  peptides_numbered <- data.frame(
    ProteinID = protein_ids[i],
    PeptideNumber = seq_along(peptides),
    Peptide = peptides,
    PeptideLength = nchar(peptides)
  )
  
  # Filter peptides by length
  peptides_numbered <- peptides_numbered[peptides_numbered$PeptideLength >= 7 & peptides_numbered$PeptideLength <= 30, ]
  
  # Append to the list
  digested_proteins[[protein_ids[i]]] <- peptides_numbered
}

# Combine all results into a single data frame
digested_proteins_df <- do.call(rbind, digested_proteins)
rownames(digested_proteins_df) <- NULL  # Remove row names

## -- Load peptides after directLFQ

pre <- read_tsv('../data/protein_quantification/scDVP_precursor-report.tsv')

pre %>%
  filter(protein == "P01009") %>%
  mutate(Peptide = str_replace(ion, "\\d$", "")) %>%
  filter(Peptide %in% digested_proteins_df$Peptide) %>%
  gather(sample, int, contains("OA")) %>%
  filter(int != 0) %>%
  left_join(peptides_numbered %>% select(Peptide, PeptideNumber)) %>%
  mutate(batch = str_extract(sample, "(SoSt|FAR)")) %>%
  mutate(well = ifelse(batch == "SoSt",
                       sapply(str_extract(sample, "(?<=cW-)[^_]+(?=_Box)"), transform_string),
                       str_replace_all(sample, ".*_", ""))) %>%
  mutate(well_id = paste(batch, well, sep = "_")) %>%
  left_join(meta) -> pre_filtered

pre_filtered %>%
  filter(biology == "Border", spatial %in% c("minusone", "minustwo")) %>%
  group_by(PeptideNumber) %>%
  summarise(n = n()) -> pre_filtered_stats

ggplot(data = pre_filtered, aes(x = as.factor(PeptideNumber), y = log2(int), fill = spatial))+
  geom_boxplot()

pre_filtered %>%
  filter(biology == "Border", PeptideNumber %in% c(3,39)) %>%
  drop_na(spatial) %>%
  group_by(ion, PeptideNumber, spatial) %>%
  summarise(n = n()) %>%
  spread(spatial, n) -> pre_filtered_N

pre_filtered %>%
  filter(biology == "Border", PeptideNumber %in% c(3,39)) %>%
  drop_na(spatial) %>%
  group_by(ion, PeptideNumber, spatial) %>%
  summarise(int_median = median(int)) -> pre_filtered_int

pre_filtered %>%
  filter(biology == "Border", PeptideNumber %in% c(3,39)) %>%
  drop_na(spatial) %>%
  mutate(
    spatial = str_replace(spatial, "minus.*", "AAT-minus"),
    spatial = str_replace(spatial, "plus.*", "AAT-plus")
  ) %>%
  filter(str_detect(ion, "2$")) %>%
  select(spatial, PeptideNumber, int, sample_short, well) %>%
  mutate(PeptideNumber = as.factor(PeptideNumber), int = log2(int)) %>%
  spread(PeptideNumber, int) %>%
  mutate(ratio = `3` - `39`) -> truncation_stats


ggplot(truncation_stats, aes(x = spatial, y = ratio))+
    geom_boxplot()+
    facet_wrap(.~ sample_short)+
    theme_bw()

ggsave(file = "../output/Figures/Reviewer_Fig-R1.pdf", width = 8, height = 5)

truncation_stats %>%
  group_by(sample_short) %>%
  filter(!is.na(ratio)) %>% # Remove rows with NA in ratio
  summarise(
    t_test = list(t.test(ratio ~ spatial))
  ) -> test



