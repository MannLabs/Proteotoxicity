## -- Figure R2

er_proteins <- c("EIF2AK3", "ATF4", "ERN1", "ATF6", "HSPA5", "HSP90B1", "DNAJB11", "HYOU1",
                 "P4HB", "PDIA3", "PDIA4", "PDIA6", "ERO1A", "ERO1B", "ERO1L", "CALR", "CANX")

# Load the FASTA file
fasta_file <- "../data/meta/ER-proteins__uniprotkb_organism_id_9606_AND_cc_scl_t_2024_12_14.fasta"  # Replace with the path to your FASTA file
fasta_sequences <- readAAStringSet(fasta_file)

# Extract headers
headers <- names(fasta_sequences)

# Extract UniProt IDs and gene names using regex
fasta_info <- data.frame(
  UniProt_ID = sub("^(\\S+).*", "\\1", headers),  # Extract first word in header (assumes UniProt ID)
  Gene_Name = sub(".* GN=(\\S+).*", "\\1", headers)  # Extract Gene Name (GN=...)
)

# Handle cases where gene names are missing
fasta_info$Gene_Name[!grepl("GN=", headers)] <- NA

fasta_info %>%
  distinct(Gene_Name) %>%
  mutate(UPR = Gene_Name %in% er_proteins) -> fasta_subset

limma_level %>%
  filter(Genes %in% fasta_subset$Gene_Name) %>%
  mutate(UPR = Genes %in% er_proteins) -> limma_er

ggplot() +
  geom_boxplot(data = limma_er, aes(x = UPR, y = logFC), outlier.shape = NA)+
  geom_jitter(data = limma_er %>% filter(group == "not significant"), aes(x = UPR, y = logFC), color = "grey80", width = 0.1, alpha = 0.5) +
  geom_jitter(data = limma_er %>% filter(group != "not significant"), aes(x = UPR, y = logFC), color = "black", width = 0.1, alpha = 0.5)+
  geom_text_repel(data = limma_er %>% filter(UPR == TRUE), aes(x = UPR, y = logFC, label = Genes), max.overlaps = 20)+
  theme_classic()

ggsave(file = "../output/Figures/Figure_R2.pdf", width = 5, height = 8)
