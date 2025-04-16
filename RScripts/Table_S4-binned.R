## -- Table S4, binned

level_of_alpha1 <- d_long %>%
  filter(protein == "P01009") %>%
  dplyr::rename(ref_int = int) %>%
  dplyr::select(well, ref_int)

level_of_alpha1 %>%
  mutate(
    AAT_Bin = cut(
      ref_int,
      breaks = seq(min(ref_int, na.rm = TRUE),
                   max(ref_int, na.rm = TRUE),
                   length.out = 11),  # Creates 10 bins
      include.lowest = TRUE,
      labels = paste0("Bin_", 1:10))) -> level_of_alpha1

d_long %>%
  rename(Protein.Group = protein) %>%
  left_join(pg) %>%
  left_join(meta) %>%
  left_join(level_of_alpha1) %>%
  select(Genes, sample, int, AAT_Bin, ref_int) -> d_binned

d_binned %>%
  group_by(Genes, AAT_Bin) %>%
  summarise(median_intensity = median(int, na.rm = TRUE), .groups = "drop") %>%
  group_by(Genes) %>%
  mutate(Fraction = median_intensity / sum(median_intensity, na.rm = TRUE)) %>%  # Normalize to sum = 1 per protein
  ungroup() %>%
  select(Genes, AAT_Bin, Fraction) %>%
  spread(AAT_Bin, Fraction) %>%
  filter(complete.cases(.)) -> d_fraction

write_tsv(d_fraction, file = "../output/Tables/Binned_mgDVP.tsv")

## -- Table S4, binned by CRP

level_of_alpha1 <- d_long %>%
  filter(protein == "P02741") %>%
  dplyr::rename(ref_int = int) %>%
  dplyr::select(well, ref_int)

level_of_alpha1 %>%
  mutate(
    AAT_Bin = cut(
      ref_int,
      breaks = seq(min(ref_int, na.rm = TRUE),
                   max(ref_int, na.rm = TRUE),
                   length.out = 11),  # Creates 10 bins
      include.lowest = TRUE,
      labels = paste0("Bin_", 1:10))) -> level_of_alpha1

d_long %>%
  rename(Protein.Group = protein) %>%
  left_join(pg) %>%
  left_join(meta) %>%
  left_join(level_of_alpha1) %>%
  select(Genes, sample, int, AAT_Bin, ref_int) -> d_binned

d_binned %>%
  group_by(Genes, AAT_Bin) %>%
  summarise(median_intensity = median(int, na.rm = TRUE), .groups = "drop") %>%
  group_by(Genes) %>%
  mutate(Fraction = median_intensity / sum(median_intensity, na.rm = TRUE)) %>%  # Normalize to sum = 1 per protein
  ungroup() %>%
  select(Genes, AAT_Bin, Fraction) %>%
  spread(AAT_Bin, Fraction) %>%
  filter(complete.cases(.)) -> d_fraction

write_tsv(d_fraction, file = "../output/Tables/Binned_mgDVP-CRP.tsv")

# Expression matrix

d_long %>%
  mutate(id = paste(slide, well, sep = "__")) %>%
  select(id, int, protein) %>%
  spread(id, int) -> d_wide

write_tsv(d_wide, file = "../output/Tables/mgDVP-expression.tsv")

