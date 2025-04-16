## -- Table S1, binned

level_of_alpha1 <- d_long %>%
  filter(Protein.Group == "P01009") %>%
  dplyr::rename(ref_int = int) %>%
  dplyr::select(well_id, ref_int)

level_of_alpha1 %>%
  mutate(ref_int = ref_int) %>%
  mutate(
    AAT_Bin = cut(
      ref_int,
      breaks = seq(min(ref_int, na.rm = TRUE),
                   max(ref_int, na.rm = TRUE),
                   length.out = 11),  # Creates 10 bins
      include.lowest = TRUE,
      labels = paste0("Bin_", 1:10))) -> level_of_alpha1
      
  d_long %>%
    left_join(meta_pg) %>%
    left_join(meta_biopsies) %>%
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
    filter(complete.cases(.)) %>%
    left_join(limma_level %>% select(Genes, adj.P.Val))-> d_fraction

  write_tsv(d_fraction, file = "../output/Tables/Binned_DVP.tsv")
  