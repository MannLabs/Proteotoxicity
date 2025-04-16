## -- Table S3, Stats

write_tsv(limma_border_firstrow, file = "../output/Tables/Stats_scDVP.tsv")

# Core proteome

d_long %>%
  filter(sample %in% SA_incl) -> d_filtered

d_filtered %>%
  select(Protein.Group, sample, int) %>%
  spread(sample, int) %>%
  column_to_rownames("Protein.Group") %>%
  filter(complete.cases(.)) %>%
  rownames_to_column("protein") -> core_proteome

write_tsv(core_proteome, file = "../output/Tables/scDVP_core-proteome.tsv")

# Expression matrix

d_long %>%
  select(well_id, int, Protein.Group) %>%
  spread(well_id, int) -> d_wide

write_tsv(d_wide, file = "../output/Tables/scDVP-expression-matrix.tsv")

# meta table

d_long_noexcl %>%
  select(sample, well, batch) %>%
  left_join(meta) %>%
  distinct(sample_short, sample, biology, spatial, well, batch, rank_id, area_gpt) %>%
  mutate(included = sample %in% SA_incl) %>%
  select(sample_short, biology, spatial, well, batch, rank_id, area_gpt, included) -> meta_supplement

write_tsv(meta_supplement, "../output/Tables/scDVP-meta.tsv")


