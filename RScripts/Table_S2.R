## -- Table S2

pg_n_detected %>%
  filter(n > 0.7*length(SA_incl)) %>%
  pull(Protein.Group) -> pg_70

limma_comparison <- function(level1, level2, coef){
  SA_stats <- meta_biopsies %>%
    filter(level %in% c(level1, level2)) %>%
    pull(well_id)
  
  d_long %>%
    filter(Protein.Group %in% pg_70) %>%
    filter(well_id %in% SA_stats) %>%
    mutate(int = log2(int)) %>%
    dplyr::select(- ms_id) %>%
    spread(well_id, int) %>%
    column_to_rownames("Protein.Group") -> d_stats
  
  patient <- as.factor((meta_biopsies %>% column_to_rownames("well_id"))[colnames(d_stats),]$sample)
  level <- as.factor((meta_biopsies %>% column_to_rownames("well_id"))[colnames(d_stats),]$level)
  
  design <- model.matrix(~ level + patient)
  
  fit <- lmFit(d_stats, design)
  fit <- eBayes(fit)
  limma <- topTable(fit, number = Inf, confint = TRUE, coef = coef, adjust.method = "fdr") %>%
    rownames_to_column("Protein.Group") %>%
    left_join(meta_pg) %>%
    mutate(comparison = coef)
  
  return(limma)
}

limma_comparison(level1 = "0_low", level2 = "1_moderate", coef = "level1_moderate") -> limma_0_1
limma_comparison(level1 = "1_moderate", level2 = "2_high", coef = "level2_high") -> limma_1_2

pg_significant_extremes <- unique(c(limma_level %>%
                                      filter(adj.P.Val < 0.05) %>%
                                      pull(Protein.Group),
                                    limma_0_1 %>%
                                      filter(adj.P.Val < 0.05) %>%
                                      pull(Protein.Group),
                                    limma_1_2 %>%
                                      filter(adj.P.Val < 0.05) %>%
                                      pull(Protein.Group)))

# Plot both categories against each other
limma_0_1 %>%
  dplyr::select(Protein.Group, logFC, adj.P.Val, B) %>%
  dplyr::rename(low_moderate_fc = logFC, low_moderate_p = adj.P.Val, low_moderate_B = B) -> limma_low_moderate

limma_1_2 %>%
  dplyr::select(Protein.Group, logFC, adj.P.Val, B) %>%
  dplyr::rename(moderate_high_fc = logFC, moderate_high_p = adj.P.Val, moderate_high_B = B) -> limma_moderate_high

full_join(limma_low_moderate, limma_moderate_high) %>%
  mutate(diff_fc = low_moderate_fc - moderate_high_fc) %>%
  mutate(diff_B = low_moderate_B - moderate_high_B) %>%
  left_join(meta_pg) %>%
  mutate(annotate = abs(diff_B) > 5) %>%
  mutate(significant = Protein.Group %in% pg_significant_extremes) -> limma_staging

poi_early_up <- limma_staging %>%
  filter(low_moderate_fc > 0) %>%
  filter(diff_B > 0) %>%
  filter(max(abs(c(low_moderate_fc, moderate_high_fc))) > 0) %>%
  filter(moderate_high_p < 0.05 | low_moderate_p < 0.05)%>%
  arrange(diff_fc) 

poi_late_up <- limma_staging %>%
  filter(moderate_high_fc > 0) %>%
  filter(diff_B < 0) %>%
  filter(max(abs(c(low_moderate_fc, moderate_high_fc))) > 0) %>%
  filter(moderate_high_p < 0.05 | low_moderate_p < 0.05) %>%
  arrange(-diff_fc)

poi_early_down <- limma_staging %>%
  filter(low_moderate_fc < 0) %>%
  filter(diff_B > 0) %>%
  filter(if_else(abs(low_moderate_fc) > abs(moderate_high_fc), low_moderate_fc, moderate_high_fc) < 0) %>%
  filter(moderate_high_p < 0.05 | low_moderate_p < 0.05) %>%
  arrange(diff_fc)

poi_late_down <- limma_staging %>%
  filter(moderate_high_fc < 0) %>%
  filter(diff_B < 0) %>%
  filter(if_else(abs(low_moderate_fc) > abs(moderate_high_fc), moderate_high_fc, low_moderate_fc) < 0) %>%
  filter(moderate_high_p < 0.05 | low_moderate_p < 0.05) %>%
  arrange(-diff_fc)

# Write tables
write_tsv(poi_early_down, file = "../output/Tables/Timing_early_down.tsv")
write_tsv(poi_early_up, file = "../output/Tables/Timing_early_up.tsv")
write_tsv(poi_late_down, file = "../output/Tables/Timing_late_down.tsv")
write_tsv(poi_late_up, file = "../output/Tables/Timing_late_up.tsv")

# How many are late?
(nrow(poi_late_up)) / (nrow(poi_early_up) + nrow(poi_late_up))
