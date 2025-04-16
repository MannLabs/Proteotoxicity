## -- Figure 2F

d_long %>%
  left_join(meta_pg) %>%
  left_join(meta_biopsies) %>%
  filter(Protein.Group %in% pg_70) %>%
  mutate(ENSEMBL = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,  
                          keys=Genes, 
                          column="ENSEMBL", 
                          keytype="SYMBOL",
                          multiVals="first")) %>%
  left_join(level_of_alpha1) %>%
  drop_na(ref_int) %>%
  filter(ENSEMBL != "NANA") -> d_ancova

d_ancova %>%
  dplyr::select(Protein.Group, Genes, int, ref_int, kleiner_score) %>%
  filter(log2(ref_int) < 25) %>%
  drop_na(Genes) -> d_slope_low

kleiner_slope_low <- c()
counter <- 0
for(i in unique(d_slope_low$Genes)){
  
  #print(counter)
  
  d_slope_low %>%
    filter(Genes == i) -> d_ancova_select
  
  for(j in c(1:4)){
    d_ancova_select %>%
      filter(kleiner_score == j) -> tmp
    
    if(dim(tmp)[1] < 5) next
    
    lm_model <- lm(log2(int) ~ log2(ref_int), data = tmp)
    summary(lm_model) -> lm_summary
    
    slope <- lm_summary[["coefficients"]][2]
    
    result <- data.frame(Genes = i, Stage = j, slope = slope, df = lm_summary[["df"]][2])
    
    kleiner_slope_low <- rbind(kleiner_slope_low, result)
  }
  counter = counter + 1
}

kleiner_slope_low %>%
  left_join(meta_pg) %>%
  mutate(ENSEMBL = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,  
                          keys=Genes, 
                          column="ENSEMBL", 
                          keytype="SYMBOL",
                          multiVals="first")) -> kleiner_slope_low


kegg_slope_test <- c()
for(i in names(db_kegg$kg.sets)){
  
  #print(i)
  db_kegg_specific <- data.frame(ENTREZ = str_replace(db_kegg$kg.sets[[i]], ".*_", "")) %>%
    mutate(ENSEMBL = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,  
                            keys=ENTREZ, 
                            column="ENSEMBL", 
                            keytype="ENTREZID",
                            multiVals="first"))
  
  kleiner_slope_low %>%
    filter(ENSEMBL %in% db_kegg_specific$ENSEMBL)  -> tmp
  
  if(dim(tmp)[1] < 10) next
  
  test_result <- wilcox.test(tmp %>% filter(Stage %in% 1) %>% pull(slope),
                             tmp %>% filter(Stage %in% 4) %>% pull(slope))
  
  kegg_slope_test <- rbind(kegg_slope_test, data.frame(pathway = i, p.value = test_result$p.value))
}

kegg_slope_test %>%
  mutate(p_adjusted = p.adjust(kegg_slope_test$p.value, method = "fdr")) %>%
  arrange(p_adjusted) %>%
  mutate(rank = c(1:length(kegg_slope_test$pathway))) -> kegg_slope_test

ggplot(kegg_slope_test, aes(x = rank, y = -log10(p_adjusted))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), lty = "dotted") +
  theme_classic()+
  geom_text_repel(data = kegg_slope_test %>% filter(p_adjusted < 0.01), aes(x = rank, y = -log10(p_adjusted), label = pathway))

ggsave(file = "../output/Figures/Figure_2F.pdf", width = 5, height = 5)