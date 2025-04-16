## -- Figure S5A

limma_fstatus <- function(x, score){
  
  x %>%
    filter(Protein.Group %in% pg_30) %>%
    left_join(meta_biopsies) %>%
    filter(kleiner_score == score) %>%
    filter(well_id %in% SA_stats) %>%
    mutate(int = log2(int)) %>%
    select(well_id, int, Protein.Group) %>%
    spread(well_id, int) %>%
    column_to_rownames("Protein.Group") -> d_stats_f
  
  level <- unlist(d_stats_f["P01009",])
  
  design <- model.matrix(~ level)
  
  fit <- lmFit(d_stats_f, design)
  fit <- eBayes(fit)
  limma_level_f <- topTable(fit, number = Inf, confint = TRUE, coef = "level", adjust.method = "fdr")
  
  return(limma_level_f)
}

limma_fstatus(d_long, 1) -> limma_level_1
limma_fstatus(d_long, 4) -> limma_level_4

limma_level_f <- data.frame(Protein.Group = rownames(limma_level_1), F1_logFC = limma_level_1$logFC, F1_p.adj = limma_level_1$adj.P.Val) %>%
  left_join(limma_level_4 %>% rownames_to_column("Protein.Group") %>% select(Protein.Group, logFC, adj.P.Val) %>% rename(F4_logFC = logFC, F4_p.adj = adj.P.Val)) %>%
  left_join(meta_pg) %>%
  mutate(annotate = abs(F4_logFC) > 0.3 | abs(F1_logFC) > 0.35)
  

ggplot(data = limma_level_f, aes(x = F1_logFC, y = F4_logFC))+
  geom_point()+
  geom_density_2d(color = "grey60")+
  geom_text_repel(data = limma_level_f %>% filter(annotate == TRUE), aes(x = F1_logFC, y = F4_logFC, label = Genes))+
  geom_abline(yintercept = 0, slope = 1, lty = "dotted")+
  geom_hline(yintercept = 0, lty = "dotted") +
  geom_vline(xintercept = 0, lty = "dotted") +
  theme_classic()+
  coord_fixed()

ggsave(file = "../output/Figures/Figure_S5A.pdf", width = 8, height = 8)


  
