## -- Figure 3C

# Checking neighbors first and second row, A1-positive

d_norm %>%
  filter(Protein.Group %in% pg_30) %>%
  filter(biology == "Border", spatial %in% c("plusone", "plustwo")) %>%
  mutate(well = paste(batch, well, sep = "_")) %>%
  dplyr::select(well, Protein.Group, int_norm) %>%
  spread(well, int_norm) %>%
  column_to_rownames("Protein.Group") -> d_borderstats

condition <- as.factor((meta %>%
                          mutate(well = paste(batch, well, sep = "_")) %>%
                          column_to_rownames("well"))[colnames(d_borderstats),]$spatial)

sample <- as.factor((meta %>%
                       mutate(well = paste(batch, well, sep = "_")) %>%
                       column_to_rownames("well"))[colnames(d_borderstats),]$sample_short)

design <- model.matrix(~ condition + sample)

fit <- limma::lmFit(d_borderstats, design)
fit <- limma::eBayes(fit)

limma_border_pos12 <- limma::topTable(fit, number = Inf, confint = TRUE, coef = 2, adjust.method = "fdr") %>%
  rownames_to_column("Protein.Group") %>%
  left_join(meta_pg)

# Checking neighbors first and second row, A1-negative

d_norm %>%
  filter(Protein.Group %in% pg_30) %>%
  filter(biology == "Border", spatial %in% c("minusone", "minustwo")) %>%
  mutate(well = paste(batch, well, sep = "_")) %>%
  dplyr::select(well, Protein.Group, int_norm) %>%
  spread(well, int_norm) %>%
  column_to_rownames("Protein.Group") -> d_borderstats

condition <- as.factor((meta %>%
                          mutate(well = paste(batch, well, sep = "_")) %>%
                          column_to_rownames("well"))[colnames(d_borderstats),]$spatial)

sample <- as.factor((meta %>%
                       mutate(well = paste(batch, well, sep = "_")) %>%
                       column_to_rownames("well"))[colnames(d_borderstats),]$sample_short)

design <- model.matrix(~ condition + sample)

fit <- limma::lmFit(d_borderstats, design)
fit <- limma::eBayes(fit)

limma_border_neg12 <- limma::topTable(fit, number = Inf, confint = TRUE, coef = 2, adjust.method = "fdr") %>%
  rownames_to_column("Protein.Group") %>%
  left_join(meta_pg)

# Plot p-values

data.frame(comparison = c(rep("firstrow", length(limma_border_firstrow$adj.P.Val)),
                          rep("neg12", length(limma_border_neg12$adj.P.Val)),
                          rep("pos12", length(limma_border_pos12$adj.P.Val))),
           pvalue = c(limma_border_firstrow$adj.P.Val,
                      limma_border_neg12$adj.P.Val,
                      limma_border_pos12$adj.P.Val)) -> p_hist

ggplot(p_hist, aes(x = pvalue))+
  geom_histogram(bins = 50)+
  facet_grid(comparison~., scales = "free")+
  theme_classic()

ggsave("../output/Figures/Figure_3C.pdf", width = 4, height = 5)