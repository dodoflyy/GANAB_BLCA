library(tidyverse)


exprs_dir <- file.path("/home/pengguoyu/Workshop/MultiOmics/GANAB2020/20220119UMUC3_PKD1_GANAB/expression")
exprs_path <- file.path(exprs_dir, "RLog.csv")

exprs_data <- read_csv(exprs_path) %>% 
  select(- entrezgene_id, - hgnc_symbol) %>% 
  as.data.frame()
rownames(exprs_data) <- exprs_data$ensembl_gene_id
exprs_data$ensembl_gene_id <- NULL

pca_data <- t(as.matrix(exprs_data))
pca <- prcomp(pca_data)
summary(pca)

new_data <- pca$x[, c(1, 2)] %>% 
  as_tibble(new_data, rownames = "Sample") %>% 
  mutate(Treatment = c(rep_len("CTR", 3), rep_len("GANAB KD", 3), rep_len("PKD1 KD", 3)))
pca_plot <- ggplot(new_data) +
  geom_point(mapping = aes(x = PC1, y = PC2, color = Treatment))
pca_plot