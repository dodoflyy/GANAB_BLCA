library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)

mainDir <- file.path("..")
degDir <- file.path(mainDir, "DEGs")
enrichDir <- file.path(mainDir, "Enrichment")
degPath <- file.path(degDir, "AllGenes.csv")

degData1 <- read_csv(degPath) %>% 
  dplyr::filter(!is.na(entrezgene_id))
degData2 <- filter(degData1, abs(log2FoldChange) >= 1 & padj < 0.05)
glimpse(degData2)

geneList1 <- degData1$entrezgene_id %>% unique() %>% as.character()
geneList2 <- degData2$entrezgene_id %>% unique()

gobpTest1 <- enrichGO(gene = geneList2, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", 
                      ont = "BP", universe = geneList1, minGSSize = 15, maxGSSize = 250, readable = TRUE)
gobpTest2 <- simplify(gobpTest1)
gobpResult <- gobpTest2@result %>% as_tibble()
gobpPath <- file.path(enrichDir, "GO-BP.csv")
write_csv(gobpResult, gobpPath)

keggTest <- enrichKEGG(gene = geneList2, organism = "hsa", keyType = "ncbi-geneid", pvalueCutoff = 0.1, 
                       universe = geneList1, minGSSize = 15, maxGSSize = 250)
keggResult <- keggTest@result %>% as_tibble()
keggPath <- file.path(enrichDir, "KEGG.csv")
write_csv(keggResult, keggPath)

imagePath <- file.path(enrichDir, "Enrichment.RData")
save.image(imagePath)
