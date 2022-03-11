library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

mainDir <- file.path("..")
degDir <- file.path(mainDir, "DEGs")
enrichDir <- file.path(mainDir, "Enrichment")
degPath <- file.path(degDir, "DEGsAll.csv")

degData <- read_csv(degPath) %>% 
  dplyr::filter(!is.na(entrezgene_id))
glimpse(degData)
bgEntrez <- degData$entrezgene_id %>% unique() %>% as.character()
degEntrez <- dplyr::filter(degData, abs(log2FoldChange) >= 1 & padj < 0.05) %>% 
  dplyr::pull(entrezgene_id) %>% unique()

length(bgEntrez)
length(degEntrez)

gobpTest1 <- enrichGO(gene = degEntrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", 
                      ont = "BP", readable = TRUE, universe = bgEntrez, 
                      minGSSize = 15, maxGSSize = 250)
gobpTest2 <- clusterProfiler::simplify(gobpTest1)
gobpResult <- gobpTest2@result %>% as_tibble()
glimpse(gobpResult)
gobpPath <- file.path(enrichDir, "GO-BP.csv")
write_csv(gobpResult, gobpPath)

keggTest <- enrichKEGG(gene = degEntrez, organism = "hsa", keyType = "ncbi-geneid", 
                       pvalueCutoff = 0.1, universe = bgEntrez, 
                       minGSSize = 15, maxGSSize = 250)
keggResult <- keggTest@result %>% as_tibble()
glimpse(keggResult)
keggPath <- file.path(enrichDir, "KEGG.csv")
write_csv(keggResult, keggPath)

keyKegg <- c("hsa04020", "hsa04010", "hsa04151", "hsa04066", "hsa04657",
             "hsa04668", "hsa04062", "hsa04024", "hsa04390", "hsa04933",
             "hsa04550", "hsa03320", "hsa04068")

# keyKegg <- c("hsa04020", "hsa04024")
keggTest2 <- enrichKEGG(gene = degEntrez, organism = "hsa", keyType = "ncbi-geneid", 
                       pvalueCutoff = 1, universe = bgEntrez, 
                       minGSSize = 15, maxGSSize = 250)
keggTest3 <- keggTest2[keyKegg, asis = TRUE]
cnetplot(keggTest3, showCategory = 13, max.overlaps = 100)
imagePath <- file.path(enrichDir, "Enrichment.RData")
save.image(imagePath)