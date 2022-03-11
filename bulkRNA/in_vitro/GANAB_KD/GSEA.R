library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

mainDir <- file.path("..")
degDir <- file.path(mainDir, "DEGs")
gseaDir <- file.path(mainDir, "GSEA")
degPath <- file.path(degDir, "DEGsAll.csv")

degData <- read_csv(degPath) %>% dplyr::filter(!is.na(entrezgene_id)) %>% 
  dplyr::arrange(desc(hgnc_symbol)) %>% 
  dplyr::distinct(entrezgene_id, .keep_all = TRUE) %>% 
  dplyr::arrange(desc(log2FoldChange))
glimpse(degData)
gseaData <- degData$log2FoldChange
names(gseaData) <- degData$entrezgene_id
head(gseaData)

gobpTest1 <- gseGO(geneList = gseaData, ont = "BP", OrgDb = org.Hs.eg.db)
gobpTest2 <- clusterProfiler::simplify(gobpTest1)
gobpResult <- gobpTest2@result %>% as_tibble()
glimpse(gobpResult)
gobpPath1 <- file.path(gseaDir, "GO-BP.csv")
write_csv(gobpResult, gobpPath1)

gobpPath2 <- file.path(gseaDir, "GO-BP_ES.pdf")
resultSets <- gobpResult$ID
resultNames <- gobpResult$Description
pdf(gobpPath2, width = 9, height = 9)
for (i in 1:nrow(gobpResult)) {
  print(enrichplot::gseaplot2(gobpTest2, geneSetID = resultSets[i], title = resultNames[i]))
}
dev.off()

# KEGG
keggTest <- gseKEGG(geneList = gseaData, organism = "hsa", keyType = "ncbi-geneid", 
                    pvalueCutoff = 0.1, eps = 0)
keggResult <- keggTest@result %>% as_tibble()
glimpse(keggResult)
keggPath1 <- file.path(gseaDir, "KEGG.csv")
write_csv(keggResult, keggPath1)

keggPath2 <- file.path(gseaDir, "KEGG_ES.pdf")
keggSets <- keggResult$ID
keggNames <- keggResult$Description
pdf(keggPath2, width = 9, height = 9)
for (i in 1:nrow(keggResult)) {
  print(enrichplot::gseaplot2(keggTest, geneSetID = keggSets[i], title = keggNames[i]))
}
dev.off()

imagePath <- file.path(gseaDir, "GSEA.RData")
save.image(imagePath)