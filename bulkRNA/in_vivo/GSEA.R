library(tidyverse, quietly = TRUE)
library(org.Hs.eg.db, quietly = TRUE)
library(clusterProfiler, quietly = TRUE)
library(enrichplot, quietly = TRUE)

deg <- read_csv("/datapool/pengguoyu/ATACseq/20200916ATAC/KO_vs_WT/DEG/AllGenes.csv") %>% 
  dplyr::filter(!is.na(entrezgene_id)) %>% 
  dplyr::distinct(entrezgene_id, .keep_all=TRUE) %>% 
  dplyr::arrange(desc(log2FoldChange))
rankList <- deg$log2FoldChange
names(rankList) <- deg$entrezgene_id
head(rankList)
tail(rankList)

KEGG <-  gseKEGG(geneList = rankList, organism = "hsa", keyType = "ncbi-geneid", nPerm = 5000, pvalueCutoff = 1)
GO_BP <- gseGO(geneList = rankList, ont = "BP", OrgDb = org.Hs.eg.db, keyType = "ENTREZID", nPerm = 5000, pvalueCutoff = 1)
saveRDS(KEGG, file = "/datapool/pengguoyu/ATACseq/20200916ATAC/KO_vs_WT/GSEA/KEGG.rds")
saveRDS(GO_BP, file = "/datapool/pengguoyu/ATACseq/20200916ATAC/KO_vs_WT/GSEA/GO-BP.rds")
as_tibble(KEGG@result) %>% write_csv("/datapool/pengguoyu/ATACseq/20200916ATAC/KO_vs_WT/GSEA/KEGG.csv")
as_tibble(GO_BP@result) %>% write_csv("/datapool/pengguoyu/ATACseq/20200916ATAC/KO_vs_WT/GSEA/GO-BP.csv")

pdf("/datapool/pengguoyu/ATACseq/20200916ATAC/KO_vs_WT/GSEA/KEGG_ES1.pdf")
for (i in 1:nrow(KEGG)) {
  print(gseaplot(x = KEGG, geneSetID = i, title = KEGG$Description[i]))
}
dev.off()

pdf("/datapool/pengguoyu/ATACseq/20200916ATAC/KO_vs_WT/GSEA/KEGG_ES2.pdf")
for (i in 1:nrow(KEGG)) {
  print(gseaplot2(x = KEGG, geneSetID = i, title = KEGG$Description[i]))
}
dev.off()
