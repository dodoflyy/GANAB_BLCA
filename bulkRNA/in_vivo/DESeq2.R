library(tidyverse, quietly = TRUE)
library(GenomicFeatures, quietly = TRUE)
library(tximport, quietly = TRUE)
library(DESeq2, quietly = TRUE)

SalmonDir="../Salmon"

txd <- loadDb("../GRCh38/Gencode_V33_TxDb.sqlite")
txToGene <- AnnotationDbi::select(txd, keys=keys(txd, "TXNAME"), keytype="TXNAME", columns=c("TXNAME", "GENEID"))
write_csv(txToGene, "../GRCh38/Gencode_V33_txTogene.csv")
sampleGroup <- read.csv("../SampleGroup.csv", header = TRUE, row.names = 1)
sampleGroup$Group <- factor(sampleGroup$Group, levels = c("WT", "KO"))
print(sampleGroup)

# mapping gene id
geneMap <- read_csv("../GRCh38/GRCh38.p13_V100_20200629.csv") %>% 
  dplyr::select(`Gene stable ID`, `HGNC symbol`, `NCBI gene (formerly Entrezgene) ID`) %>% 
  dplyr::rename(ensembl_gene_id=`Gene stable ID`, hgnc_symbol=`HGNC symbol`, entrezgene_id=`NCBI gene (formerly Entrezgene) ID`) %>% 
  dplyr::filter(!(is.na(entrezgene_id) & is.na(hgnc_symbol))) %>% 
  dplyr::arrange(entrezgene_id, desc(hgnc_symbol)) %>% 
  dplyr::distinct(ensembl_gene_id, .keep_all=TRUE)
glimpse(geneMap)

sampleList <- c("WT1", "WT2_2", "WT3", "KO1", "KO2", "KO3")
fileList <- file.path(SalmonDir, sampleList, "quant.sf")
names(fileList) <- c("WT1", "WT2", "WT3", "KO1", "KO2", "KO3")
txi <- tximport(fileList, type = "salmon", tx2gene = txToGene)
# get TPM matrix
txi$abundance %>% 
  as_tibble(rownames = "gene_id") %>% 
  tidyr::separate(col = gene_id, into = c("ensembl_gene_id", "ensembl_gene_version"), sep = "\\.", remove = TRUE, extra = "drop", fill = "right") %>% 
  dplyr::select(- ensembl_gene_version) %>% dplyr::left_join(geneMap, by="ensembl_gene_id") %>% 
  dplyr::select(ensembl_gene_id, entrezgene_id, hgnc_symbol, everything()) %>% 
  write_csv("../KO_vs_WT/Expression/TPM.csv")

dds <- DESeqDataSetFromTximport(txi, colData = sampleGroup, design = ~ Group)
saveRDS(dds, file = "../KO_vs_WT/DEG/DESeq2Dds.rds")

# filter genes
keep <- rowSums(counts(dds) >= 5) >= 3
dds <- dds[keep,]
dds <- DESeq(dds)

resultsNames(dds)
readCounts <- counts(dds) %>% as_tibble(rownames="gene_id") %>% 
  tidyr::separate(col = gene_id, into = c("ensembl_gene_id", "ensembl_gene_version"), sep = "\\.", remove = TRUE, extra = "drop", fill = "right") %>% 
  dplyr::select(- ensembl_gene_version) %>% dplyr::left_join(geneMap, by="ensembl_gene_id") %>% 
  dplyr::select(ensembl_gene_id, entrezgene_id, hgnc_symbol, everything())
norCounts <- counts(dds, normalized=TRUE) %>% as_tibble(rownames="gene_id") %>% 
  tidyr::separate(col = gene_id, into = c("ensembl_gene_id", "ensembl_gene_version"), sep = "\\.", remove = TRUE, extra = "drop", fill = "right") %>% 
  dplyr::select(- ensembl_gene_version) %>% dplyr::left_join(geneMap, by="ensembl_gene_id") %>% 
  dplyr::select(ensembl_gene_id, entrezgene_id, hgnc_symbol, everything())

rlogB <- rlog(dds, blind = TRUE)
rlogN <- rlog(dds, blind = FALSE)
rlogBT <- assay(rlogB) %>% as_tibble(rownames="gene_id") %>% 
  tidyr::separate(col = gene_id, into = c("ensembl_gene_id", "ensembl_gene_version"), sep = "\\.", remove = TRUE, extra = "drop", fill = "right") %>% 
  dplyr::select(- ensembl_gene_version) %>% dplyr::left_join(geneMap, by="ensembl_gene_id") %>% 
  dplyr::select(ensembl_gene_id, entrezgene_id, hgnc_symbol, everything())
rlogNT1 <- assay(rlogN) %>% as_tibble(rownames="gene_id") %>% 
  tidyr::separate(col = gene_id, into = c("ensembl_gene_id", "ensembl_gene_version"), sep = "\\.", remove = TRUE, extra = "drop", fill = "right") %>% 
  dplyr::select(- ensembl_gene_version) %>% dplyr::left_join(geneMap, by="ensembl_gene_id") %>% 
  dplyr::select(ensembl_gene_id, entrezgene_id, hgnc_symbol, everything())

write_csv(readCounts, path="../KO_vs_WT/Expression/ReadCounts.csv")
write_csv(norCounts, path="../KO_vs_WT/Expression/NormalizedCounts.csv")
write_csv(rlogBT, path="../KO_vs_WT/Expression/RlogBlind.csv")
write_csv(rlogNT1, path="../KO_vs_WT/Expression/RlogNotBlind.csv")

res <- lfcShrink(dds, coef = "Group_KO_vs_WT", type = "apeglm")
result <- as_tibble(res, rownames="gene_id") %>% 
  tidyr::separate(col = gene_id, into = c("ensembl_gene_id", "ensembl_gene_version"), sep = "\\.", remove = TRUE, extra = "drop", fill = "right") %>% 
  dplyr::select(- ensembl_gene_version) %>% dplyr::left_join(geneMap, by="ensembl_gene_id") %>% 
  dplyr::select(ensembl_gene_id, entrezgene_id, hgnc_symbol, everything())
resFilter <- dplyr::filter(result, abs(log2FoldChange) >= 1 & padj < 0.05)
sigFilter <- dplyr::filter(result, padj < 0.05)

write_csv(result, path="../KO_vs_WT/DEG/AllGenes.csv")
write_csv(resFilter, path="../KO_vs_WT/DEG/FilteredDEGs.csv")
write_csv(sigFilter, path="../KO_vs_WT/DEG/DEGs.csv")