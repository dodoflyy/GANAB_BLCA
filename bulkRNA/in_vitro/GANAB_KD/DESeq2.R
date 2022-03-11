library(tidyverse)
library(tximport)
library(DESeq2)

mainDir <- file.path("/home/pengguoyu/Workshop/MultiOmics/GANAB2020/20210813UMUC3_KD/RNAseq")
exprDir <- file.path(mainDir, "Expression")
degDir <- file.path(mainDir, "DEGs")
SalmonDir <- file.path(mainDir, "Salmon")
sampleGPath <- file.path(mainDir, "SampleGroup.csv")
biomartPath <- file.path("/home/pengguoyu/Database/GRCh38/ENSEMBL/Release-101/Biomart/GRCh38.p13_Release101_20201109.csv")
txiPath <- file.path("/home/pengguoyu/Database/GRCh38/GENCODE/gencode.v35.annotation.TxToGene.csv")

geneMap <- read_csv(biomartPath) %>% 
  dplyr::select(`Gene stable ID`, `NCBI gene (formerly Entrezgene) ID`, `HGNC symbol`) %>% 
  dplyr::rename(ensembl_gene_id = `Gene stable ID`, entrezgene_id = `NCBI gene (formerly Entrezgene) ID`, hgnc_symbol = `HGNC symbol`) %>% 
  dplyr::filter(!(is.na(entrezgene_id) & is.na(hgnc_symbol))) %>% 
  dplyr::arrange(entrezgene_id, desc(hgnc_symbol)) %>% dplyr::distinct(ensembl_gene_id, .keep_all=TRUE)

glimpse(geneMap)

txToGene <- read.csv(txiPath, header = TRUE, quote = "", stringsAsFactors = FALSE)
head(txToGene)

sampleList <- c("CTR_1", "CTR_2", "CTR_3", "KD_1", "KD_2", "KD_3")
fileList <- file.path(SalmonDir, sampleList, "quant.sf")
names(fileList) <- sampleList
print(fileList)

txi <- tximport(fileList, type = "salmon", tx2gene = txToGene)
tpm1 <- txi$abundance
tpm2 <- tpm1[rowSums(tpm1) > 0, ]
tpm3 <- as_tibble(tpm2, rownames = "gene_id") %>% 
  tidyr::separate(col = gene_id, into = c("ensembl_gene_id", "ensembl_gene_version"), sep="\\.", remove = TRUE, extra = "drop", fill = "right") %>% 
  dplyr::select(-ensembl_gene_version) %>% dplyr::left_join(geneMap, by="ensembl_gene_id") %>% 
  dplyr::select(ensembl_gene_id, entrezgene_id, hgnc_symbol, everything())

tpmPath <- file.path(exprDir, "TPM.csv")
write_csv(tpm3, tpmPath)


# DEGs
sampleGroup <- read.csv(sampleGPath, header = TRUE, row.names = 1, quote = "", stringsAsFactors = TRUE)
print(sampleGroup)

dds1 <- DESeqDataSetFromTximport(txi, colData = sampleGroup, design = ~ Condition)
keep <- rowSums(counts(dds1) >= 5) >= 2
dds2 <- dds1[keep, ]
dds3 <- DESeq(dds2)
print(resultsNames(dds3))

readCounts <- DESeq2::counts(dds3) %>% as_tibble(rownames = "gene_id") %>% 
  tidyr::separate(col = gene_id, into = c("ensembl_gene_id", "ensembl_gene_version"), sep="\\.", remove = TRUE, extra = "drop", fill = "right") %>% 
  dplyr::select(-ensembl_gene_version) %>% dplyr::left_join(geneMap, by="ensembl_gene_id") %>% 
  dplyr::select(ensembl_gene_id, entrezgene_id, hgnc_symbol, everything())

normalizedCounts <- DESeq2::counts(dds3, normalized = TRUE) %>% as_tibble(rownames = "gene_id") %>% 
  tidyr::separate(col = gene_id, into = c("ensembl_gene_id", "ensembl_gene_version"), sep="\\.", remove = TRUE, extra = "drop", fill = "right") %>% 
  dplyr::select(-ensembl_gene_version) %>% dplyr::left_join(geneMap, by="ensembl_gene_id") %>% 
  dplyr::select(ensembl_gene_id, entrezgene_id, hgnc_symbol, everything())

rLog1 <- rlog(dds3, blind = TRUE) %>% assay() %>% as_tibble(rownames = "gene_id") %>% 
  tidyr::separate(col = gene_id, into = c("ensembl_gene_id", "ensembl_gene_version"), sep="\\.", remove = TRUE, extra = "drop", fill = "right") %>% 
  dplyr::select(-ensembl_gene_version) %>% dplyr::left_join(geneMap, by="ensembl_gene_id") %>% 
  dplyr::select(ensembl_gene_id, entrezgene_id, hgnc_symbol, everything())

rLog2 <- rlog(dds3, blind = FALSE) %>% assay() %>% as_tibble(rownames = "gene_id") %>% 
  tidyr::separate(col = gene_id, into = c("ensembl_gene_id", "ensembl_gene_version"), sep="\\.", remove = TRUE, extra = "drop", fill = "right") %>% 
  dplyr::select(-ensembl_gene_version) %>% dplyr::left_join(geneMap, by="ensembl_gene_id") %>% 
  dplyr::select(ensembl_gene_id, entrezgene_id, hgnc_symbol, everything())

exprPath1 <- file.path(exprDir, "ReadCounts.csv")
exprPath2 <- file.path(exprDir, "NormalizedCounts.csv")
exprPath3 <- file.path(exprDir, "RLogBlind.csv")
exprPath4 <- file.path(exprDir, "RLog.csv")
write_csv(readCounts, exprPath1)
write_csv(normalizedCounts, exprPath2)
write_csv(rLog1, exprPath3)
write_csv(rLog2, exprPath4)

rawRes <- results(dds3, contrast = c("Condition", "KD", "CTR"))
lfcRes <- lfcShrink(dds3, coef = "Condition_KD_vs_CTR", type = "apeglm", res = rawRes)
summary(rawRes)
maPath <- file.path(exprDir, "MA_Plot.pdf")
pdf(maPath, width = 9)
plotMA(rawRes, ylim = c(-6, 6), main="RAW")
plotMA(lfcRes, ylim = c(-6, 6), main="SHRINK")
dev.off()

rawDEGs1 <- as_tibble(rawRes, rownames = "gene_id") %>% 
  tidyr::separate(col = gene_id, into = c("ensembl_gene_id", "ensembl_gene_version"), sep="\\.", remove = TRUE, extra = "drop", fill = "right") %>% 
  dplyr::select(-ensembl_gene_version) %>% dplyr::left_join(geneMap, by="ensembl_gene_id") %>%
  dplyr::select(ensembl_gene_id, entrezgene_id, hgnc_symbol, everything())
rawDEGs2 <- filter(rawDEGs1, !is.na(padj), abs(log2FoldChange) >= 1, padj < 0.05)

lfcDEGs1 <- as_tibble(lfcRes, rownames = "gene_id") %>% 
  tidyr::separate(col = gene_id, into = c("ensembl_gene_id", "ensembl_gene_version"), sep="\\.", remove = TRUE, extra = "drop", fill = "right") %>% 
  dplyr::select(-ensembl_gene_version) %>% dplyr::left_join(geneMap, by="ensembl_gene_id") %>%
  dplyr::select(ensembl_gene_id, entrezgene_id, hgnc_symbol, everything())
lfcDEGs2 <- filter(lfcDEGs1, !is.na(padj), abs(log2FoldChange) >= 1, padj < 0.05)

rawDEGPath1 <- file.path(degDir, "KD_Raw_All.csv")
rawDEGPath2 <- file.path(degDir, "KD_Raw_Filter.csv")
write_csv(rawDEGs1, rawDEGPath1)
write_csv(rawDEGs2, rawDEGPath2)

lfcDEGPath1 <- file.path(degDir, "KD_Shrink_All.csv")
lfcDEGPath2 <- file.path(degDir, "KD_Shrink_Filter.csv")
write_csv(lfcDEGs1, lfcDEGPath1)
write_csv(lfcDEGs2, lfcDEGPath2)

imagePath <- file.path(exprDir, "DESeq2.RData")
save.image(imagePath)