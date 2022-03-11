library(tximport)
library(DESeq2)
library(biomaRt)
library(tidyverse)

main_dir <- file.path("/home/pengguoyu/Workshop/MultiOmics/GANAB2020/20220119UMUC3_PKD1_GANAB")
ref_dir <- file.path("/home/pengguoyu/Database/GRCh38/GENCODE")
salmon_dir <- file.path(main_dir, "salmon")
exprs_dir <- file.path(main_dir, "expression")
degs_dir <- file.path(main_dir, "DEGs")

tx2gene_path <- file.path(ref_dir, "gencode.v35.annotation.TxToGene.csv")
group_path <- file.path(main_dir, "sample_group.csv")

group_info <- read.csv(group_path, header = TRUE, row.names = 1)
sample_list <- rownames(group_info)
quant_list <- list.files(salmon_dir, pattern = "quant.sf", full.names = TRUE, 
                         recursive = TRUE)
names(quant_list) <- sample_list

tx_gene <- read.csv(tx2gene_path, header = TRUE)
txp <- tximport(files = quant_list, type = "salmon", tx2gene = tx_gene)
tpm1 <- txp$abundance
tpm2 <- tpm1[rowSums(tpm1) > 0,]

short_id <- function(long_id) {
  id_version <- strsplit(long_id, split = ".", fixed = TRUE) %>% unlist()
  ensembl_id <- id_version[1]
  return(ensembl_id)
}

ensemblgene <- sapply(rownames(tpm2), short_id, USE.NAMES = FALSE)
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(mart = ensembl, attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"), 
                  values = ensemblgene, filters = "ensembl_gene_id") %>% 
  as_tibble() %>% 
  arrange(entrezgene_id, desc(hgnc_symbol)) %>% 
  distinct(ensembl_gene_id, .keep_all = TRUE)

tpm3 <- as_tibble(tpm2, rownames = "gene_id") %>% 
  separate(col = gene_id, into = c("ensembl_gene_id", "version"), sep = "\\.") %>% 
  select(- version) %>% 
  left_join(gene_map, by = "ensembl_gene_id") %>% 
  select(ensembl_gene_id, entrezgene_id, hgnc_symbol, everything())
tpm_path <- file.path(exprs_dir, "TPM.csv")
write_csv(tpm3, tpm_path)


dds1 <- DESeqDataSetFromTximport(txi = txp, colData = group_info, design = ~ kd_gene)
keep_genes <- rowSums(counts(dds1) >= 5) >= 2
dds2 <- dds1[keep_genes,]
dds3 <- DESeq(dds2)
resultsNames(dds3)

dds_assay <- assays(dds3)
cooks <- dds_assay$cooks
cooks_path <- file.path(exprs_dir, "cooks.pdf")
pdf(cooks_path, width = 9)
boxplot(log10(cooks), main = "log10(cooks)")
dev.off()

raw_counts <- counts(dds3, normalized = FALSE) %>% 
  as_tibble(rownames = "gene_id") %>% 
  separate(col = gene_id, into = c("ensembl_gene_id", "version"), sep = "\\.") %>% 
  select(- version) %>% 
  left_join(gene_map, by = "ensembl_gene_id") %>% 
  select(ensembl_gene_id, entrezgene_id, hgnc_symbol, everything())

norm_counts <- counts(dds3, normalized = TRUE) %>% 
  as_tibble(rownames = "gene_id") %>% 
  separate(col = gene_id, into = c("ensembl_gene_id", "version"), sep = "\\.") %>% 
  select(- version) %>% 
  left_join(gene_map, by = "ensembl_gene_id") %>% 
  select(ensembl_gene_id, entrezgene_id, hgnc_symbol, everything())

rlog_design <- assay(rlog(dds3, blind = FALSE)) %>% 
  as_tibble(rownames = "gene_id") %>% 
  separate(col = gene_id, into = c("ensembl_gene_id", "version"), sep = "\\.") %>% 
  select(- version) %>% 
  left_join(gene_map, by = "ensembl_gene_id") %>% 
  select(ensembl_gene_id, entrezgene_id, hgnc_symbol, everything())

rlog_blind <- assay(rlog(dds3, blind = TRUE)) %>% 
  as_tibble(rownames = "gene_id") %>% 
  separate(col = gene_id, into = c("ensembl_gene_id", "version"), sep = "\\.") %>% 
  select(- version) %>% 
  left_join(gene_map, by = "ensembl_gene_id") %>% 
  select(ensembl_gene_id, entrezgene_id, hgnc_symbol, everything())

exprs_path1 <- file.path(exprs_dir, "read_counts.csv")
write_csv(raw_counts, exprs_path1)
exprs_path2 <- file.path(exprs_dir, "normalized_counts.csv")
write_csv(norm_counts, exprs_path2)
exprs_path3 <- file.path(exprs_dir, "RLog.csv")
write_csv(rlog_design, exprs_path3)
exprs_path4 <- file.path(exprs_dir, "RLog_blind.csv")
write_csv(rlog_blind, exprs_path4)

# GANAB DEGs
raw_res1 <- results(dds3, contrast = c("kd_gene", "GANAB", "CTR"))
shrink_res1 <- lfcShrink(dds3, coef = "kd_gene_GANAB_vs_CTR", 
                         type = "apeglm", res = raw_res1)
ma_path <- file.path(exprs_dir, "GANAB_MA.pdf")
pdf(ma_path, width = 9)
plotMA(raw_res1, ylim = c(-5, 5), main = "Raw")
plotMA(shrink_res1, ylim = c(-5, 5), main = "Shrink")
dev.off()

raw_degs1 <- as_tibble(raw_res1, rownames = "gene_id") %>% 
  separate(col = gene_id, into = c("ensembl_gene_id", "version"), sep = "\\.") %>% 
  select(- version) %>% 
  left_join(gene_map, by = "ensembl_gene_id") %>% 
  select(ensembl_gene_id, entrezgene_id, hgnc_symbol, everything())
degs_path1 <- file.path(degs_dir, "GANAB_RawAll.csv")
write_csv(raw_degs1, degs_path1)

shrink_degs1 <- as_tibble(shrink_res1, rownames = "gene_id") %>% 
  separate(col = gene_id, into = c("ensembl_gene_id", "version"), sep = "\\.") %>% 
  select(- version) %>% 
  left_join(gene_map, by = "ensembl_gene_id") %>% 
  select(ensembl_gene_id, entrezgene_id, hgnc_symbol, everything())
degs_path2 <- file.path(degs_dir, "GANAB_ShrinkAll.csv")
write_csv(shrink_degs1, degs_path2)

shrink_degs2 <- filter(shrink_degs1, !is.na(padj), padj < 0.05)
degs_path3 <- file.path(degs_dir, "GANAB_ShrinkSig.csv")
write_csv(shrink_degs2, degs_path3)
shrink_degs3 <- filter(shrink_degs2, abs(log2FoldChange) >= 1)
degs_path4 <- file.path(degs_dir, "GANAB_ShrinkPass.csv")
write_csv(shrink_degs3, degs_path4)

# PKD1
raw_res1 <- results(dds3, contrast = c("kd_gene", "PKD1", "CTR"))
shrink_res1 <- lfcShrink(dds3, coef = "kd_gene_PKD1_vs_CTR", 
                         type = "apeglm", res = raw_res1)
ma_path <- file.path(exprs_dir, "PKD1_MA.pdf")
pdf(ma_path, width = 9)
plotMA(raw_res1, ylim = c(-5, 5), main = "Raw")
plotMA(shrink_res1, ylim = c(-5, 5), main = "Shrink")
dev.off()

raw_degs1 <- as_tibble(raw_res1, rownames = "gene_id") %>% 
  separate(col = gene_id, into = c("ensembl_gene_id", "version"), sep = "\\.") %>% 
  select(- version) %>% 
  left_join(gene_map, by = "ensembl_gene_id") %>% 
  select(ensembl_gene_id, entrezgene_id, hgnc_symbol, everything())
degs_path1 <- file.path(degs_dir, "PKD1_RawAll.csv")
write_csv(raw_degs1, degs_path1)

shrink_degs1 <- as_tibble(shrink_res1, rownames = "gene_id") %>% 
  separate(col = gene_id, into = c("ensembl_gene_id", "version"), sep = "\\.") %>% 
  select(- version) %>% 
  left_join(gene_map, by = "ensembl_gene_id") %>% 
  select(ensembl_gene_id, entrezgene_id, hgnc_symbol, everything())
degs_path2 <- file.path(degs_dir, "PKD1_ShrinkAll.csv")
write_csv(shrink_degs1, degs_path2)

shrink_degs2 <- filter(shrink_degs1, !is.na(padj), padj < 0.05)
degs_path3 <- file.path(degs_dir, "PKD1_ShrinkSig.csv")
write_csv(shrink_degs2, degs_path3)
shrink_degs3 <- filter(shrink_degs2, abs(log2FoldChange) >= 1)
degs_path4 <- file.path(degs_dir, "PKD1_ShrinkPass.csv")
write_csv(shrink_degs3, degs_path4)