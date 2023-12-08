### Seurat V3
library(BoutrosLab.plotting.general)
source('tune_tsne.R');
source('find_clusters_snn.R');
library(dplyr)
library("scater")
library("scran")
library(monocle)
packageVersion("monocle")
library(canprot)
source("../functions.seurat.R")
library(Seurat)
packageVersion("Seurat")
# [1] '3.2.0'
# BiocManager::available()


# saveRDS(tumor, "daqingge.2.tumor.SCTpure.cycleRgrsOut.rds")
# tumor=readRDS("daqingge.2.tumor.SCTpure.cycleRgrsOut.rds")


### read 10X results
WT <- Read10X(data.dir = "./WT/outs/raw_feature_bc_matrix")
KO <- Read10X(data.dir = "./KO/outs/raw_feature_bc_matrix")
inte.list=list()# the list prepared for further dataset integration
inte.list[["WT"]] <- CreateSeuratObject(counts = WT, project = "WT", min.cells = 3, min.features = 100)
inte.list[["KO"]] <- CreateSeuratObject(counts = KO, project = "KO", min.cells = 3, min.features = 100)
for (i in 1:length(inte.list))
{
  inte.list[[i]][["percent.mt"]] <- PercentageFeatureSet(inte.list[[i]], pattern = "^MT-")
}

######################## check quality --->>>
hist(inte.list[[1]]@meta.data$"percent.mt")
dev.off()
hist(inte.list[[2]]@meta.data$"percent.mt")
dev.off()
quantile(inte.list[[1]]@meta.data$"percent.mt", 0.9)
# 9.7005
quantile(inte.list[[2]]@meta.data$"percent.mt", 0.9)
# 9.411921
summary(inte.list[[1]]@meta.data$"percent.mt")
summary(inte.list[[2]]@meta.data$"percent.mt")
i=1
i=2
pdf(paste("quality_", i,".pdf", sep=""))
VlnPlot(inte.list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,  pt.size = 0.01)
dev.off()
i=1
i=2
pdf(paste("quality_corr_", i,".pdf", sep=""))
plot1 <- FeatureScatter(inte.list[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(inte.list[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

### choose good-quality cells
str(colnames(inte.list[[1]]))
# chr [1:17558] "AAACCCAAGAACGTGC-1" "AAACCCAAGATTGCGG-1" ...
str(colnames(inte.list[[2]]))
# chr [1:12895] "AAACCCAAGCACTAAA-1" "AAACCCACAAGTATCC-1" ...
for (i in 1:length(inte.list))
{
  inte.list[[i]]<- subset(inte.list[[i]], subset = nFeature_RNA > 2500  & percent.mt < 20)
}
str(colnames(inte.list[[1]]))
# chr [1:4628] "AAACCCACAAGCCCAC-1" "AAACCCACAGTTGTTG-1" ...
str(colnames(inte.list[[2]]))
# chr [1:2946] "AAACCCAAGCACTAAA-1" "AAACCCACAAGTATCC-1" ...

summary(inte.list[[1]]@meta.data$nFeature_RNA
        # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
        # 2502    4796    5810    5729    6773    9826
        summary(inte.list[[2]]@meta.data$nFeature_RNA)
        # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
        # 2501    4348    5452    5473    6572   10318
        table(inte.list[[1]]@meta.data$nFeature_RNA>3000)
        # WT
        # FALSE  TRUE
        #   195  4433
        table(inte.list[[2]]@meta.data$nFeature_RNA>3000)
        # KO
        # FALSE  TRUE
        #   210  2736
        ######################## <<<--- check quality 
        integrated <- merge(inte.list[["WT"]], y = inte.list[["KO"]], add.cell.ids = c("WT", "KO"), project = "daqingge.2")
        str(integrated)
        table(integrated@meta.data$orig.ident)
        #   KO   WT
        # 2946 4628
        s.genes <- cc.genes$s.genes
        g2m.genes <- cc.genes$g2m.genes
        integrated<- CellCycleScoring(integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
        
        num.tab=table(integrated@meta.data$Phase, integrated@meta.data$phenotype)
        num.tab
        num.sum=apply(num.tab,2,sum)
        for(i in 1:nrow(num.tab))
        {
          num.tab[i, ]=num.tab[i, ]/num.sum
        }
        round(num.tab,2)
        #       KO   WT
        # G1  0.69 0.49
        # G2M 0.14 0.20
        # S   0.17 0.31
        
        tumor<- SCTransform(integrated, vars.to.regress = c("percent.mt", "Phase"), verbose = FALSE, return.only.var.genes=FALSE)
        DefaultAssay(tumor) <- "SCT"
        str(tumor@assays$SCT@var.features)
        # chr [1:3000] "SPANXB1" "MALAT1" "
        tumor <- RunPCA(tumor, assay="SCT", verbose = FALSE)
        ElbowPlot(tumor,  ndims = 50)
        dev.off()
        tumor <- RunUMAP(tumor, dims = 1:50, verbose = FALSE)
        tumor <- FindNeighbors(tumor, dims = 1:50, verbose = FALSE)
        tumor <- FindClusters(tumor, verbose = FALSE)
        DefaultAssay(tumor)="RNA"
        tumor <- NormalizeData(tumor, normalization.method = "LogNormalize", scale.factor = 10000)
        tumor <- ScaleData(tumor, features = rownames(tumor), assay="RNA" )
        DefaultAssay(tumor)="SCT"
        
        pdf("dimplot.pdf")
        DimPlot(tumor, label = TRUE) + NoLegend()
        dev.off()
        tumor@meta.data$fig.cluster=Idents(tumor)
        tumor@meta.data$phenotype=tumor@meta.data$orig.ident
        tumor@meta.data$phenotype=factor(tumor@meta.data$phenotype, levels=c("WT", "KO"))
        table(tumor@meta.data$fig.cluster)
        
        # Fig.A
        pdf("dimplot.cluster.pdf",width=5, height=4.5)
        DimPlot(tumor, label = TRUE, group.by="fig.cluster")
        dev.off()
        
        pdf("dimplot.phenotype.pdf",width=5, height=5)
        DimPlot(tumor, label = TRUE, group.by="phenotype")
        dev.off()
        
        pdf("dimplot.phase.pdf",width=5, height=5)
        DimPlot(tumor, label = TRUE, group.by="Phase")
        dev.off()
        
        # eFig.A
        pdf("dimplot.split.by.phenotype.pdf",width=9, height=5)
        DimPlot(tumor, label = TRUE, split.by = 'phenotype', group.by="fig.cluster") + NoLegend()
        dev.off()
        
        # eFig.C
        pdf("vlnplot.nGene.pdf",width=10, height=2)
        tocheck=c( "nFeature_RNA" )
        temp=VlnPlot(tumor, features = tocheck, group.by="fig.cluster",assay = "RNA", pt.size=0)
        for(i in 1:length(tocheck)) 
        {
          temp[[i]] <- temp[[i]]+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)+geom_hline(yintercept=median(tumor@meta.data[, tocheck[i]]), linetype="dashed")
        }
        print(temp)
        dev.off()
        
        
        ### identify DEGs as markers of clusters
        tumor <- SetIdent(tumor, value = "phenocluster")
        # wt.tumor=subset(x=tumor, subset =(phenotype=="WT"))
        tumor.markers <- FindAllMarkers(tumor, assay="SCT")
        tumor.markers=tumor.markers[order(tumor.markers$cluster, tumor.markers$avg_logFC, decreasing=c("FALSE", "TRUE")), ]
        # saveRDS(tumor.markers, "daqingge.2.tumor.markers.SCTcycleRMed.rds")
        # tumor.markers=readRDS("daqingge.2.tumor.markers.SCTcycleRMed.rds")
        
        
        ### cell state transition trajectory analysis
        ### use Monocle2
        # refer to: http://cole-trapnell-lab.github.io/monocle-release/docs/#installing-monocle
        # BiocManager::install("monocle")
        require(monocle)
        # tumor=readRDS("daqingge.2.tumor.SCTpure.cycleRgrsOut.rds")
        ### check cell lineage signatures in trajectory of different phenotypes
        {
          # norme=readRDS("yulu.normal.bladder.seurat.rds")
          # tumor=readRDS("daqingge.2.tumor.SCTpure.cycleRgrsOut.V2.rds")
          
          # get signature genes
          norme <- SetIdent(norme, value = "celltype")
          norme.celltype.markers <- FindAllMarkers(norme, assay="SCT")
          # saveRDS(norme.celltype.markers, "norme.celltype.markers.rds")
          # norme.celltype.markers =readRDS("norme.celltype.markers.rds")
          
          markers= norme.celltype.markers
          markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
          markers$foldchange=2^(markers$avg_log2FC)
          markers=markers[markers$p_val_adj<0.05 & markers$avg_log2FC>log2(1.5),]
          write.csv(markers[,c("gene","foldchange","pct.1", "pct.2", "p_val_adj", "cluster")], file=paste("norme.cluster.markers","csv",sep="."), row.name=T)
          
          celltypes=levels(markers$cluster)
          cellsig=list()
          for(i in celltypes)
          {
            cellsig[[i]]=markers[markers$cluster==i, ]$gene
          }
          
          # wt=readRDS("monocle.wt.V2.rds")
          # ko=readRDS("monocle.ko.V2.rds")
          
          # get signature
          for(i in celltypes)
          {
            temp=cellsig[[i]][cellsig[[i]]%in%rownames(tumor)]
            tumor@meta.data[, i]=colMeans(tumor@assays$SCT@scale.data[temp, ])
          }
          
        ### --->>> WT part evolution 
        temp=subset( tumor, subset=(phenotype=="WT") )
        temp<-temp[, temp$fig.cluster %in% c( "1","2","3","4","5","7","9","10")]
        gene_metadata=data.frame(gene_short_name=rownames(temp@assays$SCT@counts))
        rownames(gene_metadata)=rownames(temp@assays$SCT@counts)
        wt <- newCellDataSet(  temp@assays$SCT@counts,
                               phenoData = new("AnnotatedDataFrame", temp@meta.data),
                               featureData = new("AnnotatedDataFrame", gene_metadata), 
                               expressionFamily=negbinomial.size() )
        table(wt@phenoData@data$phenocluster)
        table(wt@phenoData@data$phenotype)
        wt <- estimateSizeFactors(wt)
        wt <- estimateDispersions(wt)
        
        wt <- detectGenes(wt, min_expr = 0.1)
        ### only keep expressed genes
        expressed_genes <- row.names(wt)[wt@featureData@data$num_cells_expressed>= 10]
        wt <- wt[expressed_genes,]
        
        ### use all significant markers of clusters as ordering genes
        tumor.markers=readRDS("daqingge.2.tumor.markers.SCTcycleRMed.rds") 
        
        markers=tumor.markers
        markers=markers[markers$p_val_adj<0.05, ]
        markers=markers[abs(markers$avg_logFC)>1.09, ]
        markers$foldChange=exp(markers$avg_logFC)
        # markers=markers[order(markers$cluster, markers$foldChange, decreasing=c("FALSE", "TRUE")), ]
        markers=markers[order(markers$foldChange, decreasing=TRUE), ]
        ordering.genes=unique(markers$gene)
        # ordering.genes=unique(rownames(markers))
        
        wt <- setOrderingFilter(wt,  ordering.genes) # 
        pdf("plot_ordering_genes.pdf")
        plot_ordering_genes(wt)
        dev.off()
        
        wt <- reduceDimension(wt, max_components = 2,
                              method = 'DDRTree')
        
        wt <- orderCells(wt)
        
        plot_cell_trajectory(wt, color_by = "phenocluster")
        dev.off()
        
        # set a indicator for time, and sort again
        wt@phenoData@data$timeIset=rep(1, ncol(wt))
        wt@phenoData@data$timeIset[wt@phenoData@data$phenocluster=="WT.3"]=0
        wt <- orderCells(wt, root_state = wt@phenoData@data$timeIset)
        
        # eFig.D.1
        pdf("plot_cell_trajectory_byPseudotime_WT.pdf", width=5, height=5)
        # png("plot_cell_trajectory_byPseudotime_WT.png")
        plot_cell_trajectory(wt, color_by = "Pseudotime", cell_size=0.4)
        dev.off()
        
        # eFig.D.2
        pdf("plot_cell_trajectory_fig.cluster_WT.pdf", width=5, height=5)
        # png("plot_cell_trajectory_byPseudotime_WT.png")
        plot_cell_trajectory(wt, color_by = "fig.cluster", cell_size=0.4)
        dev.off()
        
        wt.time=wt@phenoData@data$Pseudotime
        names(wt.time)=rownames(wt@phenoData@data)
        
        png("plot_cell_trajectory_allIn1_WT.png")
        plot_cell_trajectory(wt, color_by = "phenocluster")
        dev.off()
        
        png("plot_cell_trajectory_details_phenocluster_WT.png", width=30*100, height=6*100)
        temp=plot_cell_trajectory(wt, color_by = "phenocluster") +
          facet_wrap(~fig.cluster, nrow = 1)
        print(temp)
        dev.off()
        
        # eFig.D.3
        png("plot_cell_trajectory_details_cluster_WT.png", width=30*100, height=5*100)
        temp=plot_cell_trajectory(wt, color_by = "fig.cluster", cell_size = 0.8) +
          facet_wrap(~fig.cluster, nrow = 1)
        print(temp)
        dev.off()
        pdf("plot_cell_trajectory_details_cluster_WT.pdf", width=30, height=5)
        temp=plot_cell_trajectory(wt, color_by = "fig.cluster", cell_size = 0.4) + facet_wrap(~fig.cluster, nrow = 1)
        print(temp)
        dev.off()
        
        saveRDS(wt, "monocle.wt.rds")
        wt=readRDS("monocle.wt.rds")
        
        # find DEGs along Pseudotime 
        diff_test_res.wt <- differentialGeneTest(wt, fullModelFormulaStr = "~sm.ns(Pseudotime)")
        # saveRDS(diff_test_res.wt, "diff_test_res.wt.rds")
        diff_test_res.wt=readRDS("diff_test_res.wt.rds")
        # diff_test_res.bak=diff_test_res
        diff_test_res =  diff_test_res.wt
        str(diff_test_res[,c("gene_short_name", "pval", "qval")])
        diff_test_res=subset(diff_test_res, qval < 0.1 & pval < 0.01)
        sig_gene_names=intersect(diff_test_res[,"gene_short_name"], subtype_markers$"TF_fromLieBing")
        imp.genes=c()
        imp.genes=c(imp.genes, c("BMPR2", "CCND2", "SOX4"))
        imp.genes=c(imp.genes, c("CDH2", "TWIST2", "HHIP"))
        imp.genes=c(imp.genes, c("COL5A1", "ITGA2", "MATN2", "SEMA3A"))
        imp.genes=c(imp.genes, sig_gene_names)
        imp.gen es=intersect(rownames(wt), imp.genes)
        pdf("heatmap.trajectory.WT.pdf")
        plot_pseudotime_heatmap(wt[imp.genes,],
                                num_clusters = 3,
                                cores = 1,
                                show_rownames = T)
        dev.off()
        
        pdf("exp.along.time.WT.pdf")
        imp.genes=c()
        imp.genes=c(imp.genes, c("BMPR2", "CCND2", "SOX4"))
        imp.genes=c(imp.genes, c("CDH2",  "HHIP"))# "TWIST2",
        imp.genes=c(imp.genes, c("COL5A1", "ITGA2", "MATN2", "SEMA3A"))
        genes <- imp.genes
        plot_genes_branched_pseudotime(wt[genes,],
                                       branch_point = 1,
                                       color_by = "fig.cluster",
                                       ncol = 1)
        dev.off()
        
        ### <<<--- WT part evolution 
        
        ### --->>> KO part evolution
        temp=subset( tumor, subset=(phenotype=="KO") )
        gene_metadata=data.frame(gene_short_name=rownames(temp@assays$SCT@counts))
        rownames(gene_metadata)=rownames(temp@assays$SCT@counts)
        ko <- newCellDataSet(  temp@assays$SCT@counts,
                               phenoData = new("AnnotatedDataFrame", temp@meta.data),
                               featureData = new("AnnotatedDataFrame", gene_metadata), 
                               expressionFamily=negbinomial.size() )
        
        table(ko@phenoData@data$phenocluster)
        table(ko@phenoData@data$phenotype)
        ko <- estimateSizeFactors(ko)
        ko <- estimateDispersions(ko)
        
        ko <- detectGenes(ko, min_expr = 0.1)
        ### only keep expressed genes
        expressed_genes <- row.names(ko)[ko@featureData@data$num_cells_expressed>= 10]
        ko <- ko[expressed_genes,]
        
        ### use all significant markers of clusters as ordering genes
        tumor.markers=readRDS("daqingge.2.tumor.markers.SCTcycleRMed.rds") 
        
        markers=tumor.markers
        markers=markers[markers$p_val_adj<0.05, ]
        markers=markers[abs(markers$avg_logFC)>1.09, ]
        markers$foldChange=exp(markers$avg_logFC)
        # markers=markers[order(markers$cluster, markers$foldChange, decreasing=c("FALSE", "TRUE")), ]
        markers=markers[order(markers$foldChange, decreasing=TRUE), ]
        ordering.genes=unique(markers$gene)
        # ordering.genes=unique(rownames(markers))
        
        ko <- setOrderingFilter(ko,  ordering.genes) # 
        pdf("plot_ordering_genes.pdf")
        plot_ordering_genes(ko)
        dev.off()
        
        ko <- reduceDimension(ko, max_components = 2,
                              method = 'DDRTree')
        
        ko <- orderCells(ko)
        
        # set a indicator for time
        ko@phenoData@data$timeIset=rep(1, ncol(ko))
        ko@phenoData@data$timeIset[ko@phenoData@data$phenocluster=="KO.3"]=0
        ko <- orderCells(ko, root_state = ko@phenoData@data$timeIset)
        
        # eFig.D.1
        pdf("plot_cell_trajectory_byPseudotime_KO.pdf", width=5, height=5)
        plot_cell_trajectory(ko, color_by = "Pseudotime", cell_size=0.4)
        dev.off()
        
        # eFig.D.2
        pdf("plot_cell_trajectory_fig.cluster_KO.pdf", width=5, height=5)
        plot_cell_trajectory(ko, color_by = "fig.cluster", cell_size=0.4)
        dev.off()
        
        ko.time=ko@phenoData@data$Pseudotime
        names(ko.time)=rownames(ko@phenoData@data)
        
        png("plot_cell_trajectory_allIn1_KO.png")
        plot_cell_trajectory(ko, color_by = "phenocluster")
        dev.off()
        
        # pdf("plot_cell_trajectory_details.pdf", width=30, height=6)
        png("plot_cell_trajectory_details_phenocluster_KO.png", width=30*100, height=6*100)
        temp=plot_cell_trajectory(ko, color_by = "phenocluster") +
          facet_wrap(~fig.cluster, nrow = 1)
        print(temp)
        dev.off()
        
        pdf("plot_cell_trajectory_details_KO.pdf", width=30, height=5)
        temp=plot_cell_trajectory(ko, color_by = "fig.cluster", cell_size = 0.4) +
          facet_wrap(~fig.cluster, nrow = 1)
        print(temp)
        dev.off()
        saveRDS(ko, "monocle.ko.rds")
        ko=readRDS("monocle.ko.rds")
        ### plot smoothed expression along pseudotime
        ### refer to https://stackoverflow.com/questions/11014804/plotting-multiple-smooth-lines-from-a-dataframe
        # tumor@meta.data$time
        
        temp=subset(tumor, phenotype=="KO")
        # temp=subset(tumor, phenotype=="WT")
        dataset <- data.frame( xval = temp@meta.data$time, Umbrella = temp@meta.data$umbrella, Intermediate = temp@meta.data$intermediate, Basal = temp@meta.data$basal )
        #convert data to long format
        library(reshape)
        Molten <- melt(dataset, id.vars = "xval")
        #plot it
        library(ggplot2)
        # ggplot(Molten, aes(x = xval, y = value, colour = variable)) + geom_smooth() + geom_point()
        #some tweaking
        temp=ggplot(Molten , aes(x = xval, y = value, colour = variable)) + 
          geom_smooth(se = FALSE) + geom_point(size =0.1) + theme_bw() + 
          scale_x_continuous("Pseudotime") + scale_y_continuous("Signature of cell lineage") +
          scale_colour_discrete("")
        print(temp)
        dev.off()
        }
        
        }

        ### --->>> identify phenotype-time DEGs
        ### find genes with expression correlated with time and different across phenotypes
        temp=tumor
        temp@meta.data$time=rep(-1, ncol(temp))
        temp@meta.data[names(wt.time), ]$time=wt.time
        temp@meta.data[names(ko.time), ]$time=ko.time
        temp@meta.data$phenotype1=rep(-1, ncol(temp))
        temp@meta.data[names(wt.time), ]$phenotype1="Inva"
        temp@meta.data[names(ko.time), ]$phenotype1="Non-Inv"
        temp=tumor
        # KO's genes
        ko.seurat=subset( tumor, subset=(phenotype1=="Non-Inv") )
        # change @data to @scale.data 20210524
        ko.cor=apply(ko.seurat@assays$SCT@scale.data, 1, function(x) cor(x, ko.seurat@meta.data$time, method = c("spearman")))
        names(ko.cor)=rownames(ko.seurat@assays$SCT@scale.data)
        ko.cor=ko.cor[!is.na(ko.cor)]
        ko.cor=sort(ko.cor, decreasing=TRUE)
        # saveRDS(ko.cor, "ko.cor.rds")
        # ko.cor=readRDS("ko.cor.rds")
        
        ko.p=apply(ko.seurat@assays$SCT@scale.data, 1, function(x) cor.test(x, ko.seurat@meta.data$time, method = c("spearman"))$p.value)
        names(ko.p)=rownames(ko.seurat@assays$SCT@scale.data)
        ko.p=ko.p[!is.na(ko.p)]
        ko.p.adj=p.adjust(ko.p)
        table(ko.p<0.01)
        table(ko.p.adj<0.01)
        # saveRDS(ko.p, "ko.p.rds")
        # ko.p=readRDS("ko.p.rds")
        
        # change wt.cor>0.25 to >0.2
        ko.cor.genes=names(ko.cor[ko.cor>0.2])
        ko.p.genes=names(ko.p.adj[ko.p.adj<0.01])
        ko.cor.genes=intersect(ko.cor.genes, ko.p.genes)
        str(ko.cor.genes)
        # chr [1:3132] "B2M" "BST2" "IFI6" "HLA-C" "ITM2B" "HLA-B" "HIST1H2AC" ...
        
        # WT's genes
        wt.seurat=subset( tumor, subset=(phenotype1=="Inva") )
        wt.cor=apply(wt.seurat@assays$SCT@data, 1, function(x) cor(x, wt.seurat@meta.data$time, method = c("spearman")))
        names(wt.cor)=rownames(wt.seurat@assays$SCT@data)
        wt.cor=wt.cor[!is.na(wt.cor)]
        wt.cor=sort(wt.cor, decreasing=TRUE)
        # saveRDS(wt.cor, "wt.cor.rds")
        # wt.cor=readRDS("wt.cor.rds")
        
        wt.p=apply(wt.seurat@assays$SCT@data, 1, function(x) cor.test(x, wt.seurat@meta.data$time, method = c("spearman"))$p.value)
        names(wt.p)=rownames(wt.seurat@assays$SCT@data)
        wt.p=wt.p[!is.na(wt.p)]
        wt.p.adj=p.adjust(wt.p)
        table(wt.p<0.01)
        table(wt.p.adj<0.01)
        # saveRDS(wt.p, "wt.p.rds")
        # wt.p=readRDS("wt.p.rds")
        
        # change wt.cor>0.25 to >0.2
        wt.cor.genes=names(wt.cor[wt.cor>0.2])
        wt.p.genes=names(wt.p.adj[wt.p.adj<0.01])
        wt.cor.genes=intersect(wt.cor.genes, wt.p.genes)
        str(wt.cor.genes)
        # chr [1:500] "TSC22D1" "PSAP" "IFITM3" "GRN" "PFDN5" "ALKBH7" "RPL34" ...
        
        # eFig.D.3.alternate
        pdf("WT.vlnplot.time.pdf", width=10, height=5)
        VlnPlot(wt.seurat, features = c("Pseudotime"),group.by="fig.cluster",assay = "RNA", pt.size=0)+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
        dev.off()
        
        # eFig.D.6.alternate
        pdf("KO.vlnplot.time.pdf", width=10, height=5)
        VlnPlot(ko.seurat, features = c("Pseudotime"),group.by="fig.cluster",assay = "RNA", pt.size=0)+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
        dev.off()
        
        temp=c("APBB2", "COL5A1", "HSPG2", "MATN2", "SEMA3A", "SEMA3C", "BMPR2", "CCND2", "FOXN3", "FSTL1", "BST2", "PODXL", "CAV1")
        ko.p.adj[temp]
        ko.cor[temp]
        
        temp=c("EDIL3", "HES1", "NNMT")
        wt.p.adj[temp]
        wt.cor[temp]
        
        temp=c("CDH2", "TWIST1")
        ko.p.adj[temp]
        ko.cor[temp]
        wt.p.adj[temp]
        wt.cor[temp]
        # check cell lineage markers in clusters. 
        # markers are obtained from "Single-Cell Transcriptomic Map of the Human and Mouse Bladders"
        # basal, intermediate, TNNT1+, umbrella
        tocheck=c("KRT5", "KRT17",     "KRT13",     "TNNT1",     "UPK1A", "UPK2")
        # tocheck=c("ITGA2")
        pdf("celllineage.marker.vlnplot.across.pheno.pdf", width=10, height=5)
        VlnPlot(tumor, features = tocheck, group.by="fig.cluster",assay = "SCT", pt.size=0)+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)
        dev.off()
        
        tocheck=c("KRT5", "KRT17",     "KRT13",     "TNNT1",     "UPK1A", "UPK2")
        # tumor@meta.data$phenotype=factor(tumor@meta.data$phenotype, levels=c("WT", "KO"))
        # tumor <- ScaleData(tumor, features = rownames(tumor), assay="RNA" )
        pdf("heatmap.DEG.phenotype.pdf")
        DoHeatmap( tumor, features = tocheck, group.by = "fig.cluster", slot="scale.data", assay="RNA", angle = 0 ) 
        dev.off()
        
        
        tumor.cluster.markers <- FindAllMarkers(tumor, assay="SCT")
        saveRDS(tumor.cluster.markers, "tumor.cluster.markers.rds")
        # tumor.cluster.markers =readRDS("tumor.cluster.markers.rds")
        
        markers=tumor.cluster.markers
        markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
        markers$foldchange=2^(markers$avg_log2FC)
        markers=markers[markers$p_val_adj<0.05 & abs(markers$avg_log2FC)>log2(1.5) & ( markers$pct.1>0.1 | markers$pct.2>0.1 ),]
        write.csv(markers[,c("gene","foldchange","pct.1", "pct.2", "p_val_adj", "cluster")], file=paste("tumor.cluster.markers","csv",sep="."), row.name=T)
       ##3
        ### plot some signatures
        # read signatures from daqingge
        # differentiation related gene sets
        sigset=list()
        dif.sig=read.csv("./KamounPathwayGenes.20201204.csv", stringsAsFactors=FALSE)
        sig.names=unique(dif.sig$ont)
        for(i in sig.names)
        {
          sigset[[i]]=dif.sig$gene[dif.sig$ont==i]
        }
        # epithelial signature (etan) and mesenchymal (mtan) assocaited with TAN dataset
        sigset[["E.TAN"]]=c("KRT19", "RAB25", "EPCAM", "CDH1")
        sigset[["M.TAN"]]=c("VIM", "ZEB1", "ZEB2", "CDH2", "TWIST1", "SNAI2")
        # epithelial signature (etan) and mesenchymal (mtan) assocaited with TCGA dataset
        sigset[["E.TCGA"]]=c("OCLN", "DSP", "CDH1")
        sigset[["M.TCGA"]]=c("VIM", "MMP2", "FN1", "SNAI1", "MMP9", "FOXC2", "CDH2", "GSC", "TWIST1", "SNAI2", "MMP3")
        
        temp.matrix=matrix(NA, 0, ncol(tumor@assays$RNA@scale.data))
        for(i in names(sigset))
        {
          sigset[[i]]=intersect(sigset[[i]], rownames(tumor@assays$RNA@scale.data))
          temp.row=colMeans(tumor@assays$RNA@scale.data[sigset[[i]], ])
          temp.matrix=rbind(temp.matrix, temp.row)
        }
        rownames(temp.matrix)=names(sigset)
        temp.matrix=rbind(temp.matrix, temp.matrix["M.TAN",]-temp.matrix["E.TAN",])
        temp.matrix=rbind(temp.matrix, temp.matrix["M.TCGA",]-temp.matrix["E.TCGA",])
        rownames(temp.matrix)[(nrow(temp.matrix)-1):nrow(temp.matrix)]=c("EMT.TAN", "EMT.TCGA")
        temp.matrix[1:5, 1:5]
        ### DO NOT save "tumor" for further use, after running following 2 rows 
        ### tumor@assays$RNA@scale.data=rbind(tumor@assays$RNA@scale.data, temp.matrix)
        ### tumor@assays$RNA@data=tumor@assays$RNA@scale.data
        
        # Fig.E
        tocheck=c(names(sigset), c("EMT.TAN", "EMT.TCGA"))
        tocheck=c("EMT.TAN", "EMT.TCGA")
        for(i in tocheck)
        {
          pdf(paste("vlnplot.", i, ".pdf", sep=""), width=5, height=3)
          temp=VlnPlot(tumor, features = i, group.by="pheno.time.group",assay = "RNA", pt.size=0, cols=c("#99CCFF", "#9999FF", "#000099", "#FFCCCC", "#FF6666", "#660000"))
          temp=temp+stat_summary(fun.y = mean, geom='point', size = 8, colour = "black", shape = 95)
          print(temp)
          dev.off()
        }
        
        # Fig.F
        tocheck=c(names(sigset), c("EMT.TAN", "EMT.TCGA"))
        tocheck=c("EMT.TAN", "EMT.TCGA")
        for(i in tocheck)
        {
          pdf(paste("featureplot.", i, ".pdf", sep=""), width=5, height=5)
          temp=FeaturePlot(tumor, features = i, cols=c("lightblue", "red"))
          print(temp)
          dev.off()
        }
        
        ### normal bladder monocle trajectory
        {
          require(monocle)
          # norme=readRDS("yulu.normal.bladder.seurat.rds")
          # tumor.cluster.markers =readRDS("tumor.cluster.markers.rds")
          
          ### --->>> normem part evolution 
          # temp=subset( norme, subset=(phenotype=="normem") )
          temp=norme
          
          gene_metadata=data.frame(gene_short_name=rownames(temp@assays$SCT@counts))
          rownames(gene_metadata)=rownames(temp@assays$SCT@counts)
          normem <- newCellDataSet(  temp@assays$SCT@counts,
                                     phenoData = new("AnnotatedDataFrame", temp@meta.data),
                                     featureData = new("AnnotatedDataFrame", gene_metadata), 
                                     expressionFamily=negbinomial.size() )
          
          table(normem@phenoData@data$phenocluster)
          table(normem@phenoData@data$phenotype)
          normem <- estimateSizeFactors(normem)
          normem <- estimateDispersions(normem)
          
          normem <- detectGenes(normem, min_expr = 0.1)
          ### only keep expressed genes
          expressed_genes <- row.names(normem)[normem@featureData@data$num_cells_expressed>= 10]
          normem <- normem[expressed_genes,]
          
          ### use all significant markers of clusters as ordering genes
          # norme <- SetIdent(norme, value = "celltype")
          # norme.celltype.markers <- FindAllMarkers(norme, assay="SCT")
          # saveRDS(norme.celltype.markers, "norme.celltype.markers.rds")
          norme.celltype.markers = readRDS("norme.celltype.markers.rds")
          
          markers=norme.celltype.markers
          markers=markers[markers$p_val_adj<0.05, ]
          markers=markers[abs(markers$avg_log2FC)>log2(3), ]
          markers$foldChange=2^(markers$avg_log2FC)
          # markers=markers[order(markers$cluster, markers$foldChange, decreasing=c("FALSE", "TRUE")), ]
          markers=markers[order(markers$foldChange, decreasing=TRUE), ]
          ordering.genes=unique(markers$gene)
          # ordering.genes=unique(rownames(markers))
          
          normem <- setOrderingFilter(normem,  ordering.genes) # 
          pdf("plot_ordering_genes.pdf")
          plot_ordering_genes(normem)
          dev.off()
          
          normem <- reduceDimension(normem, max_components = 2,
                                    method = 'DDRTree')
          
          normem <- orderCells(normem)
          
          plot_cell_trajectory(normem, color_by = "celltype")
          dev.off()
          
          
          plot_cell_trajectory(normem, color_by = "State")
          dev.off()
          
          
          # check cell types distribution and monocle's states in norm epithelial cells
          table(normem.state=normem@phenoData@data$State)
          # normem.state
          #    1    2    3    4    5    6    7
          # 3495  441   19  922  585 1420  647
          normem.state=normem@phenoData@data$State
          names(normem.state)=rownames(normem@phenoData@data)
          norme@meta.data[names(normem.state), "state"]=normem.state
          
          pdf("dimplot.celltype.pdf",width=5, height=4.5)
          DimPlot(norme, label = TRUE, group.by="celltype")
          dev.off()
          
          pdf("dimplot.state.pdf",width=5, height=4.5)
          DimPlot(norme, label = TRUE, group.by="state")
          dev.off()
          
          tocheck=c("nFeature_RNA")
          for(i in tocheck)
          {
            pdf(paste("vlnplot.", i, ".pdf", sep=""), width=5, height=3)
            temp=VlnPlot(norme, features = i, group.by="state",assay = "RNA", pt.size=0)
            temp=temp+stat_summary(fun.y = mean, geom='point', size = 8, colour = "black", shape = 95)
            print(temp)
            dev.off()
          }
          
          # Fig.F
          tocheck=c("nFeature_RNA")
          for(i in tocheck)
          {
            pdf(paste("featureplot.", i, ".pdf", sep=""), width=5, height=5)
            temp=FeaturePlot(norme, features = i, cols=c("lightblue", "red"))
            print(temp)
            dev.off()
          }
          
          # set a indicator for time, and sort again
          table(normem@phenoData@data$State)
          #    1    2    3    4    5    6    7
          # 3495  441   19  922  585 1420  647
          str(normem@phenoData@data$State)
          # Factor w/ 7 levels "1","2","3","4",..: 7 7 7 6 7 7 7 6 7 7 ...
          # as state 6 has the highest expressed gene number
          
          normem <- orderCells(normem, root_state = "1")
          
          # eFig.D.1
          pdf("plot_cell_trajectory_byPseudotime_normem.pdf", width=5, height=5)
          # png("plot_cell_trajectory_byPseudotime_normem.png")
          plot_cell_trajectory(normem, color_by = "Pseudotime", cell_size=0.4)
          dev.off()
          
          # eFig.D.2
          pdf("plot_cell_trajectory_fig.cluster_normem.pdf", width=5, height=5)
          # png("plot_cell_trajectory_byPseudotime_normem.png")
          plot_cell_trajectory(normem, color_by = "fig.cluster", cell_size=0.4)
          dev.off()
          
          # eFig.D.2
          pdf("plot_cell_trajectory_celltype_normem.pdf", width=5, height=5)
          # png("plot_cell_trajectory_byPseudotime_normem.png")
          plot_cell_trajectory(normem, color_by = "celltype", cell_size=0.4)
          dev.off()
          
          normem.time=normem@phenoData@data$Pseudotime
          names(normem.time)=rownames(normem@phenoData@data)
          
          norme@meta.data$time=rep(-1, ncol(norme))
          norme@meta.data[names(normem.time), ]$time=normem.time
          norme@meta.data$time=round(norme@meta.data$time, 1 )
          
          png("plot_cell_trajectory_allIn1_normem.png")
          plot_cell_trajectory(normem, color_by = "phenocluster")
          dev.off()
          
          png("plot_cell_trajectory_details_phenocluster_normem.png", width=30*100, height=6*100)
          temp=plot_cell_trajectory(normem, color_by = "phenocluster") +
            facet_wrap(~fig.cluster, nrow = 1)
          print(temp)
          dev.off()
          
          # eFig.D.3
          png("plot_cell_trajectory_details_cluster_normem.png", width=25*100, height=5*100)
          temp=plot_cell_trajectory(normem, color_by = "fig.cluster", cell_size = 0.8) +
            facet_wrap(~fig.cluster, nrow = 1)
          print(temp)
          dev.off()
          pdf("plot_cell_trajectory_details_cluster_normem.pdf", width=25, height=5)
          temp=plot_cell_trajectory(normem, color_by = "fig.cluster", cell_size = 0.4) + facet_wrap(~fig.cluster, nrow = 1)
          print(temp)
          dev.off()
          
          saveRDS(normem, "monocle.normem.V2.rds")
          # normem=readRDS("monocle.normem.V2.rds")
          
          ### plot celltype signature along pseudotime
          {
            norme.celltype.markers =readRDS("norme.celltype.markers.rds")
            
            markers= norme.celltype.markers
            markers=markers[order(markers$cluster, markers$avg_log2FC,decreasing = c(FALSE,TRUE)), ]
            markers$foldchange=2^(markers$avg_log2FC)
            markers=markers[markers$p_val_adj<0.05 & markers$avg_log2FC>log2(1.5),]
            write.csv(markers[,c("gene","foldchange","pct.1", "pct.2", "p_val_adj", "cluster")], file=paste("norme.cluster.markers","csv",sep="."), row.name=T)
            
            celltypes=levels(markers$cluster)
            cellsig=list()
            for(i in celltypes)
            {
              cellsig[[i]]=markers[markers$cluster==i, ]$gene
            }
            
            # wt=readRDS("monocle.wt.V2.rds")
            # ko=readRDS("monocle.ko.V2.rds")
            
            # get signature
            for(i in celltypes)
            {
              temp=cellsig[[i]][cellsig[[i]]%in%rownames(norme@assays$SCT@scale.data)]
              norme@meta.data[, i]=colMeans(norme@assays$SCT@scale.data[temp, ])
            }
            
            ## plot smoothed expression along pseudotime
            # norme@meta.data$time
            temp=norme
            dataset <- data.frame( xval = temp@meta.data$time, Umbrella = temp@meta.data$umbrella, Intermediate = temp@meta.data$intermediate, Basal = temp@meta.data$basal )
            #convert data to long format
            library(reshape)
            Molten <- melt(dataset, id.vars = "xval")
            #plot it
            library(ggplot2)
            # ggplot(Molten, aes(x = xval, y = value, colour = variable)) + geom_smooth() + geom_point()
            #some tweaking
            temp=ggplot(Molten , aes(x = xval, y = value, colour = variable)) + 
              geom_smooth(se = FALSE) + geom_point(size =0.1) + theme_bw() + 
              scale_x_continuous("Pseudotime") + scale_y_continuous("Signature of cell lineage") +
              scale_colour_discrete("")
            print(temp)
            dev.off()
          }
          
        }      
       