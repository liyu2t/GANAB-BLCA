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

# saveRDS(tumor, "daqingge.2.tumor.SCTpure.cycleRgrsOut.rds")
# tumor=readRDS("daqingge.2.tumor.SCTpure.cycleRgrsOut.rds")


### identify DEGs as markers of clusters
tumor <- SetIdent(tumor, value = "phenocluster")
# wt.tumor=subset(x=tumor, subset =(phenotype=="WT"))
tumor.markers <- FindAllMarkers(tumor, assay="SCT")
tumor.markers=tumor.markers[order(tumor.markers$cluster, tumor.markers$avg_logFC, decreasing=c("FALSE", "TRUE")), ]
# saveRDS(tumor.markers, "daqingge.2.tumor.markers.SCTcycleRMed.rds")
# tumor.markers=readRDS("daqingge.2.tumor.markers.SCTcycleRMed.rds")

