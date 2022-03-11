### Seurat V3
library(BoutrosLab.plotting.general)
library(dplyr)
library(monocle)
packageVersion("monocle")
library(canprot)
source("../functions.seurat.R")
library(Seurat)
packageVersion("Seurat")
# [1] '3.2.2'
# BiocManager::available()
setwd("")

# saveRDS(tumor, "daqingge.2.tumor.SCTpure.cycleRgrsOut.V2.rds")
# tumor=readRDS("daqingge.2.tumor.SCTpure.cycleRgrsOut.V2.rds")

### read RDS produced by "scRNA.step1.R"
# tumor=readRDS("daqingge.2.tumor.SCTpure.cycleRgrsOut.rds")

pdf("dimplot.cluster.pdf",width=5, height=4.5)
DimPlot(tumor, label = TRUE, group.by="fig.cluster")
dev.off()

pdf("dimplot.phenotype.pdf",width=5, height=5)
DimPlot(tumor, label = TRUE, group.by="phenotype")
dev.off()

dim(tumor)
# [1] 20130  7574
iden=!(tumor@meta.data$fig.cluster%in%c("0", "3", "6", "8") & tumor@meta.data$phenotype=="WT")
tumor=subset(tumor, cells=colnames(tumor)[iden])
dim(tumor)
# [1] 20130  7286
iden=!(tumor@meta.data$fig.cluster%in%c("1", "2", "4", "5", "7", "9") & tumor@meta.data$phenotype=="KO")
tumor=subset(tumor, cells=colnames(tumor)[iden])
dim(tumor)
# [1] 20130  7162

### cell state transition trajectory analysis
### use Monocle2
# refer to: http://cole-trapnell-lab.github.io/monocle-release/docs/#installing-monocle
# BiocManager::install("monocle")
require(monocle)
# tumor=readRDS("daqingge.2.tumor.SCTpure.cycleRgrsOut.V2.rds")

### --->>> WT part evolution 
  temp=subset( tumor, subset=(phenotype=="WT") )
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
  markers=markers[abs(markers$avg_logFC)>log(3), ]
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
  wt@phenoData@data$timeIset[wt@phenoData@data$phenocluster=="WT.7"]=0
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
  png("plot_cell_trajectory_details_cluster_WT.png", width=25*100, height=5*100)
  temp=plot_cell_trajectory(wt, color_by = "fig.cluster", cell_size = 0.8) +
      facet_wrap(~fig.cluster, nrow = 1)
  print(temp)
  dev.off()
  pdf("plot_cell_trajectory_details_cluster_WT.pdf", width=25, height=5)
  temp=plot_cell_trajectory(wt, color_by = "fig.cluster", cell_size = 0.4) + facet_wrap(~fig.cluster, nrow = 1)
  print(temp)
  dev.off()

  saveRDS(wt, "monocle.wt.V2.rds")
  # wt=readRDS("monocle.wt.V2.rds")

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
  png("plot_cell_trajectory_details_phenocluster_KO.png", width=20*100, height=5*100)
  temp=plot_cell_trajectory(ko, color_by = "phenocluster") +
      facet_wrap(~fig.cluster, nrow = 1)
  print(temp)
  dev.off()

  pdf("plot_cell_trajectory_details_KO.pdf", width=20, height=5)
  temp=plot_cell_trajectory(ko, color_by = "fig.cluster", cell_size = 0.4) +
      facet_wrap(~fig.cluster, nrow = 1)
  print(temp)
  dev.off()

  saveRDS(ko, "monocle.ko.V2.rds")
  # ko=readRDS("monocle.ko.V2.rds")

### <<<--- KO evolution


### give predicted pseudotime to seurat obj
tumor@meta.data$time=rep(-1, ncol(tumor))
tumor@meta.data[names(wt.time), ]$time=wt.time
tumor@meta.data[names(ko.time), ]$time=ko.time
tumor@meta.data$time=round(tumor@meta.data$time, 1 )

# pdf("dimplot.time.pdf",width=5, height=5)

# Fig.B
pdf("fplot.time.pdf", width=5, height=5)
FeaturePlot(tumor, features = c("time"), cols=c("#FFFFCC", "#CC3300"))
# LightGreen 2 DarkGreen
dev.off()
# cols=c("grey", "blue")

pdf("fplot.time.split.by.phenotype.pdf", width=10, height=5)
FeaturePlot(tumor, features = c("Pseudotime"), split.by="phenotype", cols=c("lightblue", "blue"))
dev.off()

# eFig.E
pdf("ALL.vlnplot.time.pdf", width=10, height=2)
VlnPlot(tumor, features = c("Pseudotime"),group.by="fig.cluster",assay = "RNA", pt.size=0)+stat_summary(fun.y = median, geom='point', size = 8, colour = "black", shape = 95)+geom_hline(yintercept=median(tumor@meta.data$time), linetype="dashed")
dev.off()

# plot gene expression along time
ptime=tumor@meta.data$time
quantile.time=c(0, quantile(ptime, 1/3), quantile(ptime, 2/3))
quantile.group=c("Early", "Intermediate", "Late")
tumor@meta.data$time.group=rep("NA", ncol(tumor))
for(i in 1:length(quantile.time))
{
  tumor@meta.data$time.group[ptime>=quantile.time[i]]=quantile.group[i]
}
table(tumor@meta.data$time.group)
tumor@meta.data$pheno.time.group=paste(tumor@meta.data$time.group, " ", tumor@meta.data$phenotype, sep="")
table(tumor@meta.data$pheno.time.group)
tumor@meta.data$pheno.time.group=factor(tumor@meta.data$pheno.time.group, levels=c("Early WT", "Intermediate WT", "Late WT", "Early KO", "Intermediate KO", "Late KO"))
table(tumor@meta.data$pheno.time.group)

# eFig.F
pdf("dimplot.pheno.time.group.pdf",width=5.8, height=4.5)
DimPlot(tumor, label = FALSE, group.by="pheno.time.group", cols=c("#99CCFF", "#9999FF", "#000099", "#FFCCCC", "#FF6666", "#660000"))
# lightBlue 2 darkBlue, lightRed 2 darkRed
dev.off()

pdf("dimplot.time.group.pdf",width=5.5, height=4.5)
DimPlot(tumor, label = FALSE, group.by="time.group", cols=c("grey", "blue", "darkblue"))
dev.off()

### --->>> identify phenotype-time DEGs
### find genes with expression correlated with time and different across phenotypes

  # KO's genes
  ko.seurat=subset( tumor, subset=(phenotype=="KO") )
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
  wt.seurat=subset( tumor, subset=(phenotype=="WT") )
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

  # find DEGs between KO and WT
  # tumor.markers <- FindMarkers(tumor, ident.1 = "KO", ident.2= "WT", min.pct = 0.1, group.by="phenotype", assay="SCT")
  # saveRDS(tumor.markers, "daqingge.1.tumor.markers.KOvsWT.V2.rds")
  tumor.markers=readRDS("daqingge.1.tumor.markers.KOvsWT.V2.rds")
  markers=tumor.markers
  markers=markers[markers$p_val_adj<0.05, ]
  markers=markers[order(markers$avg_logFC, decreasing = TRUE), ]
  markers$FoldChange=exp(markers$avg_logFC)
  write.csv(markers[,c( "FoldChange","pct.1", "pct.2", "p_val_adj")], file=paste("markers.KOvsWT.V2",".csv",sep=""), row.name=T)

  marker.wt.spec =  markers[markers$pct.1<0.05, ]
  deg.wt.spec=rownames(marker.wt.spec)
  marker.ko.spec =  markers[markers$pct.2<0.05, ]
  deg.ko.spec=rownames(marker.ko.spec)

  # eFig.B  :  heatmap of phenotype specific DEG sets
  tocheck=c(deg.wt.spec, deg.ko.spec)
  # tumor@meta.data$phenotype=factor(tumor@meta.data$phenotype, levels=c("WT", "KO"))
  # tumor <- ScaleData(tumor, features = rownames(tumor), assay="RNA" )
  pdf("heatmap.DEG.phenotype.pdf")
  DoHeatmap( tumor, features = tocheck, group.by = "phenotype", slot="scale.data", assay="RNA", angle = 0 ) 
  dev.off()

  str(tumor@meta.data$phenotype)


  enrichment.gene.set(rownames(markers)[markers$avg_logFC>log(1.5)], prefix="ko.DEGs.FC1p5")
  enrichment.gene.set(rownames(markers)[markers$avg_logFC<log(0.5)], prefix="wt.DEGs.FC1p5")

  ko.up.cor.genes=intersect(ko.cor.genes, rownames(markers)[markers$avg_logFC>log(1.5)])
  wt.up.cor.genes=intersect(wt.cor.genes, rownames(markers)[markers$avg_logFC<log(1/1.5)])

  write.table(ko.up.cor.genes, file="ko.up.cor.genes.txt", col.name=F, row.name=F, quote=F, sep='\t')
  write.table(wt.up.cor.genes, file="wt.up.cor.genes.txt", col.name=F, row.name=F, quote=F, sep='\t')
  enrichment.gene.set(ko.up.cor.genes, prefix="ko")
  enrichment.gene.set(wt.up.cor.genes, prefix="wt")

  # eFig.G  :  heatmap of pseudotime-correlated & phenotype-related DEG sets
  # order as pseudotime
  pseudotime=tumor@meta.data$time
  names(pseudotime)=rownames(tumor@meta.data)
  pseudotime=sort(pseudotime)
  tocheck=c(ko.up.cor.genes, wt.up.cor.genes)
  pdf("heatmap.DEG.as.pseudotime.phenotype.pdf", width=8, height=14)
  DoHeatmap( tumor, features = tocheck, cells= names(pseudotime), group.by = "phenotype", slot="scale.data", assay="RNA", angle = 0 ) # 
  dev.off()


  tocheck=c(ko.up.cor.genes, wt.up.cor.genes)
  # tocheck=c("ITGA2")
  for(i in tocheck)
  {
    pdf(paste("vlnplot.", i, ".pdf", sep=""), width=5, height=5)
    temp=VlnPlot(tumor, features = i, group.by="pheno.time.group",assay = "RNA", pt.size=0, cols=c("#99CCFF", "#9999FF", "#000099", "#FFCCCC", "#FF6666", "#660000"))
    temp=temp+stat_summary(fun.y = mean, geom='point', size = 8, colour = "black", shape = 95)
    print(temp)
    dev.off()
  }

  # Fig.C
  tocheck=c("TWIST1", "CDH2", "EDIL3", "HES1", "NNMT", "APBB2", "COL5A1", "HSPG2", "MATN2", "SEMA3A", "SEMA3C", "BMPR2", "CCND2", "FOXN3", "FSTL1", "BST2", "PODXL", "CAV1")
  # tocheck=c("ITGA2")
  for(i in tocheck)
  {
    pdf(paste("vlnplot.", i, ".pdf", sep=""), width=5, height=3)
    temp=VlnPlot(tumor, features = i, group.by="pheno.time.group",assay = "RNA", pt.size=0, cols=c("#99CCFF", "#9999FF", "#000099", "#FFCCCC", "#FF6666", "#660000"))
    temp=temp+stat_summary(fun.y = mean, geom='point', size = 8, colour = "black", shape = 95)
    print(temp)
    dev.off()
  }

  # Fig.D
  tocheck=c("TWIST1", "CDH2", "EDIL3", "HES1", "NNMT", "APBB2", "COL5A1", "HSPG2", "MATN2", "SEMA3A", "SEMA3C", "BMPR2", "CCND2", "FOXN3", "FSTL1", "BST2", "PODXL", "CAV1")
  for(i in tocheck)
  {
    pdf(paste("featureplot.", i, ".pdf", sep=""), width=5, height=5)
    temp=FeaturePlot(tumor, features = i, cols=c("lightblue", "red"))
    print(temp)
    dev.off()
  }

  # eFig.H
  tocheck=c("TWIST1", "CDH2", "EDIL3", "HES1", "NNMT", "APBB2", "COL5A1", "HSPG2", "MATN2", "SEMA3A", "SEMA3C", "BMPR2", "CCND2", "FOXN3", "FSTL1", "BST2", "PODXL", "CAV1")
  for(i in tocheck)
  {
    pdf(paste("featureplot.split.", i, ".pdf", sep=""), width=10, height=5)
    temp=FeaturePlot(tumor, features = i, cols=c("lightblue", "red"), split.by="phenotype")
    print(temp)
    dev.off()
  }



  # check some genes
  wt.cor["CDH2"]
  "ITGA2"%in%ko.cor.genes
  "COL5A1"%in%ko.cor.genes
  markers["COL5A1", ]
### <<<--- identify phenotype-time DEGs

# check functions across pheno.time.group by qusage
  packageVersion("qusage")
  # [1] '2.24.0'
  require(qusage)

  require(msigdbr)
  m_df = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
  m_t2g = m_df %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
  m_df.reactome = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
  m_t2g.reactome = m_df.reactome %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
  m_df.kegg = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
  m_t2g.kegg = m_df.kegg %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
  m_df.hall = msigdbr(species = "Homo sapiens", category = "H")
  m_t2g.hall = m_df.hall %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
  m_df.regulate = msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:GTRD")
  m_t2g.regulate = m_df.regulate %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()

  # convert annotation from data.frame to list format
  kegg.list=list()
  kegg.name=unique(m_t2g.kegg$gs_name)
  for(i in kegg.name)
  {
    kegg.list[[i]]=m_t2g.kegg$gene_symbol[m_t2g.kegg$gs_name==i]
  }

  hall.list=list()
  hall.name=unique(m_t2g.hall$gs_name)
  for(i in hall.name)
  {
    hall.list[[i]]=m_t2g.hall$gene_symbol[m_t2g.hall$gs_name==i]
  }

  regulate.list=list()
  regulate.name=unique(m_t2g.regulate$gs_name)
  for(i in regulate.name)
  {
    regulate.list[[i]]=m_t2g.regulate$gene_symbol[m_t2g.regulate$gs_name==i]
  }

  tumor@meta.data$for_temp_color=as.numeric(tumor@meta.data$pheno.time.group)
  run_qusage_heatmap.seurat3(tumor, nm = 'hall', hall.list, my.seed=100)
  run_qusage_heatmap.seurat3(tumor, nm = 'kegg', kegg.list, my.seed=100)
  run_qusage_heatmap.seurat3(tumor, nm = 'regulate', regulate.list, my.seed=100)


### --->>> compare cluster 3 and 7 (the 2 origining clusters)
  tumor <- SetIdent(tumor, value = "fig.cluster")
  tumor.markers <- FindMarkers(tumor, ident.1 = "3", ident.2= "7", min.pct = 0.1, group.by="fig.cluster", assay="SCT") # , test.use="DESeq2" 
  # saveRDS(tumor.markers, "daqingge.1.tumor.markers.3vs7.rds")
  # tumor.markers=readRDS("daqingge.1.tumor.markers.3vs7.rds")
  markers=tumor.markers
  markers=markers[markers$p_val_adj<0.05, ]
  markers=markers[order(markers$avg_logFC,decreasing = TRUE),]
  markers$FoldChange=exp(markers$avg_logFC)
  write.csv(markers[,c( "FoldChange","pct.1", "pct.2", "p_val_adj")], file=paste("markers.3vs7.all",".csv",sep=""), row.name=T)

  # heatmap
  markers$gene=rownames(markers)
  top10 <- markers  %>% top_n(n = 20, wt = abs(avg_logFC) )
  DefaultAssay(tumor)="RNA"
  pdf("heatmap.cluster3vs7.pdf")
  DoHeatmap( ScaleData( subset( tumor, subset=(fig.cluster%in%c("3", "7")))), features = top10$gene, group.by = "fig.cluster", slot="scale.data", assay="RNA", angle = 0 ) 
  dev.off()
  top10[1:15,]

  # enrichment
  enrichment.gene.set(markers$gene[markers$avg_logFC>log(1.5)], prefix="3vs7.Up")
  enrichment.gene.set(markers$gene[markers$avg_logFC<log(0.66)], prefix="3vs7.Down")

  enrichment.gene.set.plot(markers$gene[markers$avg_logFC<log(0.66)])

### <<<--- compare cluster 3 and 7 (the 2 origining clusters)


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

### check cell lineage across phenotype
{
  ### CCA to identify epithelial cell types WithiN luminal tumor cells
  # norme=readRDS("yulu.normal.bladder.seurat.rds")
  # tumor=readRDS("daqingge.2.tumor.SCTpure.cycleRgrsOut.V2.rds")
  DefaultAssay(norme)="SCT"
  DefaultAssay(tumor)="SCT"
  anchors <- FindTransferAnchors(reference = norme, query = tumor, 
      dims = 1:30, normalization.method="LogNormalize", reference.assay="SCT")
  predictions <- TransferData(anchorset = anchors, refdata = norme$celltype, 
      dims = 1:30)
  tumor@meta.data$predtype=predictions$predicted.id 
  pdf("dimplot.predtype.pdf", height=5, width=6)
  DimPlot(tumor, label = FALSE, group.by="predtype", raster = TRUE)
  dev.off()
  table(norme@meta.data$celltype)

  tocheck=c("KRT5", "KRT17",     "KRT13",     "TNNT1",     "UPK1A", "UPK2")
  # tumor@meta.data$phenotype=factor(tumor@meta.data$phenotype, levels=c("WT", "KO"))
  # tumor <- ScaleData(tumor, features = rownames(tumor), assay="RNA" )
  pdf("heatmap.DEG.phenotype.pdf")
  DoHeatmap( tumor, features = tocheck, group.by = "predtype", slot="scale.data", assay="RNA", angle = 0 ) 
  dev.off()

  ### check predicted cell lineage percentage in phenotypes
  {
    # celltype.sta=table(incidente@meta.data$epitype.pheno, incidente@meta.data$fig.zone.side.p)
    celltype.sta=table(tumor@meta.data$predtype, tumor@meta.data$phenotype)
    celltype.sta=apply(celltype.sta, 2, function(x) 100*x/sum(x))
    print(unique(colnames(celltype.sta)))
    celltype.sta=celltype.sta[, c( "WT", "KO" ) ]
    barN=ncol(celltype.sta)
    n_top = nrow(celltype.sta) 
    pdf(paste("barplot.celltypeAcrossPheno.pdf",sep=""))
    # the width between paper edge and plot ( down,left,up,right : four directions)
    par(mar = c(4, 4, 1, 12), xpd=T)
    # coll=DiscretePalette(n_top, palette = "polychrome")
    require("scales")
    coll=hue_pal()(n_top)
    barplot( celltype.sta, xlab="", ylab="Percentage", col=coll, names.arg=colnames(celltype.sta), las=2, cex.names=1) # , xaxt="n"
    # facts=barN+1
    end_point = barN #0.5 + n_top *facts-1
    # text(x=(c(1:barN)*1.2-1-0.35), y=0-0.35, cex=1, labels=colnames(celltype.sta), srt=45) # par("usr")[3] - 0, xpd=TRUE
    legend( barN+0.5,50, rownames(celltype.sta), fill=coll );
    # text( x=seq(1,end_point,by=1)*1.2-0.2, y = celltype.sta[1,]/2, adj= c(1,1), xpd = TRUE, labels = paste(round(celltype.sta[1,],1),"%") , cex=0.7 )
    # text( x=seq(1,end_point,by=1)*1.2-0.2, y = 100-celltype.sta[2,]/2, adj= c(1,1), xpd = TRUE, labels = paste(round(celltype.sta[2,],1),"%") , cex=0.7 )
    dev.off()
  }

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