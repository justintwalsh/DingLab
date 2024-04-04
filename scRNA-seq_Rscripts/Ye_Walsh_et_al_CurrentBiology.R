#This script will allow you to recreate Figure3 A, B, and F 
#from Ye, Walsh, et al. 2024 
#Changes in the cellular makeup of motor patterning circuits drive courtship song evolution in Drosophila  

#Load required packages
library(Seurat)
library(SeuratDisk)
library(ggplot2)

#Get data from Ding lab Github or GEO Accession #GSE262732
#Two Rds files, one for just TN1 neurons (Figure3 A and B)
#and one for the entire dsx dataset (Figure3 F)

#Set working directory to file location
#setwd("")
#Load in TN1 dataset
TN1<-LoadSeuratRds("TN1OnlyMelYak.Rds")

#Jackstraw analysis to see how many PCs to include in cluster analysis
DefaultAssay(TN1) <- "integrated"
TN1 <- RunPCA(TN1, npcs = 50, verbose = FALSE) 
TN1 <- JackStraw(TN1, dims=50)
TN1 <- ScoreJackStraw(TN1, dims = 1:50)
ElbowPlot(TN1, ndims=50)
JackStrawPlot(TN1, dims = 1:50)
#PCs one through seven explain a significant amount of variation

#Clustering analysis using 7 PCs and a resolution of 0.4
#Supplemental Figure4 explores the resiliency of the subclusters across different PCs and resolutions

TN1 <- RunPCA(TN1, npcs = 15, verbose = FALSE)
TN1 <- RunUMAP(TN1, reduction = "pca", dims = 1:7)
TN1 <- FindNeighbors(TN1, reduction = "pca", dims = 1:7)
TN1 <- FindClusters(TN1, resolution = 0.4)
#5 subclusters

#Plots
#Show both species together
DimPlot(TN1, reduction = "umap", label=T, label.size=10,pt.size=1)+NoLegend()

#Show species in same space, ID by color
DimPlot(TN1, reduction = "umap", group.by = "Species",pt.size=1, label=F)

#Show species on separate UMAP spaces (Figure3 A)
DimPlot(TN1, reduction = "umap", label=F,pt.size=1, split.by="Species",label.size=8)+NoLegend()

#Find markers for each TN1 subcluster
TN1markers<-FindAllMarkers(TN1,assay="RNA", only.pos = T,logfc.threshold = 0, min.pct = 0)

#CCKLR-17D3 is a great marker for the melanogaster-specific cluster 4
#plot expression of CCKLR-17D3
DefaultAssay(TN1) <- "RNA"
#Combined
FeaturePlot(TN1, features = c("CCKLR-17D3"),pt.size=1, reduction = "umap", label=T, label.size=5)

#Split by species (Figure3 B)
FeaturePlot(TN1, features = c("CCKLR-17D3"),pt.size=1, split.by="Species", reduction = "umap", label=F, label.size=5)


#Look at expression of various genes of interest (SFigure4 B)
FeaturePlot(TN1, features = c("elav","ChAT","dsx","VAChT","fru","Gad1","TfAP-2","VGAT","Antp","VGlut"),pt.size=0.1,reduction = "umap", label=F, label.size=2, split.by="Species",keep.scale = "all") + theme(legend.position = "right")

#Check for robustness of mel-specific subcluster across different parameters (SFigure4 A)

TN1 <- RunPCA(TN1, npcs = 15, verbose = FALSE)
TN1 <- RunUMAP(TN1, reduction = "pca", dims = 1:5)
TN1 <- FindNeighbors(TN1, reduction = "pca", dims = 1:5)
TN1 <- FindClusters(TN1, resolution = 0.2)

p5_02<-DimPlot(TN1, reduction = "umap", label=F,pt.size=0.1, split.by="Species",label.size=8)+NoLegend()

TN1 <- RunPCA(TN1, npcs = 15, verbose = FALSE)
TN1 <- RunUMAP(TN1, reduction = "pca", dims = 1:5)
TN1 <- FindNeighbors(TN1, reduction = "pca", dims = 1:5)
TN1 <- FindClusters(TN1, resolution = 0.4)

p5_04<-DimPlot(TN1, reduction = "umap", label=F,pt.size=0.1, split.by="Species",label.size=8)+NoLegend()

TN1 <- RunPCA(TN1, npcs = 15, verbose = FALSE)
TN1 <- RunUMAP(TN1, reduction = "pca", dims = 1:5)
TN1 <- FindNeighbors(TN1, reduction = "pca", dims = 1:5)
TN1 <- FindClusters(TN1, resolution = 0.6)

p5_06<-DimPlot(TN1, reduction = "umap", label=F,pt.size=0.1, split.by="Species",label.size=8)+NoLegend()

TN1 <- RunPCA(TN1, npcs = 15, verbose = FALSE)
TN1 <- RunUMAP(TN1, reduction = "pca", dims = 1:7)
TN1 <- FindNeighbors(TN1, reduction = "pca", dims = 1:7)
TN1 <- FindClusters(TN1, resolution = 0.2)

p7_02<-DimPlot(TN1, reduction = "umap", label=F,pt.size=0.1, split.by="Species",label.size=8)+NoLegend()

TN1 <- RunPCA(TN1, npcs = 15, verbose = FALSE)
TN1 <- RunUMAP(TN1, reduction = "pca", dims = 1:7)
TN1 <- FindNeighbors(TN1, reduction = "pca", dims = 1:7)
TN1 <- FindClusters(TN1, resolution = 0.4)

p7_04<-DimPlot(TN1, reduction = "umap", label=F,pt.size=0.1, split.by="Species",label.size=8)+NoLegend()

TN1 <- RunPCA(TN1, npcs = 15, verbose = FALSE)
TN1 <- RunUMAP(TN1, reduction = "pca", dims = 1:7)
TN1 <- FindNeighbors(TN1, reduction = "pca", dims = 1:7)
TN1 <- FindClusters(TN1, resolution = 0.6)

p7_06<-DimPlot(TN1, reduction = "umap", label=F,pt.size=0.1, split.by="Species",label.size=8)+NoLegend()

TN1 <- RunPCA(TN1, npcs = 15, verbose = FALSE)
TN1 <- RunUMAP(TN1, reduction = "pca", dims = 1:9)
TN1 <- FindNeighbors(TN1, reduction = "pca", dims = 1:9)
TN1 <- FindClusters(TN1, resolution = 0.2)

p9_02<-DimPlot(TN1, reduction = "umap", label=F,pt.size=0.1, split.by="Species",label.size=8)+NoLegend()

TN1 <- RunPCA(TN1, npcs = 15, verbose = FALSE)
TN1 <- RunUMAP(TN1, reduction = "pca", dims = 1:9)
TN1 <- FindNeighbors(TN1, reduction = "pca", dims = 1:9)
TN1 <- FindClusters(TN1, resolution = 0.4)

p9_04<-DimPlot(TN1, reduction = "umap", label=F,pt.size=0.1, split.by="Species",label.size=8)+NoLegend()

TN1 <- RunPCA(TN1, npcs = 15, verbose = FALSE)
TN1 <- RunUMAP(TN1, reduction = "pca", dims = 1:9)
TN1 <- FindNeighbors(TN1, reduction = "pca", dims = 1:9)
TN1 <- FindClusters(TN1, resolution = 0.6)

p9_06<-DimPlot(TN1, reduction = "umap", label=F,pt.size=0.1, split.by="Species",label.size=8)+NoLegend()

CombinePlots(plots=list(p5_02,p7_02,p9_02,p5_04,p7_04,p9_04,p5_06,p7_06,p9_06),ncol=3)

#####################
#Load the entire trimmed dsx dataset to look for the presence of the melanogaster-specific cluster
MelYakFinal <- LoadSeuratRds("Global_Trimmed_MelYak.Rds")

#Jackstraw analysis to see how many PCs to include in cluster analysis
DefaultAssay(MelYakFinal) <- "integrated"
MelYakFinal <- RunPCA(MelYakFinal, npcs = 100, verbose = FALSE) 

MelYakFinal <- JackStraw(MelYakFinal, num.replicate = 200, dims=100)
MelYakFinal <- ScoreJackStraw(MelYakFinal, dims = 1:100)
JackStrawPlot(MelYakFinal, dims = 1:100)
#Keep 88 PCs

#Clustering analysis
DefaultAssay(MelYakFinal) <- "integrated"
MelYakFinal <- RunPCA(MelYakFinal, npcs = 100, verbose = FALSE) #
MelYakFinal <- RunUMAP(MelYakFinal, reduction = "pca", dims = 1:88)
MelYakFinal <- RunTSNE(MelYakFinal,dims = 1:88, perplexity=30)
MelYakFinal <- FindNeighbors(MelYakFinal, reduction = "pca", dims = 1:88)
MelYakFinal <- FindClusters(MelYakFinal, resolution = 4,graph.name='integrated_snn')
#97 clusters

#Plot
DimPlot(MelYakFinal, reduction = "umap", label=T,pt.size=1,label.size=5)+NoLegend()
#Where is TN1? TN1 can be identified by the expression of Antp and TfAP-2

DefaultAssay(MelYakFinal) <- "RNA"
FeaturePlot(MelYakFinal, features = c("Antp","TfAP-2"),pt.size=1,reduction = "umap", label=F, label.size=5)
VlnPlot(MelYakFinal, features = c("Antp","TfAP-2"), pt.size=1.5) +NoLegend()
#Clusters 1, 21, and 28 represent TN1

#Plot both species in same space, ID by color (Figure3 F)
DimPlot(MelYakFinal, reduction = "umap", group.by = "Species",pt.size=0.1, cols = c("red","blue"))
#Is the red space in cluster 1, the CCKLR-17D3 cluster?

#(Figure3 F small panel)
FeaturePlot(MelYakFinal, features = c("CCKLR-17D3"),pt.size=1,reduction = "umap", label=F, label.size=5)
#Yes!
