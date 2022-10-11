#### Subcluster basal and intermediate cells from integrated ALID03 dataset
setwd("/Volumes/Reiterlab_3/NS_WtAd3_Ad3_Tdsort_agg/")
library(monocle3)
library(Seurat)
library(viridis)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(dplyr)

seur=readRDS("NS_WtAd3_Ad3_Tdsort_v1_magic.rds")
dir.create("basal_to_intermed_v2")

DimPlot(seur, reduction="monocle",group.by="monocle_ord")

seur$paper_names <- seur$monocle_ord
seur$paper_names[seur$monocle_ord == 9] <- 2

DimPlot(seur, reduction="monocle",group.by="paper_names")
Idents(seur) <- "paper_names"
subseur = subset(seur, idents=c(1,2))
DimPlot(subseur, reduction="monocle")

DefaultAssay(subseur) <- "RNA"
subseur=FindVariableFeatures(subseur,assay = "RNA")

DefaultAssay(subseur) <- "integrated"
subseur=ScaleData(subseur)
subseur=RunPCA(subseur)
ElbowPlot(subseur,ndims=40)
subseur=FindNeighbors(subseur,assay = "integrated",dims = 1:20)
subseur=FindClusters(subseur,resolution = 0.5)
subseur=RunUMAP(subseur,dims = 1:20)

DimPlot(subseur, group.by="orig.ident")
DimPlot(subseur,label=T, reduction="umap", group.by="integrated_snn_res.0.5")

DefaultAssay(subseur) <- "RNA"

### rename res 0.5 clusters

subseur$res5_clus <- NA
subseur$res5_clus[subseur$integrated_snn_res.0.5 == 1] <- 1
subseur$res5_clus[subseur$integrated_snn_res.0.5 == 2] <- 2
subseur$res5_clus[subseur$integrated_snn_res.0.5 == 3] <- 3
subseur$res5_clus[subseur$integrated_snn_res.0.5 == 0] <- 4
subseur$res5_clus[subseur$integrated_snn_res.0.5 == 4] <- 5
subseur$res5_clus[subseur$integrated_snn_res.0.5 == 6] <- 6
subseur$res5_clus[subseur$integrated_snn_res.0.5 == 7] <- 7
subseur$res5_clus[subseur$integrated_snn_res.0.5 == 8] <- 8
subseur$res5_clus[subseur$integrated_snn_res.0.5 == 9] <- 9
subseur$res5_clus[subseur$integrated_snn_res.0.5 == 10] <- 10
subseur$res5_clus[subseur$integrated_snn_res.0.5 == 11] <- 11
subseur$res5_clus[subseur$integrated_snn_res.0.5 == 5] <- 12

DimPlot(subseur, group.by="res5_clus",label=T,reduction="umap")
Idents(subseur) <- "res5_clus"

marks = FindAllMarkers(subseur, group.by="res5_clus",only.pos = T,assay="RNA")

marks_filt =marks %>% filter(p_val_adj < 0.05)

write.csv(marks_filt, file="basal_to_intermed_v2/subclus_basal_to_intermed_v2_res5_markers.csv")

saveRDS(subseur, file="basal_to_intermed_v2/subcluster_basal_to_intermed_v2.rds")
subseur=readRDS("basal_to_intermed_v2/subcluster_basal_to_intermed_v2.rds")

FeaturePlot(subseur, c("Agr2","Dynlrb2","Foxj1","Muc5b","Mycl","Myb"),reduction="umap",order=T) &scale_color_viridis()

DimPlot(subseur, label=T,reduction="umap")

#### keep just main body cells
Idents(subseur) <- "res5_clus"
subseur2 = subset(subseur, idents=c(1:5))
DimPlot(subseur2,reduction="umap")

DefaultAssay(subseur2) <- "RNA"
subseur2=FindVariableFeatures(subseur2,assay = "RNA")

DefaultAssay(subseur2) <- "integrated"
subseur2=ScaleData(subseur2)
subseur2=RunPCA(subseur2)
ElbowPlot(subseur2,ndims=40)
subseur2=FindNeighbors(subseur2,assay = "integrated",dims = 1:20)
subseur2=FindClusters(subseur2,resolution = 0.4)
subseur2=RunUMAP(subseur2,dims = 1:20)

DimPlot(subseur2, group.by="orig.ident")
DimPlot(subseur2,label=T, reduction="umap", group.by="integrated_snn_res.0.4")


DefaultAssay(subseur2) <- "RNA"


## rename clusters

subseur2$ord_clus <- NA
subseur2$ord_clus[subseur2$integrated_snn_res.0.4 == 1] <- 1
subseur2$ord_clus[subseur2$integrated_snn_res.0.4 == 2] <- 2
subseur2$ord_clus[subseur2$integrated_snn_res.0.4 == 3] <- 3
subseur2$ord_clus[subseur2$integrated_snn_res.0.4 == 0] <- 4
subseur2$ord_clus[subseur2$integrated_snn_res.0.4 == 4] <- 5


DimPlot(subseur2, reduction="umap",group.by="ord_clus")

Idents(subseur2)<-"ord_clus"

marks = FindAllMarkers(subseur2, group.by="ord_clus",only.pos = T,assay="RNA")

marks_filt =marks %>% filter(p_val_adj < 0.05)

write.csv(marks_filt, file="basal_to_intermed_v2/basal_to_intermed_clean_res0.4_marks.csv")
marks_filt=read.csv(file="basal_to_intermed_v2/basal_to_intermed_clean_res0.4_marks.csv",row.names=1)

saveRDS(subseur2, file="subclus_basal_to_intermed_v2_clean.rds")
subseur2= readRDS("subclus_basal_to_intermed_v2_clean.rds")

