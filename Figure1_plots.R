### Figure 1 and supplement

library(Seurat)
library(monocle3)
library(dplyr)
library(viridis)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)

setwd("/Volumes/Reiterlab_3/NS_WtAd3_Ad3_Tdsort_agg/")

seur=readRDS("NS_WtAd3_Ad3_Tdsort_v1_magic.rds")
cds=readRDS(file=paste0("v1/NS_WtAd3_Ad3_Tdsort_agg_v1_aligned.rds"))


seur$paper_names <- factor(seur$paper_names, levels=c("Basal","Int","MCC 1","MCC 2","MCC 3","MCC 4","MCC 5","Secretory","G2M","S","Perictye-like","Neuroendocrine"))

DimPlot(seur, reduction="monocle",group.by="paper_names",label=T)

identical(rownames(colData(cds)),rownames(seur@meta.data))

colData(cds) <- cbind(colData(cds),paper_names=seur$paper_names)

plot_cells(cds, color_cells_by = "paper_names")

#### clusters with overlay of pseudotime

mcc_color <- colorRampPalette(c("lightskyblue2", "royalblue4"))

fig1_col = c("#F7A12D","#009245",rep(mcc_color(5)[3],5),"indianred","#996699","palegreen1","#998675","darkgoldenrod")

colData(cds)[,"paper_names"] <- factor(colData(cds)[,"paper_names"],levels=c("Basal","Int","MCC 1","MCC 2","MCC 3","MCC 4","MCC 5","Secretory","G2M","S","Perictye-like","Neuroendocrine"))

pdf("Figure1_v3/monocle_umap_clusters.pdf",height=8,width=10)
plot_cells(cds, 
           color_cells_by = "paper_names",
           show_trajectory_graph = T,
           label_cell_groups = T,
           cell_size = 2,
           cell_stroke = 0,
           alpha = 1,
           labels_per_group = F,
           rasterize = T, 
           label_branch_points = F, 
           label_roots = F, 
           label_leaves = F, 
           trajectory_graph_color = "black",
           trajectory_graph_segment_size = 2)+
        scale_color_manual(values=fig1_col)+
        theme_void()+
        theme(legend.position = "none")
dev.off()

##### split by dataset

pdf("Figure1_v3/NS_WtAd3_Ad3_Tdsort_split_UMAP.pdf",height=4,width=15,useDingbats = F)
DimPlot(seur, reduction="monocle",group.by="paper_names", cols = fig1_col, pt.size = 2,raster = T,split.by = "orig.ident") + theme_void() + NoLegend()
dev.off()

##### EXPRESSION MAPS
DefaultAssay(seur) <- "RNA"

pdf("Figure1_v3/Expression_maps_Deup1_Dnah5_Ascl1_Mki67.pdf",height=4,width=18,useDingbats = F)
FeaturePlot(seur, c("Deup1","Dnah5","Ascl1","Mki67"), reduction="monocle",order=T, raster=T,pt.size=1.5,ncol=4) & scale_color_viridis() & theme_void()
dev.off()


###### overlay of three gradients
library(dittoSeq)
library(ggnewscale)
scgb=dittoDimPlot(seur, "Scgb3a1",reduction="monocle",assay = "RNA",data.out = T)
trp63=dittoDimPlot(seur,"Trp63",reduction="monocle",assay="RNA",data.out = T)
foxj1=dittoDimPlot(seur,"Foxj1",reduction="monocle",assay="RNA",data.out = T)
myb=dittoDimPlot(seur,"Myb",reduction="monocle",assay="RNA",data.out = T)

identical(scgb$Target_data$X,trp63$Target_data$X)

full_dat = scgb$Target_data
colnames(full_dat) <- c("X","Y","gene_Scgb3a1")
full_dat$gene_Trp63 = trp63$Target_data$color
full_dat$gene_Foxj1 = foxj1$Target_data$color
full_dat$gene_Myb = myb$Target_data$color

full_dat$bin_scgb3a1 = 0
full_dat$bin_scgb3a1[full_dat$gene_Scgb3a1 > 0.25] = 1
full_dat$bin_trp63 = 0
full_dat$bin_trp63[full_dat$gene_Trp63>0.25] = 1
full_dat$bin_foxj1 = 0
full_dat$bin_foxj1[full_dat$gene_Foxj1 > 0.25 ] = 1
full_dat$bin_myb = 0
full_dat$bin_myb[full_dat$gene_Myb > 0.25 ] = 1


c("#F7A12D","#009245","#6589BC","indianred")

## cap top 10% expression levels so easier to see lower gradient
trp63_top=range(full_dat$gene_Trp63)[2] - range(full_dat$gene_Trp63)[2]*0.1
full_dat$gene_Trp63[full_dat$gene_Trp63 > trp63_top] = trp63_top

foxj1_top=range(full_dat$gene_Foxj1)[2] - range(full_dat$gene_Foxj1)[2]*0.1
full_dat$gene_Foxj1[full_dat$gene_Foxj1 > foxj1_top] = foxj1_top

scgb3a1_top=range(full_dat$gene_Scgb3a1)[2] - range(full_dat$gene_Scgb3a1)[2]*0.1
full_dat$gene_Scgb3a1[full_dat$gene_Scgb3a1 > scgb3a1_top] = scgb3a1_top

myb_top=range(full_dat$gene_Myb)[2] - range(full_dat$gene_Myb)[2]*0.1
full_dat$gene_Myb[full_dat$gene_Myb > myb_top] = myb_top

g1=ggplot(full_dat, aes(x=X,y=Y,col=gene_Trp63))+ggrastr::rasterise(geom_point(size=0.3))+scale_color_gradient(low="grey20",high="#F7A12D") + theme_void()
g2=ggplot(full_dat, aes(x=X,y=Y,col=gene_Scgb3a1))+ggrastr::rasterise(geom_point(size=0.3))+scale_color_gradient(low="grey20",high="indianred") + theme_void()
g3=ggplot(full_dat, aes(x=X,y=Y,col=gene_Foxj1))+ggrastr::rasterise(geom_point(size=0.3))+scale_color_gradient(low="grey20",high="#6589BC") + theme_void()
g4=ggplot(full_dat, aes(x=X,y=Y,col=gene_Myb))+ggrastr::rasterise(geom_point(size=0.3))+scale_color_gradient(low="grey20",high=mcc_color(5)[3]) + theme_void()

pdf("Figure1_v3/Trp63_Scgb3a1_Myb_Foxj1_gradients.pdf",height=6,width=9,useDingbats = F)
ggarrange(g1,g2,g4,g3)
dev.off()

##### heatmap of top DE genes
Idents(seur) <- "paper_names"

seur$combined_mcc <- as.character(seur$paper_names)
seur$combined_mcc[seur$paper_names == "MCC 1"] <- "MCC"
seur$combined_mcc[seur$paper_names == "MCC 2"] <- "MCC"
seur$combined_mcc[seur$paper_names == "MCC 3"] <- "MCC"
seur$combined_mcc[seur$paper_names == "MCC 4"] <- "MCC"
seur$combined_mcc[seur$paper_names == "MCC 5"] <- "MCC"

seur$combined_mcc <- factor(seur$combined_mcc, levels=c("Basal","Int","Secretory","MCC","G2M","S","Perictye-like","Neuroendocrine"))

DimPlot(seur, group.by="combined_mcc",reduction = "monocle",label=T)
DefaultAssay(seur) <- "RNA"
Idents(seur) <- "combined_mcc"

marks = FindAllMarkers(seur, assay="RNA",group.by="combined_mcc", only.pos = T)
write.csv(marks, file="Figure1_v3/NS_WtAd3_Ad3_Tdsort_combined_names_markers.csv")

marks_filt = marks %>% filter(p_val_adj < 0.05) %>% filter(avg_log2FC > log2(1.5))

marks_filt$cluster <- factor(marks_filt$cluster, levels=c("Basal","Int","Secretory","MCC","G2M","S","Perictye-like","Neuroendocrine"))


top5 = marks_filt %>% group_by(cluster) %>% slice_max(order_by=avg_log2FC, n=5)


DefaultAssay(seur) <- "RNA"
seur <- ScaleData(seur,features=rownames(seur[["RNA"]]@data))

heat_col = c("#F7A12D","#009245","palevioletred2",mcc_color(5)[3],"#996699","palegreen1","#998675","darkgoldenrod")


pdf("Figure1_v3/combined_mcc_clusters_heatmap.pdf",height=8,width=11,useDingbats = F)
DoHeatmap(subset(seur, downsample = 7500), features =unique(top5$gene),assay = "RNA",raster = T,group.by="combined_mcc",group.colors = heat_col) + scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()


#### label transfer
DimPlot(seur, reduction="monocle",group.by="predicted.celltype")

###################### zaragosi ####################

preds =read.csv("label_transfer/NS_WtAd3_Ad3_Tdsort_v1_magic_zaragosi_mouseintegrated_predictions.csv",row.names=1)
seur$predicted_zaragosi <-preds$predicted.id

zara_col = c("#F7A12D","#996699",mcc_color(5)[c(2,5)],"#998675","grey20","palevioletred2","hotpink4","#009245")

pdf("Figure1_v3/label_transfer_zaragosi.pdf",height=8,width=11,useDingbats = F)
DimPlot(seur, reduction="monocle",group.by="predicted_zaragosi", cols = zara_col, pt.size = 2,raster=T) + theme_void()
dev.off()


#################### Notch expression plots ####################
DefaultAssay(seur) <- "RNA"

FeaturePlot(seur, c("Heyl","Hes1","Notch3","Notch4"), reduction="monocle")

pdf("Figure1_v3/notch_genes.pdf",height=8,width=7,useDingbats = F)
FeaturePlot(seur, c("Jag2","Dll1","Notch1","Notch2","Hey1",'Hes1'),reduction="monocle", order=T,raster = T,pt.size=2) & scale_color_viridis() & theme_void()
dev.off()

############### basal to intermed v2 reclustering ############
subseur2= readRDS("subclus_basal_to_intermed_v2_clean.rds")

library(dittoSeq)
mcc_color <- colorRampPalette(c("lightskyblue2", "royalblue4"))
cols=dittoColors(reps=1)[1:9]
cols_2=c("khaki2","tan2","darkturquoise",cols[7],"indianred",mcc_color(5)[4])

pdf("Figure1_v3//basal_to_intermed_v2_umap.pdf",height=8,width=11,useDingbats = F)
DimPlot(subseur2, reduction="umap",group.by="ord_clus", cols=cols_2, pt.size=2,raster=T) + theme_void() + NoLegend()
dev.off()

####### EXPRESSION PLOTS #########

pdf("Figure1_v3/basal_to_intermed_expressionmap.pdf",height=8,width=10,useDingbats = F)
FeaturePlot(subseur2, c("Myb","Scgb3a1","Mycl","Hey1"),reduction="umap",ncol=2,order=T,raster=T,pt.size=2.5,max.cutoff = "q99") & scale_color_viridis() & theme_void()
dev.off()


### Hey1 and Mycl correlation

DefaultAssay(subseur2) <- "RNA"
matrix<-subseur2@assays$RNA@data
matrix_mod<-as.matrix(matrix)

gene<-as.numeric(matrix_mod["Mycl",])

mycl_correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})
mycl_correlations_noNA <- mycl_correlations[!is.na(mycl_correlations)]
mycl_ordered_corr <- mycl_correlations_noNA[order(mycl_correlations_noNA,decreasing = T)]

gene2<-as.numeric(matrix_mod["Hey1",])
scg_correlations<-apply(matrix_mod,1,function(x){cor(gene2,x)})
identical(names(mycl_correlations),names(scg_correlations))
scg_correlations_noNA <- scg_correlations[!is.na(scg_correlations)]
scg_ordered_corr <- scg_correlations_noNA[order(scg_correlations_noNA,decreasing = T)]

df_full = cbind(mycl_correlations,scg_correlations)
colnames(df_full) <- c("Mycl","Hey1")

df_full_thresh = as.data.frame(df_full) %>% filter(Mycl > 0.1 | Hey1 > 0.1)
df_full_thresh = df_full_thresh %>% filter(Mycl - Hey1 > 0.1 | Hey1 - Mycl > 0.1)

df_full_thresh = df_full_thresh %>% arrange(desc(Mycl))

mycl_genes = df_full_thresh %>% filter(Mycl > 0.1 & Hey1 < 0.1)
scg_genes = df_full_thresh %>% filter(Hey1 > 0.1 & Mycl < 0.1)

df_clear = df_full_thresh[match(c(rownames(mycl_genes),rownames(scg_genes)),rownames(df_full_thresh)),]

write.csv(df_clear, "basal_to_intermed_v2/basal_to_intermed_v2_clean_Mycl_vs_Hey1_correlated_genes_thresh0.1.csv")
df_clear=read.csv("basal_to_intermed_v2/basal_to_intermed_v2_clean_Mycl_vs_Hey1_correlated_genes_thresh0.1.csv",row.names=1)

df_clear_thresh = df_clear[df_clear$Mycl > 0.15 | df_clear$Hey1 > 0.15,]

mycl_arrange = df_clear_thresh %>% filter(Mycl > 0.15) %>% arrange(.,-Mycl)
hey1_arrange = df_clear_thresh %>% filter(Hey1 > 0.15) %>% arrange(.,-Hey1)

df_top=rbind(mycl_arrange, hey1_arrange)

write.csv(df_top, file="basal_to_intermed_v2/basal_to_intermed_v2_clean_Mycl_vs_Hey1_correlated_genes_thresh0.15_heatmap.csv")

paletteLength=50
cols_heatmap= colorRampPalette(brewer.pal(11,"RdBu"))
myBreaks = c(seq(-1, 0, length.out=ceiling(paletteLength/2) + 1), 
             seq(max(df_top)/paletteLength, max(df_top), length.out=floor(paletteLength/2)))

pdf("basal_to_intermed_v2/basal_to_intermed_v2_heatmap_corr_mycl_vs_hey1.pdf",height=11,width=8,useDingbats = F)
pheatmap(df_top,cluster_cols = F,cluster_rows = F,color = rev(cols_heatmap(paletteLength)),breaks=myBreaks)
dev.off()


pdf("basal_to_intermed_v2/basal_to_intermed_v2_heatmap_corr_mycl_vs_hey1_myclonly.pdf",height=11,width=8,useDingbats = F)
pheatmap(df_top[1:15,],cluster_cols = F,cluster_rows = F,color = rev(cols_heatmap(paletteLength)),breaks=myBreaks)
dev.off()

pdf("basal_to_intermed_v2/basal_to_intermed_v2_heatmap_corr_mycl_vs_hey1_hey1only.pdf",height=11,width=8,useDingbats = F)
pheatmap(df_top[51:65,],cluster_cols = F,cluster_rows = F,color = rev(cols_heatmap(paletteLength)),breaks=myBreaks)
dev.off()
