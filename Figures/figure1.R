### Figure 1 and supplement

library(Seurat)
library(monocle3)
library(dplyr)
library(viridis)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(multiClust)
library(reshape2)

setwd("/Volumes/Reiterlab_3/NS_WtAd3_Ad3_Tdsort_agg/")

seur=readRDS("NS_WtAd3_Ad3_Tdsort_v1_magic.rds")
cds=readRDS(file=paste0("v1/NS_WtAd3_Ad3_Tdsort_agg_v1_aligned.rds"))
full_pseudo=pseudotime(cds)


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

g1=ggplot(full_dat %>% arrange(gene_Trp63), aes(x=X,y=Y,col=gene_Trp63))+ggrastr::rasterise(geom_point(size=0.3))+scale_color_gradient(low="grey30",high="#F7A12D") + theme_void()
g2=ggplot(full_dat %>% arrange(gene_Scgb3a1), aes(x=X,y=Y,col=gene_Scgb3a1))+ggrastr::rasterise(geom_point(size=0.3))+scale_color_gradient(low="grey30",high="indianred") + theme_void()
g3=ggplot(full_dat %>% arrange(gene_Foxj1), aes(x=X,y=Y,col=gene_Foxj1))+ggrastr::rasterise(geom_point(size=0.3))+scale_color_gradient(low="grey30",high="#6589BC") + theme_void()
g4=ggplot(full_dat %>% arrange(gene_Myb), aes(x=X,y=Y,col=gene_Myb))+ggrastr::rasterise(geom_point(size=0.3))+scale_color_gradient(low="grey30",high=mcc_color(5)[3]) + theme_void()

pdf("Figure1_v3/Trp63_Scgb3a1_Myb_Foxj1_gradients.pdf",height=6,width=9,useDingbats = F)
ggarrange(g1,g2,g4,g3)
dev.off()

##### Pseudotime Plots #######
cds_unbias=readRDS(file=paste0("v1/basal_to_endmcc_v1_cds_unbias.rds"))
pseudos=pseudotime(cds_unbias)
cds_coldata=colData(cds)
cds_coldata$pseudo_basal_endmcc <-NA
cds_coldata$pseudo_basal_endmcc[match(names(pseudos),rownames(cds_coldata))] <- pseudos

colData(cds) <- cds_coldata
plot_cells(cds, color_cells_by = "pseudo_basal_endmcc")

pseudo_vals=pseudotime(cds_unbias)
sub=subset(seur, cells=names(pseudo_vals))
identical(rownames(sub@meta.data),names(pseudo_vals))
sub$pseudo <- pseudo_vals
sub$Pseudo_bin <- as.numeric(cut_interval(sub$pseudo,20))

avg_exp=AverageExpression(sub, assay="RNA",features=c("Myb","Foxj1","Fgfr1op"),group.by = "Pseudo_bin")
avg_norm = t(apply(avg_exp$RNA, 1, nor.min.max))
colnames(avg_norm) <- paste0("Bin",colnames(avg_norm))
plot_val = melt(avg_norm)


pdf("Figure1_v3//pseudotime_bin_basal_to_mcc_Myb_Foxj1_Fgfr1op.pdf",height=8,width=11,useDingbats = F)
ggplot(plot_val, aes(x=Var2,y=value,col=Var1,group=Var1))+
        geom_line(size=2,stat="identity")+
        theme_classic()+
        xlab("Pseudotime")+
        ylab("Gene Expression")+
        scale_color_manual(values=c("Gold1","Magenta","Cyan3"))
dev.off()

### pseudotime just for MCCs

cds_two=choose_graph_segments(cds,clear_cds = FALSE)
cds_two=order_cells(cds_two)
plot_cells(cds_two, color_cells_by = "pseudotime")
pseudo_mcc =pseudotime(cds_two)

sub=subset(seur, cells=names(pseudo_mcc))
identical(rownames(sub@meta.data),names(pseudo_vals))
sub$pseudo <- pseudo_vals
sub$Pseudo_bin <- as.numeric(cut_interval(sub$pseudo,10))

avg_exp=AverageExpression(sub, assay="RNA",features=c("Myb","Foxj1","Fgfr1op"),group.by = "Pseudo_bin")
avg_norm = t(apply(avg_exp$RNA, 1, nor.min.max))
colnames(avg_norm) <- paste0("Bin",colnames(avg_norm))
plot_val = melt(avg_norm)

pdf("Figure1_v3//pseudotime_bin_int_to_mcc_Myb_Foxj1_Fgfr1op.pdf",height=8,width=11,useDingbats = F)
ggplot(plot_val, aes(x=Var2,y=value,col=Var1,group=Var1))+
        geom_line(size=2,stat="identity")+
        theme_classic()+
        xlab("Pseudotime")+
        ylab("Gene Expression")+
        scale_color_manual(values=c("Gold1","Magenta","Cyan3"))
dev.off()

seur$int_to_mcc_pseudo <- NA
seur$int_to_mcc_pseudo[match(names(pseudo_mcc),rownames(seur@meta.data))] <- pseudo_mcc


pdf("Figure1_v3/int_to_mcc_pseudovalue.pdf",height=6,width=9,useDingbats = F)
FeaturePlot(seur, "int_to_mcc_pseudo",reduction="monocle",raster = T) + scale_color_viridis_c(option="inferno",direction = -1)+ theme_void()
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

marks = FindAllMarkers(seur, assay="RNA",group.by="combined_mcc", only.pos = T,test.use = "MAST",latent.vars = "nCount_RNA")
write.csv(marks, file="Figure1_v3/NS_WtAd3_Ad3_Tdsort_combined_names_markers.csv")
marks=read.csv(file="Figure1_v3/NS_WtAd3_Ad3_Tdsort_combined_names_markers.csv")

marks_filt = marks %>% filter(p_val_adj < 0.05) %>% filter(avg_log2FC > log2(1.5))

marks_filt$cluster <- factor(marks_filt$cluster, levels=c("Basal","Int","Secretory","MCC","G2M","S","Perictye-like","Neuroendocrine"))


top5 = marks_filt %>% group_by(cluster) %>% slice_max(order_by=avg_log2FC, n=5)


DefaultAssay(seur) <- "RNA"
seur <- ScaleData(seur,features=rownames(seur[["RNA"]]@data))

heat_col = c("#F7A12D","#009245","palevioletred2",mcc_color(5)[3],"#996699","palegreen1","#998675","darkgoldenrod")


pdf("Figure1_v3/combined_mcc_clusters_heatmap.pdf",height=8,width=11,useDingbats = F)
DoHeatmap(subset(seur, downsample = 7500), features =unique(top5$gene),assay = "RNA",raster = T,group.by="combined_mcc",group.colors = heat_col) + scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()

##### subset heatmap for main clusters #####
Idents(seur) <- "combined_mcc"

sub = subset(seur, idents = c("Basal","Int","Secretory","MCC"))
marks_sub = marks %>% filter(cluster %in% c("Basal","Int","Secretory","MCC"))
marks_sub$cluster <- factor(marks_sub$cluster, levels=c("Basal","Int","Secretory","MCC"))
top5 = marks_sub %>% group_by(cluster) %>% slice_max(order_by=avg_log2FC, n=5)

sub$combined_mcc <- factor(sub$combined_mcc, levels=c("Basal","Int","Secretory","MCC"))

pdf("Figure1_v3/reduced_clusters_heatmap.pdf",height=8,width=11,useDingbats = F)
DoHeatmap(subset(sub, downsample = 7500), features =unique(top5$gene),assay = "RNA",raster = T,group.by="combined_mcc",group.colors = heat_col) + scale_fill_gradientn(colors = c("blue", "white", "red"))
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
cols_2=c("khaki2","tan2","darkturquoise","indianred","indianred") ## combine two secretory clusters
##cols[7] = "#CC79A7"

pdf("Figure1_v3//basal_to_intermed_v2_umap.pdf",height=8,width=11,useDingbats = F)
DimPlot(subseur2, reduction="umap",group.by="ord_clus", cols=cols_2, pt.size=2,raster=T) + theme_void() + NoLegend()
dev.off()

####### EXPRESSION PLOTS #########

DefaultAssay(subseur2) <- "RNA"

pdf("Figure1_v3/basal_to_intermed_expressionmap.pdf",height=10,width=9.5,useDingbats = F)
FeaturePlot(subseur2, c("Trp63","Krt5","Myb","Scgb3a1","Mycl","Hey1"),reduction="umap",ncol=2,order=T,raster=T,pt.size=2.5,max.cutoff = "q99") & scale_color_viridis() & theme_void()
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

### DE expression

DimPlot(subseur2,label=T,reduction="umap")

subseur2$fig_clus <- subseur2$ord_clus
subseur2$fig_clus[subseur2$ord_clus == 5] <- 4

Idents(subseur2) <- "fig_clus"

marks=FindAllMarkers(subseur2,assay="RNA",slot="data",only.pos = T,pseudocount.use = 1)

write.csv(marks,"Figure1_v3/basal_to_intermed_clus3_vs_clus4_pseudocount1.csv")

clus3_vs_4 <- FindMarkers(subseur2, ident.1=3,ident.2=c(4,5),assay="RNA",slot="data",pseudocount.use = 0.5)

clus3_vs_4$gene <- rownames(clus3_vs_4)
View(clus3_vs_4)
write.csv(clus3_vs_4,"Figure1_v3/basal_to_intermed_clus3_vs_clus4_pseudocount0.5.csv")


#### Notch and Krt8 expression
DefaultAssay(subseur2) <-"RNA"
pdf("Figure1_v3/basal_to_intermed_Notch3_Krt8.pdf",height=5,width=11,useDingbats = F)
FeaturePlot(subseur2, c("Notch3","Krt8"),order=T,max.cutoff = "q99",reduction="umap",pt.size = 2.5,raster=T) & scale_color_viridis()& theme_void()
dev.off()

pdf("Figure1_v3/basal_to_intermed_Scgb1a1_Muc5ac.pdf",height=8,width=12,useDingbats = F)
FeaturePlot(subseur2, c("Notch3","Krt8","Scgb1a1","Muc5ac"),order=T,max.cutoff = "q99",reduction="umap",pt.size = 2.5,raster=T) & scale_color_viridis()& theme_void()
dev.off()

### Hey1 and Mycl correlation
DefaultAssay(subseur2) <- "MAGIC_RNA"

cells=WhichCells(subseur2, idents = c(3,4,5))

group_cols=c("indianred","darkturquoise","indianred")

pdf("Figure1_v3/basal_to_intermed_mycl_hey1_anticorr_MAGIC.pdf",height=4,width=5.5,useDingbats = F)
FeatureScatter(subseur2, "Mycl","Hey1",cells = cells, cols = group_cols,raster = T,pt.size=3)
dev.off()

DefaultAssay(seur) <- "RNA"
FeaturePlot(seur, "Notch3",reduction="monocle",order=T)

pdf("Figure1_v3/Notch1_Notch2_expression.pdf",height=5,width=11,useDingbats = F)
FeaturePlot(seur, c("Notch1","Notch2"),reduction="monocle",pt.size=2,raster = T,order=T,max.cutoff = "q99") & scale_color_viridis() & theme_void()
dev.off()




########### z-score top genes ###########
library(multiClust)
marks =read.csv("Figure1_v3/NS_WtAd3_Ad3_Tdsort_combined_names_markers.csv")

marks_filt = marks %>% filter(p_val_adj < 0.05) %>% filter(avg_log2FC > log2(2))

basal_genes = marks_filt %>% filter(cluster == "Basal")
int_genes = marks_filt %>% filter(cluster == "Int")
sec_genes = marks_filt %>% filter(cluster == "Secretory")
mcc_genes = marks_filt %>% filter(cluster == "MCC")

num_genes = length(int_genes$X)

full_genes = unique(c(basal_genes$gene,int_genes$gene,sec_genes$gene,mcc_genes$gene))

basal_top = basal_genes[order(basal_genes$avg_log2FC,decreasing=T),][1:num_genes,]
sec_top = sec_genes[order(sec_genes$avg_log2FC,decreasing=T),][1:num_genes,]
mcc_top = mcc_genes[order(mcc_genes$avg_log2FC,decreasing=T),][1:num_genes,]

all_genes = unique(c(int_genes$gene,basal_top$gene,sec_top$gene,mcc_top$gene))

gene_exp = seur[["RNA"]]@data[match(all_genes, rownames(seur[["RNA"]]@data)),]

exp_gene_exp =expm1(gene_exp)

##top 1%
cap_gene_exp=exp_gene_exp
for(i in 1:nrow(exp_gene_exp)){
        rng=range(sort(exp_gene_exp[i,], decreasing=TRUE)[1:(length(exp_gene_exp[i,])/100)])
        cap_gene_exp[i,][cap_gene_exp[i,] > rng[1]] <- rng[1]
}
gene_exp_scale = apply(cap_gene_exp, 1, nor.min.max)
gene_exp_sum = rowSums(gene_exp_scale)

gene_exp_scale_nocap = apply(exp_gene_exp, 1, nor.min.max)
gene_exp_nocap_sum = rowSums(gene_exp_scale_nocap)

seur$stick_shift_scale_cap <- gene_exp_sum
seur$stick_shift_scale <- gene_exp_nocap_sum

pdf("Figure1_v3/stick_shift_scale_nocap.pdf",height=8,width=11,useDingbats = F)
FeaturePlot(seur, "stick_shift_scale",reduction="monocle",order=T,pt.size=1) +scale_color_viridis()
dev.off()
VlnPlot(seur, "stick_shift_scale",group.by = "combined_mcc",pt.size=0)

pdf("Figure1_v3/stick_shift_scale_cap.pdf",height=8,width=11,useDingbats = F)
FeaturePlot(seur, "stick_shift_scale_cap",reduction="monocle",order=T,pt.size=1) +scale_color_viridis()
dev.off()

meta=seur@meta.data

library(dittoSeq)
library(multiClust)
library(ggnewscale)

scgb=dittoDimPlot(seur, "Scgb3a1",reduction="monocle",assay = "RNA",data.out = T)
trp63=dittoDimPlot(seur,"Trp63",reduction="monocle",assay="RNA",data.out = T)
foxj1=dittoDimPlot(seur,"Foxj1",reduction="monocle",assay="RNA",data.out = T)
myb=dittoDimPlot(seur,"Myb",reduction="monocle",assay="RNA",data.out = T)


### threshold genes at 1 SD below average for cluster - binarize cells as 0 or 1 - sum and plot
## seur object should be default RNA and idents set to clusters

non_log_exp = expm1(gene_exp)
avg_exp=rowMeans(non_log_exp)
sd_results=apply(non_log_exp,1,sd)

DefaultAssay(seur) <- "RNA"
Idents(seur) <- "combined_mcc"

cluster="Int"
genes=int_genes$gene
find_threshold = function(genes, seur,cluster,mat){
        cells=WhichCells(seur, idents=cluster)
        sub_data=mat[match(genes,rownames(mat)),match(cells,colnames(mat))]
        
        ##average expression
        med_exp= as.data.frame(rowMedians(sub_data))
        colnames(med_exp) <- cluster
        med_exp$gene<-rownames(sub_data)
        med_exp[,"thresh_50"] <- med_exp[,cluster]*0.5
        return(med_exp[,c("thresh_50","gene"),drop=F])
}

Int=find_threshold(int_genes$gene, seur, "Int",exp_gene_exp)
Sec=find_threshold(sec_top$gene, seur, "Secretory",exp_gene_exp)
MCC=find_threshold(mcc_top$gene, seur, "MCC",exp_gene_exp)
Basal=find_threshold(basal_top$gene, seur, "Basal",exp_gene_exp)

df = do.call(rbind, list(Int, Sec, MCC, Basal))

bin_gene_exp <- exp_gene_exp

for(x in df$gene){
        bin_row=bin_gene_exp[match(x, rownames(bin_gene_exp)),]
        bin_row[bin_row <= df$thresh_50[match(x, df$gene)]] = 0
        bin_row[bin_row > df$thresh_50[match(x, df$gene)]] = 1
        bin_gene_exp[match(x, rownames(bin_gene_exp)),] = bin_row
}

colSums(bin_gene_exp)


extra=colSums(bin_gene_exp)
identical(names(extra),rownames(seur@meta.data))
seur$median_50_thresh <- extra

FeaturePlot(seur, "median_50_thresh",order=T,reduction="monocle")+scale_color_viridis()


### use all DE genes
gene_exp_full = seur[["RNA"]]@data[match(full_genes, rownames(seur[["RNA"]]@data)),]

full_gene_exp =expm1(gene_exp_full)

Int_f=find_threshold(int_genes$gene, seur, "Int",full_gene_exp)
Sec_f=find_threshold(sec_genes$gene, seur, "Secretory",full_gene_exp)
MCC_f=find_threshold(mcc_genes$gene, seur, "MCC",full_gene_exp)
Basal_f=find_threshold(basal_genes$gene, seur, "Basal",full_gene_exp)

df_f = do.call(rbind, list(Int_f, Sec_f, MCC_f, Basal_f))

bin_exp_full <- full_gene_exp

for(x in df_f$gene){
        bin_row=bin_exp_full[match(x, rownames(bin_exp_full)),]
        bin_row[bin_row <= df_f$thresh_50[match(x, df_f$gene)]] = 0
        bin_row[bin_row > df_f$thresh_50[match(x, df_f$gene)]] = 1
        bin_exp_full[match(x, rownames(bin_exp_full)),] = bin_row
}

colSums(bin_exp_full)

extra=colSums(bin_exp_full)
identical(names(extra),rownames(seur@meta.data))
seur$median_50_thresh_allDE <- extra

FeaturePlot(seur, "median_50_thresh_allDE",order=T,reduction="monocle")+scale_color_viridis()

VlnPlot(seur,"median_50_thresh_allDE",group.by="combined_mcc",pt.size=0)

####### print out info for GEO
library(Matrix)
meta=seur@meta.data
predictions=grep("prediction.score",colnames(meta))

meta_thin=meta[,-predictions]

meta_final=meta_thin[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","paper_names","combined_mcc","pseudotime","predicted_zaragosi")]

umaps=as.data.frame(Embeddings(seur,reduction = "monocle"))

meta_final$UMAP_1 <- umaps$monocle_1
meta_final$UMAP_2 <- umaps$monocle_2

exp_mat=seur[["RNA"]]@data
exp_mat_counts=seur[["RNA"]]@counts

write.csv(meta_final,file="~/Box Sync/Mycl_paper/ALID03_combined_metadata.csv")
writeMM(exp_mat, file="~/Box Sync/Mycl_paper/ALID03_combined_normdata.mtx")
writeMM(exp_mat_counts, file="~/Box Sync/Mycl_paper/ALID03_combined_counts.mtx")

library("sceasy")

new_meta = seur@meta.data
new_meta=new_meta[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","paper_names","pseudotime","combined_mcc","predicted_zaragosi")]
new_meta$paper_names <- as.character(new_meta$paper_names)
new_meta$paper_names[which(new_meta$paper_names == "Perictye-like")] <- "Tnfrsf12a+"
new_meta$combined_mcc <- as.character(new_meta$combined_mcc)
new_meta$combined_mcc[which(new_meta$combined_mcc =="Perictye-like")]<- "Tnfrsf12a+"
new_meta$combined_mcc[which(new_meta$combined_mcc =="Int")]<- "Intermediate"
new_meta$combined_mcc[which(new_meta$combined_mcc =="G2M")]<- "Cycling (G2M)"
new_meta$combined_mcc[which(new_meta$combined_mcc =="S")]<- "Cycling (S)"
new_meta$combined_mcc[which(new_meta$combined_mcc =="MCC")]<- "Multiciliated"

new_meta$predicted_zaragosi[which(new_meta$predicted_zaragosi == "Precursor")]<-"Undefined rare"

new_meta$orig.ident[which(new_meta$orig.ident == "NS")]<-"unsorted (Foxj1CreERT RosatdTomato)"
new_meta$orig.ident[which(new_meta$orig.ident == "Tdsort")]<-"tdTomato-enriched (Foxj1CreERT RosatdTomato)"
new_meta$orig.ident[which(new_meta$orig.ident == "Ad3")]<-"wild type replicate 1"

new_meta$orig.ident[which(new_meta$orig.ident == "WtAd3")]<-"wild type replicate 2"

## add ontology terms
new_meta$organism_ontology_term_id <-"NCBITaxon:10090"
new_meta$tissue_ontology_term_id <- "UBERON:0001901"
new_meta$assay_ontology_term_id <- "EFO:0008913" 
new_meta$disease_ontology_term_id <- "PATO:0000461"
new_meta$cell_type_ontology_term_id <- "CL:0000307"
new_meta$self_reported_ethnicity_ontology_term_id <- "na"
new_meta$development_stage_ontology_term_id <- "MmusDv:0000110"
new_meta$sex_ontology_term_id <- "PATO:0000384 and PATO:0000383"
new_meta$donor_id <- "na"
new_meta$suspension_type <- "cell"

seur_temp=seur
seur_temp[["RNA"]]@counts <- round(seur[["RNA"]]@counts)

seur_temp@meta.data = new_meta
sceasy::convertFormat(seur_temp, assay='RNA', from="seurat", to="anndata", main_layer='data', transfer_layers='counts', drop_single_values=FALSE, outFile='ALID03_dataset_v2.h5ad')


### SPRING output


