### Figure 2 plots

library(monocle3)
library(Seurat)
library(ggplot2)
library(multiClust)
library(dplyr)
library(reshape2)
library(viridis)
setwd("/Volumes/Reiterlab_3/NS_WtAd3_Ad3_Tdsort_agg/")

seur=readRDS("NS_WtAd3_Ad3_Tdsort_v1_magic.rds")
cds=readRDS(file=paste0("v1/NS_WtAd3_Ad3_Tdsort_agg_v1_aligned.rds"))

seur$paper_names <- factor(seur$paper_names, levels=c("Basal","Int","MCC 1","MCC 2","MCC 3","MCC 4","MCC 5","Secretory","G2M","S","Perictye-like","Neuroendocrine"))

DimPlot(seur, reduction="monocle",group.by="paper_names",label=T)

identical(rownames(colData(cds)),rownames(seur@meta.data))

colData(cds) <- cbind(colData(cds),paper_names=seur$paper_names)

plot_cells(cds, color_cells_by = "paper_names")

plot_cells(cds)

##### Add pseudotime values from basal to MCC
cds_unbias=readRDS(file=paste0("v1/basal_to_endmcc_v1_cds_unbias.rds"))
pseudos=pseudotime(cds_unbias)
cds_coldata=colData(cds)
cds_coldata$pseudo_basal_endmcc <-NA
cds_coldata$pseudo_basal_endmcc[match(names(pseudos),rownames(cds_coldata))] <- pseudos

colData(cds) <- cds_coldata
plot_cells(cds, color_cells_by = "pseudo_basal_endmcc")
plot_cells(cds, color_cells_by = "paper_names",)

pdf("Fig_2/pseudotime_monocle.pdf",height=8,width=11,useDingbats = F)
plot_cells(cds,
           color_cells_by = "pseudo_basal_endmcc",
           show_trajectory_graph = T,
           label_cell_groups = T,
           cell_size = 2,
           cell_stroke = 0.1,
           alpha = 1,
           labels_per_group = F,
           rasterize = T, 
           label_branch_points = F, 
           label_roots = F, 
           label_leaves = F, 
           trajectory_graph_color = "black",
           trajectory_graph_segment_size = 2)+
        scale_color_viridis_c(option="inferno",direction = -1)+
        theme_void()
dev.off()

seur$pseudo_basal_endmcc <- cds$pseudo_basal_endmcc

pdf("Fig_2/pseudotime_umap.pdf",height=8,width=10,useDingbats = F)
FeaturePlot(seur, "pseudo_basal_endmcc",reduction="monocle",pt.size=1.5,raster = T)+
        scale_color_viridis_c(option="inferno",direction = -1)+
        theme_void()
dev.off()

pap=colData(cds)[,"paper_names",drop=F]
pap_unbias = pap[rownames(pap) %in% rownames(colData(cds_unbias)),,drop=F]
identical(rownames(pap_unbias),rownames(colData(cds_unbias)))
cds_unbias$paper_names <- pap_unbias$paper_names

#### pseudo-bin plots

mcc_color <- colorRampPalette(c("lightskyblue2", "royalblue4"))

fig1_col = c("#F7A12D","#009245",rep(mcc_color(5)[3],5),"indianred","#996699","palegreen1","#998675","darkgoldenrod")

pseudo_vals=pseudotime(cds_unbias)
sub=subset(seur, cells=names(pseudo_vals))
identical(rownames(sub@meta.data),names(pseudo_vals))
sub$pseudo <- pseudo_vals
sub$Pseudo_bin <- as.numeric(cut_interval(sub$pseudo,10))

avg_exp=AverageExpression(sub, assay="RNA",features=c("Gmnc","Mcidas","Mycl"),group.by = "Pseudo_bin")
avg_norm = t(apply(avg_exp$RNA, 1, nor.min.max))
colnames(avg_norm) <- paste0("Bin",colnames(avg_norm))
plot_val = melt(avg_norm)


pdf("Fig_2/pseudotime_bin_basal_to_mcc.pdf",height=8,width=11,useDingbats = F)
ggplot(plot_val, aes(x=Var2,y=value,col=Var1,group=Var1))+
        geom_line(size=2,stat="identity")+
        theme_classic()+
        xlab("Pseudotime")+
        ylab("Gene Expression")+
        scale_color_manual(values=c("Gold1","Magenta","Cyan3"))
dev.off()

### Pigr expression
DefaultAssay(seur) <- "RNA"
pdf("Fig_2/pigr_foxj1_mycl_pigr_expression.pdf",height=3,width=11,useDingbats = F)
FeaturePlot(seur, c("Pigr","Foxj1","Mycl"),ncol=3,reduction="monocle",pt.size=2, raster=T,order=T)&scale_color_viridis() & theme_void()
dev.off()

pdf("Fig_2/gmnc_mcidas_mycl_expression.pdf",height=3,width=11,useDingbats = F)
FeaturePlot(seur, c("Gmnc","Mcidas","Mycl"),ncol=3,reduction="monocle",pt.size=2, raster=T,order=T)&scale_color_viridis() & theme_void()
dev.off()
