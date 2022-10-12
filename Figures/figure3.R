## early_int figure
## many of the plots in Figure 3 were printed from early_int_slingshot.R or early_int_proportion_testing.R scripts.

library(Seurat)
library(monocle3)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(pheatmap)

vers="v1"
basename="early_int"
seur=readRDS(paste0(vers,"/",basename,"_",vers,".rds"))
DimPlot(seur)


dir.create("Figure")
mcc_color <- colorRampPalette(c("lightskyblue2", "royalblue4"))
fig1_col = c("#F7A12D","#009245",mcc_color(5)[3],"indianred","palevioletred2","darkgoldenrod","#998675","#996699")

### combine MCCs

seur$combined_mccs = seur$ord_clus
seur$combined_mccs[seur$ord_clus == 4] <- 3
seur$combined_mccs[seur$ord_clus == 5] <- 3
seur$combined_mccs[seur$ord_clus ==10 ] <- 11 ### cluster 10 doesn't have any DE gene vs. cluster 11.
DimPlot(seur, group.by="combined_mccs",label=T)

Idents(seur) <- "combined_mccs"
marks=FindAllMarkers(seur, assay="RNA",only.pos = T, test.use = "MAST",latent.vars = "nCount_RNA")

write.csv(marks, file="early_int_combined_mccs_markers.csv")

DimPlot(seur, group.by="combined_mccs")
umap_plot=DimPlot(seur, cols = fig1_col,group.by="combined_mccs",pt.size=2,raster=T)+ theme_void() + theme(legend.position = "none") +ggtitle("")

umap2 = DimPlot(seur, cols = fig1_col,group.by="combined_mccs",pt.size=2,raster=T,split.by="orig.ident")+ theme_void() +theme(legend.position="none") + ggtitle("")

pdf("Figure/early_int_v1_umap_big_split.pdf",height=4.5,width=12,useDingbats = F)
ggarrange(umap_plot, umap2, nrow=1,widths = c(1,3))
dev.off()

### umap 


umap_plot=DimPlot(seur, cols = fig1_col,group.by="ord_clus",pt.size=2,raster=T)+ theme_void() + ggtitle("") & NoLegend()

DefaultAssay(seur) <- "RNA"
feat_plot=FeaturePlot(seur, features=c("Trp63","Scgb3a1","Myb","Foxj1"),ncol = 2,order=T,max.cutoff = "q99",raster = T,pt.size=2) & theme_void() & scale_color_viridis()

pdf("Figure/early_int_v1_umap_feature_nolegend.pdf",height=4.5,width=16,useDingbats = F)
ggarrange(umap_plot, feat_plot, nrow=1,widths = c(1,4))
dev.off()

pdf("Figure/early_int_v1_feature_nolegend.pdf",height=8,width=6,useDingbats = F)
feat_plot
dev.off()



### split UMAP


pdf("Figure/early_int_v1_split_UMAP.pdf",height=4.5,width=11,useDingbats = F)
DimPlot(seur, group.by="ord_clus",split.by = "orig.ident",cols = fig1_col,label=F,pt.size = 2,raster = T)+
        theme_void()+
        theme(legend.position = "none",strip.text = element_text(size=14))+
        ggtitle("")
dev.off()

#### just main body 

smaller_cds=readRDS("slingshot/v2/smaller_cds_sling.rds")
cells_plot = rownames(colData(smaller_cds))


umap_small=DimPlot(seur, group.by="ord_clus",cols=fig1_col,pt.size = 2,raster=T,cells = cells_plot) + 
        theme_void() + 
        theme(legend.position = "none",strip.text = element_text(size=14))+
        ggtitle("")

DefaultAssay(seur) <- "RNA"
feat_plot=FeaturePlot(seur, features=c("Trp63","Foxj1","Scgb3a1"),ncol =3,order=T,max.cutoff = "q99",raster = T,pt.size=2,cells=cells_plot) & theme_void() & scale_color_viridis() & NoLegend()

pdf("Figure/early_int_v1_mainbody_umap_feature_nolegend.pdf",height=4.5,width=16,useDingbats = F)
ggarrange(umap_small, feat_plot, nrow=1,widths = c(1.5,4))
dev.off()

feat_plot2=FeaturePlot(seur, features=c("Mcidas","Foxj1","Rfx2"),ncol =3,order=T,max.cutoff = "q99",raster = T,pt.size=2,cells=cells_plot) & theme_void() & scale_color_viridis() & NoLegend()

pdf("Figure/early_int_v1_mainbody_umap_MCCs_nolegend.pdf",height=4.5,width=16,useDingbats = F)
ggarrange(umap_small, feat_plot2, nrow=1,widths = c(1.5,4))
dev.off()

small_seur = subset(seur, cells=cells_plot)
small_seur$orig.ident <- factor(small_seur$orig.ident, levels=c("NT1_KO1_soupx", "Gmnc_soupx","Mcidas_soupx"))
split_umap <-DimPlot(small_seur, group.by="ord_clus",cols=fig1_col,pt.size = 2,raster=T, split.by="orig.ident") + 
        theme_void() + 
        theme(legend.position = "none",strip.text = element_text(size=14))+
        ggtitle("")

pdf("Figure/early_int_v1_mainbody_umap_split_nolegend.pdf",height=4.5,width=16,useDingbats = F)
ggarrange(umap_small, split_umap, nrow=1,widths = c(1.5,4))
dev.off()

####################### LABEL TRANSFER ##################

##zaragosi human dataset predictions frmo wynton

pred = read.csv("../label_transfer_datasets/zaragosi_human/early_int_v1_predictions.csv",row.names = 1)

pred_meta =pred[,c("predicted.id","prediction.score.max"),drop=F]

seur=AddMetaData(seur, pred)

library(dittoSeq)
cols=dittoColors(reps=1)
cols_ord = c(cols[1],cols[4],cols[12],cols[2],cols[3],cols[9],cols[10],cols[8],cols[7],cols[5],cols[6],cols[11])

pdf("early_int_predictions_zaragosi.pdf",height=8,width=11,useDingbats = F)
DimPlot(seur, group.by="predicted.id", label = T, cols = cols_ord,pt.size=1,raster = T, split.by = "orig.ident") + theme_void()
dev.off()


##### print out dataset

library(Matrix)

exp_counts = seur[["RNA"]]@counts
exp_norm=seur[["RNA"]]@data

writeMM(exp_counts, file="KO_combined_counts.mtx")
writeMM(exp_norm, file="KO_combined_normdata.mtx")

meta=seur@meta.data
meta_final = meta[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","ord_clus")]

embed=as.data.frame(Embeddings(seur, reduction="umap"))

meta_final$UMAP_1 <- embed$UMAP_1
meta_final$UMAP_2 <- embed$UMAP_2

sce=readRDS(file="slingshot/v2/early_int_v2_slingshot.rds")
sling=as.data.frame(slingPseudotime(sce))
curve = as.data.frame(slingCurveWeights(sce))

meta_final$MCC_sling_pseudotime <- NA
meta_final$MCC_sling_pseudotime[match(rownames(sling),rownames(meta_final))] <-sling$Lineage2 
meta_final$Sec_sling_pseudotime <- NA
meta_final$Sec_sling_pseudotime[match(rownames(sling),rownames(meta_final))] <-sling$Lineage1 

meta_final$MCC_sling_curveweights <- NA
meta_final$MCC_sling_curveweights[match(rownames(curve),rownames(meta_final))]<- curve$Lineage2
meta_final$Sec_sling_curveweights <- NA
meta_final$Sec_sling_curveweights[match(rownames(curve),rownames(meta_final))]<- curve$Lineage1


write.csv(meta_final, file="KO_combined_meta.csv")

meta_final = read.csv("KO_combined_meta.csv",row.names=1)

library(Seurat)
library("sceasy")

new_meta = seur@meta.data

new_meta=new_meta[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt")]
new_meta$ord_clus <- meta_final$ord_clus
new_meta$MCC_sling_pseudotime <- meta_final$MCC_sling_pseudotime
new_meta$MCC_sling_curveweights <-meta_final$MCC_sling_curveweights
new_meta$Sec_sling_pseudotime <- meta_final$Sec_sling_pseudotime
new_meta$Sec_sling_curveweights <- meta_final$Sec_sling_curveweights

new_meta$ord_clus[new_meta$ord_clus == 10]
new_meta$paper_names <- new_meta$ord_clus
new_meta$paper_names[new_meta$ord_clus == 1] <- "Basal"
new_meta$paper_names[new_meta$ord_clus == 2] <- "Intermediate"
new_meta$paper_names[new_meta$ord_clus == 3] <- "Multiciliated"
new_meta$paper_names[new_meta$ord_clus == 4] <- "Multiciliated"
new_meta$paper_names[new_meta$ord_clus == 5] <- "Multiciliated"
new_meta$paper_names[new_meta$ord_clus == 6] <- "Secretory"
new_meta$paper_names[new_meta$ord_clus == 7] <- "Isg15+"
new_meta$paper_names[new_meta$ord_clus == 8] <- "Neuroendocrine"
new_meta$paper_names[new_meta$ord_clus == 9] <- "Tnfrsf12a+"
new_meta$paper_names[new_meta$ord_clus == 11] <- "Cycling"

new_meta$orig.ident[which(new_meta$orig.ident == "NT1_KO1_soupx")]<-"nontargeted control"
new_meta$orig.ident[which(new_meta$orig.ident == "Gmnc_soupx")]<-"Gmnc KO"
new_meta$orig.ident[which(new_meta$orig.ident == "Mcidas_soupx")]<-"Mcidas KO"

## add ontology terms
new_meta$organism_ontology_term_id <-"NCBITaxon:10090"
new_meta$tissue_ontology_term_id <- "UBERON:0001901"
new_meta$assay_ontology_term_id <- "EFO:0008913" 
new_meta$disease_ontology_term_id <- "MONDO:0016575"
new_meta$cell_type_ontology_term_id <- "CL:0000307"
new_meta$self_reported_ethnicity_ontology_term_id <- "na"
new_meta$development_stage_ontology_term_id <- "MmusDv:0000110"
new_meta$sex_ontology_term_id <- "PATO:0000384 and PATO:0000383"
new_meta$donor_id <- "na"
new_meta$suspension_type <- "cell"

seur_temp=seur
seur_temp[["RNA"]]@counts <- round(seur[["RNA"]]@counts)

seur_temp@meta.data <- new_meta


sceasy::convertFormat(seur_temp, assay='RNA', from="seurat", to="anndata", main_layer='data', transfer_layers='counts', drop_single_values=FALSE, outFile='KO_dataset_v2.h5ad')
