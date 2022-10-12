library(slingshot)
library(condiments)
library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(tidyr)
library(SingleCellExperiment)
library(bioc2020trajectories)
library(tradeSeq)
library(viridis)
library(ggrastr)
library(VennDiagram)

cds=readRDS("v1/monocle/early_int_v1_monocle.rds")
seur=readRDS("v1/early_int_v1.rds")

dir.create("slingshot")
dir.create("slingshot/v2")

small_cds=choose_cells(cds) ## choose only the main body of cells
plot_cells(small_cds,label_cell_groups = T,graph_label_size = 0,labels_per_group = 5,color_cells_by = "ord_clus")

## simplify cluster labels for slingshot

smaller_cds=small_cds[,rownames(colData(small_cds))[colData(small_cds)$ord_clus %in% c(1:6)]]
colData(smaller_cds)$sling_clus <- colData(smaller_cds)$ord_clus
colData(smaller_cds)$sling_clus[colData(smaller_cds)$ord_clus %in% (3:5)] <- 3
colData(smaller_cds)$sling_clus[colData(smaller_cds)$ord_clus %in% (6)] <- 4

plot_cells(smaller_cds,label_cell_groups = T,graph_label_size = 0,labels_per_group = 5,color_cells_by = "sling_clus")
saveRDS(smaller_cds, file="slingshot/v2/smaller_cds_sling.rds")
smaller_cds=readRDS("slingshot/v2/smaller_cds_sling.rds")

sce=slingshot(smaller_cds, clusterLabels = colData(smaller_cds)$sling_clus,start.clus=1, end.clus=c(3,4))

library(grDevices)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')

saveRDS(sce, file="slingshot/v2/early_int_v2_slingshot.rds")
sce=readRDS(file="slingshot/v2/early_int_v2_slingshot.rds")

crv = SlingshotDataSet(sce)
saveRDS(crv, file="slingshot/v2/early_int_v2_crv_sling.rds")


#####################  plot pseudotime values for two lineages ################ ################ 

df = as.data.frame(reducedDims(sce)$UMAP)
df=cbind(df,slingPseudotime(sce))
df=cbind(df,slingCurveWeights(sce))
colnames(df) <- c(colnames(df)[1:4],paste0(colnames(df[5:6]),"_cellweights"))
df=cbind(df, orig.ident=sce$orig.ident)
df$orig.ident <- factor(df$orig.ident, levels=c("NT1_KO1_soupx","Gmnc_soupx","Mcidas_soupx"))
df=cbind(df, ord_clus = sce$ord_clus)
df$sling_clus <- sce$sling_clus


curv_df=slingCurves(sce, as.df=T)
mcc_curv=filter(curv_df, Lineage == 2)
sec_curv=filter(curv_df, Lineage == 1)

mcc_color <- colorRampPalette(c("lightskyblue2", "royalblue4"))
fig1_col = c("#F7A12D","#009245",rep(mcc_color(5)[3],3),"indianred","tomato3","darkgoldenrod","#998675","#996699","palegreen1")


pdf("slingshot/v2/early_int_v2_slingshot_bothpath.pdf",height=7.5,width=8,useDingbats = F)
ggplot(df, aes(x=UMAP_1,y=UMAP_2, color=ord_clus))+
        ggrastr::rasterize(geom_point())+
        scale_color_manual(values=fig1_col)+
        geom_path(data=mcc_curv %>% arrange(Order), aes(group=Lineage),col="black",size=2)+
        geom_path(data=sec_curv %>% arrange(Order), aes(group=Lineage),col="black",size=2,linetype="dashed")+
        theme_void()+
        theme(legend.position = "none" )

dev.off()

pdf("slingshot/v2/early_int_v2_slingshot_mccpath.pdf",height=7.5,width=8,useDingbats = F)
ggplot(df, aes(x=UMAP_1,y=UMAP_2, color=Lineage2))+
        ggrastr::rasterize(geom_point())+
        geom_path(data=mcc_curv %>% arrange(Order), aes(group=Lineage),col="black",size=2)+
        scale_color_gradient(low="gold1",high="darkviolet",na.value="grey")+
        theme_void()+
        labs(color = "Pseudotime")+
        theme(legend.position = "top" )
        
dev.off()

pdf("slingshot/v2/early_int_v2_slingshot_secpath.pdf",height=7.5,width=8,useDingbats = F)
ggplot(df, aes(x=UMAP_1,y=UMAP_2, color=(Lineage1)))+
        ggrastr::rasterize(geom_point())+
        geom_path(data=sec_curv %>% arrange(Order), aes(group=Lineage),col="black",size=2)+
        scale_color_gradient(low="gold1",high="darkviolet",na.value="grey")+
        theme_void()+
        labs(color = "Pseudotime")+
        theme(legend.position = "top" )
dev.off()


##################  imbalance score ################ ################ 

scores <- bioc2020trajectories::imbalance_score(
        rd = reducedDims(sce)$UMAP, 
        cl = df$orig.ident,
        k = 20, smooth = 40)

df$scaled_scores = scores$scaled_scores

pdf("slingshot/v2/early_int_v2_slingshot_imbalance.pdf",height=8,width=11.5,useDingbats = F)
ggplot(df, aes(x=UMAP_1,y=UMAP_2, color=scaled_scores))+
        ggrastr::rasterise(geom_point(size=1.5))+
        scale_color_viridis(option= "D")+
        theme_void()
dev.off()

############### Topology test ################ ################ 

top_res <- topologyTest(sds = sce, conditions = df$orig.ident)

################ ################ Pseudotime distribution ################ ################ 
mcc_color <- colorRampPalette(c("lightskyblue2", "royalblue4"))
fig1_col = c("#998675","#996699","#F7A12D","#009245","tomato3","palevioletred2",mcc_color(5),"palegreen1","darkgoldenrod")
cds_col = c("#F7A12D","#009245",mcc_color(5),"palegreen1","grey17","grey17")
col_val =c("1" = "#F7A12D", "2"="#009245", "3" = mcc_color(5)[3],"4"=mcc_color(5)[3],"5"=mcc_color(5)[3],"6"="indianred","7"="palegreen1")


pdf("slingshot/v2/early_int_v2_slingshot_pseudo_mcc_distr_raster600.pdf",height=6,width=8,useDingbats = F)
ggplot(df[df$ord_clus %in% c(1:7),], aes(x=orig.ident,y=Lineage2))+
        geom_violin(na.rm = T,size=1)+
        ggrastr::rasterise(geom_jitter(aes(color=ord_clus),na.rm = T,width=0.1,cex=0.5,alpha=1),raster.dpi=600)+
        scale_color_manual(values=col_val)+
        theme_classic()+
        xlab(NULL)+
        ylab("MCC Path Pseudotime")+
        theme(legend.position = "none",axis.text = element_text(size=12), axis.title.y = element_text(size=14))+
        scale_x_discrete(labels=c("NT KO", "Gmnc KO", "Mcidas KO"))
dev.off()

pdf("slingshot/v2/early_int_v2_slingshot_pseudo_sec_distr.pdf",height=6,width=8,useDingbats = F)
ggplot(df[df$ord_clus %in% c(1:7),], aes(x=orig.ident,y=Lineage1))+
        geom_violin(na.rm = T,size=1)+
        ggrastr::rasterise(geom_jitter(aes(color=ord_clus),na.rm = T,width=0.1,cex=0.5,alpha=0.7))+
        scale_color_manual(values=col_val)+
        theme_classic()+
        xlab(NULL)+
        ylab("MCC Path Pseudotime")+
        theme(axis.text = element_text(size=12), axis.title.y = element_text(size=14))+
        scale_x_discrete(labels=c("NT KO", "Gmnc KO", "Mcidas KO"))
dev.off()


sce_nt1_gmnc = sce[,sce$orig.ident == "NT1_KO1_soupx" | sce$orig.ident == "Gmnc_soupx"]
sce_nt1_mcidas=sce[,sce$orig.ident == "NT1_KO1_soupx" | sce$orig.ident == "Mcidas_soupx"]

prog_res_gmnc <- progressionTest(slingPseudotime(sce_nt1_gmnc), cellWeights=slingCurveWeights(sce_nt1_gmnc),conditions = sce_nt1_gmnc$orig.ident, global = TRUE, lineages = TRUE)

prog_res_mcidas <- progressionTest(slingPseudotime(sce_nt1_mcidas), cellWeights=slingCurveWeights(sce_nt1_mcidas),conditions = sce_nt1_mcidas$orig.ident, global = TRUE, lineages = TRUE)


################# Cluster 2 DE genes #################
Idents(seur) <- "ord_clus"

gmnc_marks = FindMarkers(seur, assay = "RNA",ident.1 = "Gmnc_soupx",ident.2 = "NT1_KO1_soupx",group.by = "orig.ident",subset.ident = "2")
mcidas_marks =  FindMarkers(seur, assay = "RNA",ident.1 = "Mcidas_soupx",ident.2 = "NT1_KO1_soupx",group.by = "orig.ident",subset.ident = "2")

write.csv(gmnc_marks, file="early_int_cluster2_Gmnc_vs_NT1_markers.csv")
write.csv(mcidas_marks, file="early_int_cluster2_Mcidas_vs_NT1_markers.csv")

gmnc_marks_bas = FindMarkers(seur, assay = "RNA",ident.1 = "Gmnc_soupx",ident.2 = "NT1_KO1_soupx",group.by = "orig.ident",subset.ident = "1")
mcidas_marks_bas =  FindMarkers(seur, assay = "RNA",ident.1 = "Mcidas_soupx",ident.2 = "NT1_KO1_soupx",group.by = "orig.ident",subset.ident = "1")

write.csv(gmnc_marks_bas, file="early_int_cluster1_Gmnc_vs_NT1_markers.csv")

### genes with > log2(1.5) FC and p_adj_val < 0.05

gmnc_marks_filt = gmnc_marks %>% filter(p_val_adj < 0.05) %>% filter(abs(avg_log2FC) > log2(1.5))
mcidas_marks_filt = mcidas_marks %>% filter(p_val_adj < 0.05) %>% filter(abs(avg_log2FC) > log2(1.5))

write.csv(gmnc_marks_bas, file="early_int_cluster1_Mcidas_vs_NT1_markers.csv")

####### DE of signature genes from Mycl and Hey1 correlation

corr_genes = read.csv("/Volumes/Reiterlab_3/NS_WtAd3_Ad3_Tdsort_agg/basal_to_intermed_v2/basal_to_intermed_v2_clean_Mycl_vs_Hey1_correlated_genes_thresh0.15_heatmap.csv")

hey1_genes = corr_genes[51:65,]

Idents(seur) <- "ord_clus"

gmnc_corr_DE = FindMarkers(seur, assay = "RNA",ident.1 = "Gmnc_soupx",ident.2 = "NT1_KO1_soupx",group.by = "orig.ident",subset.ident = "2",features=corr_genes$X,min.pct = 0,logfc.threshold = 0)

mcidas_corr_DE = FindMarkers(seur, assay = "RNA",ident.1 = "Mcidas_soupx",ident.2 = "NT1_KO1_soupx",group.by = "orig.ident",subset.ident = "2",features=corr_genes$X,min.pct = 0,logfc.threshold = 0)

mcidas_corr_DE$gene <- rownames(mcidas_corr_DE)
mcidas_corr_DE$corr_genes <- "Mycl"
mcidas_corr_DE$corr_genes[mcidas_corr_DE$gene %in% hey1_genes$X] <- "Hey1"

write.csv(mcidas_corr_DE, file="Mcidas_vs_NT_hey1_Mycl_correlated_genes_DE.csv")
m
#mcidas_corr_DE= read.csv("Mcidas_vs_NT_hey1_Mycl_correlated_genes_DE.csv",row.names = 1)

gmnc_corr_DE_order = gmnc_corr_DE[order(gmnc_corr_DE$avg_log2FC),]
gmnc_corr_DE_order$FC <- 2^gmnc_corr_DE_order$avg_log2FC
gmnc_corr_DE_order$gene = rownames(gmnc_corr_DE_order)
gmnc_corr_DE_order$corr_genes <- "Mycl"
gmnc_corr_DE_order$corr_genes[gmnc_corr_DE_order$gene %in% hey1_genes$X] <- "Hey1"
gmnc_corr_DE_order
gmnc_corr_filt = gmnc_corr_DE_order %>% filter(p_val_adj < 0.05)

write.csv(gmnc_corr_DE_order, file="Gmnc_vs_NT_hey1_Mycl_correlated_genes_DE.csv")
#gmnc_corr_DE_order=read.csv("Gmnc_vs_NT_hey1_Mycl_correlated_genes_DE.csv",row.names=1)

##### volcano plot
library('EnhancedVolcano')

mcidas_corr_mycl = mcidas_corr_DE %>% filter(corr_genes == "Mycl")

pdf("Figure/mcidas_volcano.pdf",height=6,width=8,useDingbats = F)
EnhancedVolcano(mcidas_corr_mycl,
                lab = rownames(mcidas_corr_mycl),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                FCcutoff = log2(1.5),
                col=c('black', 'black', 'black', 'red3'),
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                drawConnectors=TRUE,
                pointSize = 4,
                colAlpha = 1,
                xlim = c(-1.5,1.5))
dev.off()

gmnc_corr_mycl = gmnc_corr_DE_order %>% filter(corr_genes == "Mycl")

pdf("Figure/gmnc_volcano.pdf",height=8,width=8,useDingbats = F)
EnhancedVolcano(gmnc_corr_mycl,
                lab = rownames(gmnc_corr_mycl),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                FCcutoff = log2(1.5),
                selectLab = c('Mycl','Gmnc','Ccno','Myb'),
                col=c('black', 'black', 'black', 'red3'),
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                drawConnectors=TRUE,
                pointSize = 6,
                colAlpha = 1,
                xlim = c(-1.5,1.5))
dev.off()


