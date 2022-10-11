 ## NS_WtAd3_Ad3_Tdsort
## soupX corrected, seurat integrated

library(Seurat)
library(dplyr)
library(monocle3)

basename = "NS_WtAd3_Ad3_Tdsort_agg"
vers="v1"

seur=readRDS("v1/NS_WtAd3_Ad3_Tdsort_agg_v1.rds")

seur <- RunUMAP(seur, reduction="pca",dims=1:30)

seur = FindNeighbors(seur, dims=1:30)
seur = FindClusters(seur, resolution = 0.3)

DefaultAssay(seur) <- "RNA"
marks =FindAllMarkers(seur,only.pos = T,assay = "RNA")
marks_sort = group_by(marks, cluster) %>% filter(p_val_adj < 0.05)

write.csv(marks_sort, file=paste0(vers, "/",basename,"_",vers,"_res0.3_markers.csv"))

##### monocle3


#Expression Matrix - select only cells of interest from seurat_ob@data
exprs<- GetAssayData(seur, slot="counts",assay = "RNA")

## phenodata
pheno.data = seur@meta.data

## feature data
genes <- data.frame(gene_short_name = rownames(exprs))
rownames(genes) <- rownames(exprs)


cds <- new_cell_data_set(exprs,
                         cell_metadata = pheno.data,
                         gene_metadata = genes)


set.seed(123)
cds=preprocess_cds(cds, num_dim=30)
plot_pc_variance_explained(cds) ## check number of PCs included
cds=align_cds(cds, preprocess_method="PCA",alignment_group="orig.ident")
cds <- reduce_dimension(cds, umap.fast_sgd = FALSE,preprocess_method = "Aligned")
plot_cells(cds, label_groups_by_cluster=FALSE)

cds <- cluster_cells(cds,cluster_method="leiden",random_seed=123)

plot_cells(cds,color_cells_by = "partition")
plot_cells(cds,color_cells_by = "cluster")


cds <- learn_graph(cds,use_partition = T,learn_graph_control = list(minimal_branch_len=25))

plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE,group_label_size = 18)

cds=order_cells(cds)

#### reorder clusters
coldat = colData(cds)
coldat$ord_clus <- clusters(cds)
coldat$mon_clus <- clusters(cds)
coldat$ord_clus[coldat$mon_clus == 9] <- 1
coldat$ord_clus[coldat$mon_clus == 10] <- 2
coldat$ord_clus[coldat$mon_clus == 1] <- 3
coldat$ord_clus[coldat$mon_clus == 2] <- 4
coldat$ord_clus[coldat$mon_clus == 4] <- 5
coldat$ord_clus[coldat$mon_clus == 3] <- 6
coldat$ord_clus[coldat$mon_clus == 8] <- 7
coldat$ord_clus[coldat$mon_clus == 6] <- 8
coldat$ord_clus[coldat$mon_clus == 5] <- 9
coldat$ord_clus[coldat$mon_clus == 7] <- 10
coldat$ord_clus[coldat$mon_clus == 11] <- 11
coldat$ord_clus[coldat$mon_clus == 12] <- 12
coldat$ord_clus[coldat$mon_clus == 13] <- 13

##
saveRDS(cds, file=paste0(vers, "/",basename,"_",vers,"_aligned.rds"))
cds=readRDS(file=paste0(vers, "/",basename,"_",vers,"_aligned.rds"))


##### add monocle dimensions and clusters back into Seurat object ########

seur=readRDS("NS_WtAd3_Ad3_Tdsort_v1_magic.rds")

monocle_clus=colData(cds)[,"mon_clus"]
seur$monocle_clus <- NA
seur$monocle_clus[match(names(monocle_clus),rownames(seur@meta.data))] <- monocle_clus

mon_umap=reducedDim(cds,type = "UMAP")
mon_umap_ord = mon_umap[match(rownames(seur@meta.data),rownames(mon_umap)),]
colnames(mon_umap_ord) <- c("UMAP 1","UMAP 2")

seur[["monocle"]] <- CreateDimReducObject(embeddings = mon_umap_ord, key = "monocle_", assay = "RNA")

DimPlot(seur, group.by = "monocle_clus",reduction="monocle")

saveRDS(seur, file="NS_WtAd3_Ad3_Tdsort_v1_magic.rds")

