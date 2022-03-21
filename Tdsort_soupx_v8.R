## Tdsort with updated scrublet and soupx correction

library(Seurat)
library(Matrix)
library(SoupX)


out_folder ="/Volumes/LB_ReiterLab/Reiter_Seq/10X_042519/CELLRANGER/crTdTom/outs/"

sc=load10X(out_folder)
sc = autoEstCont(sc)
out = adjustCounts(sc)

dir.create("Tdsort_soupx")
writeMM(out, file="Tdsort_soupx/matrix.mtx.gz")

## read in Soupx corrected data

data_dir ="Tdsort_soupx/"
basename="Tdsort"
folder="v8"


data = Read10X(data.dir=data_dir)

seur <- CreateSeuratObject(counts = data, project = basename, min.cells = 3, min.features = 200)
seur[["percent.mt"]] <- PercentageFeatureSet(object = seur, pattern = "^mt-")
VlnPlot(seur,features = c("nFeature_RNA","nCount_RNA","percent.mt"),group.by="orig.ident")

seur <- subset(seur, subset = nFeature_RNA > 4000  & percent.mt < 10)
VlnPlot(seur, features =c("nFeature_RNA","percent.mt","nCount_RNA"))

ncount_mean = mean(seur$nCount_RNA)
std_count = sd(seur$nCount_RNA)
ncount_hi = ncount_mean + 2*std_count

seur <- subset(seur, subset = nCount_RNA < ncount_hi)

seur <- NormalizeData(seur, normalization.method = "LogNormalize", scale.factor = 10000)
seur <- FindVariableFeatures(seur, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(seur)
seur <- ScaleData(seur, features = all.genes)

seur <- RunPCA(seur, features = VariableFeatures(object = seur))

png(paste0(folder,"/",basename,"_",folder,"_elbow.png"),height = 800,width=1100)
ElbowPlot(seur)
dev.off()

seur <- FindNeighbors(seur, dims = 1:15)
seur <- FindClusters(seur, resolution = 0.3)

seur <- RunUMAP(seur, dims = 1:15)
DimPlot(seur, reduction = "umap")


dubs = read.csv("v8/Tdsort_scrublet_3.5perc_thresh0.14_output_table_rerun.csv")
identical(dubs$cell_barcodes,rownames(seur@meta.data))
seur$scrublet <- dubs$predicted_doublet
table(dubs$predicted_doublet)

DimPlot(seur, group.by="scrublet")
table(seur$scrublet,seur$RNA_snn_res.0.3)

##remove cluster 5 and any remaining doublets

Idents(seur) <- "scrublet"

seur_sub <- subset(seur, idents = "False")
Idents(seur_sub) <- "RNA_snn_res.0.3"
seur_sub <- subset(seur_sub, idents = c(0:4))


seur_sub <- NormalizeData(seur_sub, normalization.method = "LogNormalize", scale.factor = 10000)
seur_sub <- FindVariableFeatures(seur_sub, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(seur_sub)
seur_sub <- ScaleData(seur_sub, features = all.genes)

seur_sub <- RunPCA(seur_sub, features = VariableFeatures(object = seur_sub))

ElbowPlot(seur_sub)

seur_sub <- FindNeighbors(seur_sub, dims = 1:15)
seur_sub <- FindClusters(seur_sub, resolution = 0.3)

seur_sub <- RunUMAP(seur_sub, dims = 1:15)
DimPlot(seur_sub, reduction = "umap")
DimPlot(seur_sub,group.by="scrublet")

saveRDS(seur_sub, file="v8/Tdsort_soupx_v8.rds")
