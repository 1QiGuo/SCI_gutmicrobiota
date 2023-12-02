#change mitochondrio percent to 20 from 5 to reserve more cells
library("qs")
#install.packages("xlsx")
library("xlsx")
library("Seurat")
setwd("/fs/ess/PCON0022/guoqi/microbe_sci/rawdata")
load("./combined_raw.RData")
#calculates the percentage of counts originating from a set of features
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^mt-")
setwd("/fs/ess/PCON0022/guoqi/microbe_sci/rawdata/mito_20")
splitobject<-SplitObject(combined,split.by = "orig.ident")
sample1<-subset(splitobject[[1]], subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
vln_qc_sample1<-VlnPlot(sample1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(
  filename = "qc_vln_sample1.tiff",
  plot = vln_qc_sample1,
  device = "tiff",
  width = 10,
  height = 10,
  units = "in",
  dpi = 150
)
sample2<-subset(splitobject[[2]], subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 20)
vln_qc_sample2<-VlnPlot(sample2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(
  filename = "qc_vln_sample2.tiff",
  plot = vln_qc_sample2,
  device = "tiff",
  width = 10,
  height = 10,
  units = "in",
  dpi = 150
)
sample3<-subset(splitobject[[3]], subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 20)
vln_qc_sample3<-VlnPlot(sample3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(
  filename = "qc_vln_sample3.tiff",
  plot = vln_qc_sample3,
  device = "tiff",
  width = 10,
  height = 10,
  units = "in",
  dpi = 150
)
sample4<-subset(splitobject[[4]], subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 20)
vln_qc_sample4<-VlnPlot(sample4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(
  filename = "qc_vln_sample4.tiff",
  plot = vln_qc_sample4,
  device = "tiff",
  width = 10,
  height = 10,
  units = "in",
  dpi = 150
)
ncol(sample1)+ncol(sample2)+ncol(sample3)+ncol(sample4)
#cell number: 6160


#add condition
#injuy or not
sample1$spinalcord<-"SCI"
sample2$spinalcord<-"SCI"
sample3$spinalcord<-"SC"
sample4$spinalcord<-"SC"
#germ free or colonized
sample1$microbiota <-"Germ free"
sample2$microbiota <-"Colonized germ"
sample3$microbiota <-"Germ free"
sample4$microbiota <-"Colonized germ"
#---------------------------Normalization, we used lognormalization here.
identical(sample1@assays$RNA@counts,sample1@assays$RNA@data)#True, no any normalization
list_4sample_qc<-list(sample1=sample1,sample2=sample2,sample3=sample3,sample4=sample4)
# normalize and identify variable features for each dataset independently
list_4sample_qc <- lapply(X = list_4sample_qc, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

#-----------------------------------check whether we need to remove batch effect and integration
merge_4sample <-merge(sample1,c(sample2,sample3,sample4))
merge_4sample <- NormalizeData(merge_4sample)
merge_4sample <- FindVariableFeatures(merge_4sample, selection.method = "vst", nfeatures = 2000)
merge_4sample  <- ScaleData(merge_4sample , verbose = FALSE)
merge_4sample  <- RunPCA(merge_4sample , npcs = 30, verbose = FALSE)
merge_4sample  <- RunUMAP(merge_4sample , reduction = "pca", dims = 1:30)
merge_4sample  <- FindNeighbors(merge_4sample , reduction = "pca", dims = 1:30)
merge_4sample <- FindClusters(merge_4sample, resolution = 0.5)
merge_umap <- DimPlot(merge_4sample, reduction = "umap", group.by = "orig.ident")
merge_umap_spinalcord <- DimPlot(merge_4sample, reduction = "umap", group.by = "spinalcord")
merge_umap_microbiota <- DimPlot(merge_4sample, reduction = "umap", group.by = "microbiota")
merge_umap_cluster <- DimPlot(merge_4sample, reduction = "umap", group.by = "seurat_clusters")
setwd("/fs/ess/PCON0022/guoqi/microbe_sci/results/merge/mit20")
ggsave(
  filename = "umap_merge_origident.tiff",
  plot = merge_umap,
  device = "tiff",
  width = 10,
  height = 10,
  units = "in",
  dpi = 150
)
ggsave(
  filename = "umap_merge_microbiota.tiff",
  plot = merge_umap_microbiota,
  device = "tiff",
  width = 10,
  height = 10,
  units = "in",
  dpi = 150
)
ggsave(
  filename = "umap_merge_spinalcord.tiff",
  plot = merge_umap_spinalcord,
  device = "tiff",
  width = 10,
  height = 10,
  units = "in",
  dpi = 150
)
qsave(merge_4sample,"merge_object.qs")
#batch effect does not apparent but still has some. It can be observed by the umap labelled by spinal cord.
#Cells in different conditions do not mixed together.
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = list_4sample_qc)
anchors <- FindIntegrationAnchors(object.list = list_4sample_qc, anchor.features = features)
sci_integrate <- IntegrateData(anchorset = anchors)
DefaultAssay(sci_integrate) <- "integrated"
# Run the standard workflow for visualization and clustering
sci_integrate <- ScaleData(sci_integrate, verbose = FALSE)
sci_integrate <- RunPCA(sci_integrate, verbose = FALSE)
pca<-ElbowPlot(sci_integrate,ndims = 40)
setwd("/fs/ess/PCON0022/guoqi/microbe_sci/results/integrate/mit20")
ggsave(
  filename = "pca_elbowplot.tiff",
  plot = pca,
  device = "tiff",
  width = 10,
  height = 10,
  units = "in",
  dpi = 150
)
sci_integrate <- FindNeighbors(sci_integrate, reduction = "pca", dims = 1:30)
sci_integrate <- RunUMAP(sci_integrate, reduction = "pca", dims = 1:30)
umap_int_spinalcord <- DimPlot(temp, reduction = "umap", group.by = "spinalcord", label.size = 6)
umap_int_origident <- DimPlot(temp, reduction = "umap", group.by = "orig.ident", label.size = 6)
umap_int_microbiota <- DimPlot(temp, reduction = "umap", group.by = "microbiota", label.size = 6)

qsave(sci_integrate,"int_object.qs")
ggsave(
  filename = "umap_int_spinalcord.tiff",
  plot = umap_int_spinalcord,
  device = "tiff",
  width = 10,
  height = 10,
  units = "in",
  dpi = 150
)
ggsave(
  filename = "umap_int_origident.tiff",
  plot = umap_int_origident,
  device = "tiff",
  width = 10,
  height = 10,
  units = "in",
  dpi = 150
)
ggsave(
  filename = "umap_int_microbiota.tiff",
  plot = umap_int_microbiota,
  device = "tiff",
  width = 10,
  height = 10,
  units = "in",
  dpi = 150
)


temp <- FindClusters(sci_integrate, resolution = 0.3)
Idents(temp)<-temp$integrated_snn_res.0.3
qsave(temp,"./int_annotation.qs")
setwd('/fs/ess/PCON0022/guoqi/microbe_sci/results/integrate/mit20')
temp<-qread("./int_annotation.qs")
plot<-DimPlot(temp,reduction = "umap")
select.cells1<-CellSelector(plot=plot)
select.cells2<-CellSelector(plot=plot)
cluster20<-select.cells1
cluster21<-setdiff(select.cells2,select.cells1)
Idents(temp,cells=cluster20)<-"20"
Idents(temp,cells=cluster21)<-"21"
temp$clusters_f<-Idents(temp)
temp$clusters_f<-factor(temp$clusters_f,levels = c(levels(temp$seurat_clusters),"20","21"))
p2 <- DimPlot(temp, reduction = "umap", group.by = "clusters_f",label = T, label.size = 6)
ggsave(
  filename = "umap_int_clusters.tiff",
  plot = p2,
  device = "tiff",
  width = 10,
  height = 10,
  units = "in",
  dpi = 150
)

#Identify cluster-specific DEGs
setwd("/fs/ess/PCON0022/guoqi/microbe_sci/results/annotation/validation")
Idents(temp)<-temp$seurat_clusters
clusterDEG<-FindAllMarkers(temp)
clusterDEG_f<-clusterDEG[which(clusterDEG$p_val_adj<0.05),]
write.csv(clusterDEG_f,"cluster-specific DEGs.xlsx")
# write.xlsx(clusterDEG_f, "cluster-specific DEGs.xlsx", sheetName = "Sheet1", 
#            col.names = TRUE, row.names = TRUE, append = FALSE)
clusterDEG_f<-read.csv("cluster-specific DEGs.xlsx")
#update cluster-specific DEGs comparison of identity class to define markers for vs all other cells 
Idents(temp)<-temp$clusters_f
clusterDEG_v2<-FindAllMarkers(temp)
clusterDEG_v2_f<-clusterDEG_v2[which(clusterDEG_v2$p_val_adj<0.05),]
write.csv(clusterDEG_v2_f,"cluster-specific DEGs_v2.csv")

#----compare difference between version 1 and version 2
length(clusterDEG_f$gene[which(clusterDEG_f$cluster=='0')])
for(i in 1:19){
  logi<-identical(clusterDEG_v2_f$gene[which(clusterDEG_v2_f$cluster==i)],clusterDEG_f$gene[which(clusterDEG_f$cluster==i)])
  if(logi==F){
    print(i)
  }
}
setdiff(clusterDEG_v2_f$gene[which(clusterDEG_v2_f$cluster=='11')],clusterDEG_f$gene[which(clusterDEG_f$cluster=='11')])
identical(clusterDEG_v2_f$gene[which(clusterDEG_v2_f$cluster=='2')],clusterDEG_f$gene[which(clusterDEG_f$cluster=='2')])
