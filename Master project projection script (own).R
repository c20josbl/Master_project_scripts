## Master's project script projection

library(Seurat)
library(gridExtra)
library(ggplot2)
library(clustree)
library(ggtree)
library(Nebulosa)


ueRA = readRDS("C:/master_project/ueRA1-3-T3Dr_base")
DW=readRDS("C:/master_project/scPure2_HB6_UMAP3D.rds")

DimPlot(DW,pt.size=1,label=TRUE,label.size = 7,repel=TRUE)+
  NoLegend()

set.seed(1001)

ref <- RunUMAP(DW, 
               dims = 1:14, reduction = "pca", 
               return.model = T,n.components = 3L)


anchors <- FindTransferAnchors(reference = ref, 
                               query = ueRA, dims = 1:14,
                               reference.reduction = "pca",
                               reference.assay = "RNA",
                               query.assay = "RNA")
query <- MapQuery(anchorset = anchors,
                  query = ueRA, reference = ref,
                  refdata = list(cluster = "cluster"), 
                  reference.reduction = "pca", reduction.model = "umap",
                  reference.dims = 1:14,query.dims = 1:30)

table(query$seurat_clusters,query$predicted.cluster)

#Using our data as reference

set.seed(1001)
ref <- RunUMAP(ueRA, 
               dims = 1:30, reduction = "pca", 
               return.model = T,n.components = 3L)


anchors <- FindTransferAnchors(reference = ref, 
                               query = DW, dims = 1:30,
                               reference.reduction = "pca",
                               reference.assay = "integrated",
                               query.assay = "RNA",
                               features = row.names(ref[["integrated"]]),
                               normalization.method = "SCT")
query <- MapQuery(anchorset = anchors,
                  query = DW, reference = ref,
                  refdata = list(seurat_clusters = "seurat_clusters"), 
                  reference.reduction = "pca", reduction.model = "umap",
                  reference.dims = 1:14,query.dims = 1:30)

DimPlot(query, group.by = "predicted.seurat_clusters", 
           reduction = "ref.umap", label = T, repel = T,pt.size=1.5,
           label.size = 7)+
  ggtitle("")+NoLegend()

table(query$predicted.seurat_clusters)

#Checking expression of BAFF-R, IL-6 and TNF-alpha

Idents(query)="predicted.seurat_clusters"

query.marks = subset(FindAllMarkers(query, only.pos = TRUE, min.pct = 0.25),
                     subset = p_val_adj < 0.05)

plot_density(query,size=1,features=c("TNFRSF13C","TNFRSF13B"),
             joint=FALSE,combine=TRUE,reduction="ref.umap")

DotPlot(query,features=c("IL6","TNF"))+xlab("Gene")+ylab("Cluster")

