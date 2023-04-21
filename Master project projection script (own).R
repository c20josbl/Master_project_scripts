## Master's project script projection

#Loading the required packages

library(Seurat)
library(gridExtra)
library(ggplot2)
library(clustree)
library(ggtree)
library(Nebulosa)

#Reading the ueRA and HD Seurat objects

ueRA = readRDS("C:/master_project/ueRA1-3-T3Dr_base")
DW=readRDS("C:/master_project/scPure2_HB6_UMAP3D.rds")

#Figure 2A (HD, to the left)

DimPlot(DW,pt.size=1,label=TRUE,label.size = 7,repel=TRUE)+
  NoLegend()

#Projection of ueRA dataset onto HD dataset

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

#Projection of HD dataset onto ueRA dataset

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

#Figure 2A (HD, middle)

DimPlot(query, group.by = "predicted.seurat_clusters", 
           reduction = "ref.umap", label = T, repel = T,pt.size=1.5,
           label.size = 7)+
  ggtitle("")+NoLegend()

#Table with predicted identities of the HD when projected onto ueRA

table(query$predicted.seurat_clusters)

#Differential gene expression analysis of HD with predicted cluster identities

Idents(query)="predicted.seurat_clusters"

query.marks = subset(FindAllMarkers(query, only.pos = TRUE, min.pct = 0.25),
                     subset = p_val_adj < 0.05)

#Figure 3A (HD, top panels)

plot_density(query,size=1,features=c("TNFRSF13C","TNFRSF13B"),
             joint=FALSE,combine=TRUE,reduction="ref.umap")

#Figure 2D (Dot plot)

DotPlot(query,features=c("IL6","TNF"))+xlab("Gene")+ylab("Cluster")

