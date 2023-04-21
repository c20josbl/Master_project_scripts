##Master's project ueRA1-3-T subclustering cluster 1

#Loading the required packages

library(Seurat)
library(dplyr)
library(xlsx)
library(ggplot2)
library(Nebulosa)

#Reading the integrated ueRA object and subsetting cluster 1

ueRA = readRDS("C:/master_project/ueRA1-3-T3Dr_base")
ueRA.1 = subset(ueRA, subset = seurat_clusters == "1")

#SCTransform and principal component analysis

set.seed(1001)

ueRA.1 = SCTransform(ueRA.1, variable.features.n = dim(ueRA.1)[1], assay = "RNA")

DefaultAssay(ueRA.1)="integrated"
ueRA.1 = RunPCA(ueRA.1, features = VariableFeatures(ueRA.1),
                verbose=FALSE)
ElbowPlot(ueRA.1, reduction = "pca", ndims = 40)

#Finding the appropriate resolution

ueRA.1 = FindNeighbors(ueRA.1, dims = 1:25)

library(clustree)
res = seq.int(0.2, 0.8, 0.1) 
for (i in res){
  ueRA.1 = FindClusters(object = ueRA.1, resolution = i)
}
clustree(ueRA.1@meta.data, prefix = "integrated_snn_res.") 

#Unsupervised clustering

ueRA.1 = FindClusters(ueRA.1, resolution = 0.2)

#UMAP dimensionality reduction

ueRA.1 = RunUMAP(ueRA.1, dims = 1:25)

#Re-naming clusters

c1.0 = WhichCells(ueRA.1, idents = "1.0")
c1.1 = WhichCells(ueRA.1, idents = "1.1")

Idents(ueRA.1, cells = c1.0) = "1.0"
Idents(ueRA.1, cells = c1.1) = "1.1"

#Differential gene expression analysis of cluster 1 subclusters

ueRA1.marks = subset(FindAllMarkers(ueRA.1, only.pos = TRUE, min.pct = 0.25),
                     subset = p_val_adj < 0.05)

#Saving the subclustered cluster 1 Seurat object

saveRDS(ueRA.1, file = "C:/master_project/ueRA1-3-T3Dr_cluster1_base")

#Reading integrated ueRA Seurat object and the subclustered cluster 1
#Seurat object

ueRA = readRDS("C:/master_project/ueRA1-3-T3Dr_base")
ueRA.1 = readRDS("C:/master_project/ueRA1-3-T3Dr_cluster1_base")

#Including cluster 1 subclusters in the integrated ueRA object

ueRA.1 = ueRA
Idents(ueRA.1, cells = c1.0) = "1.0"
Idents(ueRA.1, cells = c1.1) = "1.1"

Idents(ueRA.1) = factor(Idents(ueRA.1), 
                         levels = c("0", "1.0","1.1",
                                    "2","3","4","5","6"))
#Figure 3B

DimPlot(ueRA.1,label=TRUE,pt.size=1.5,
        cols = c("#F8766D","#CD9600","steelblue","#00BE67","#00BFC4",
                          "#00A9FF","#C77CFF","#FF61CC"),
                          label.size = 7)+NoLegend()

#Figure 3C

genes=c("CD27","SELL","HOPX","PDE4D")

DoHeatmap(ueRA.1, features=genes,raster=FALSE,
        group.colors=c("#F8766D","#CD9600","steelblue","#00BE67","#00BFC4",
                                "#00A9FF","#C77CFF","#FF61CC"))+
                                  NoLegend()

#Analysing heavy chains present in cluster 1 subclusters

#Removing cells with heavy chain=NA from the analysis
toRemove=which(is.na(ueRA.1$Heavy_chain))
df=data.frame(toRemove)
ueRA2=ueRA.1

ueRA2$Clusters=Idents(ueRA2)
Idents(ueRA2) = ueRA2$predicted.id
Idents(object=ueRA2, cells = row.names(df)) = "NA"
ueRA2 = subset(ueRA2, idents = "B cell")
Idents(ueRA2)=ueRA2$Clusters

#Table with heavy chain counts
table(ueRA2$Heavy_chain,Idents(ueRA2))
