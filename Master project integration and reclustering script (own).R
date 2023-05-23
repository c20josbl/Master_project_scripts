##Master project integration and reclustering script

#Loading the required packages

library(Seurat)
library(xlsx)
library(ggplot2)
library(gridExtra)
library(plotrix)
library(dplyr)

set.seed(1001)

#Reading in separate Seurat Objects for patients 1-3

ueRA1 = readRDS("C:/master_project/ueRA1-T_base")
ueRA2 = readRDS("C:/master_project/ueRA2-T_base")
ueRA3 = readRDS("C:/master_project/ueRA3-T_base")

#Adding patient information to metadata

ueRA1 = AddMetaData(ueRA1,metadata="ueRA1",col.name = "patient")
ueRA2 = AddMetaData(ueRA2,metadata="ueRA2",col.name = "patient")
ueRA3 = AddMetaData(ueRA3,metadata="ueRA3",col.name = "patient")

#Adding heavy chain from VDJ-sequencing analysis for each patient 

vdj1 = read.csv("C:/master_project/RA1/RA1_filtered_contig_annotations.csv",sep=",") 
vdj1 = vdj1 %>% filter(productive == "true", chain == "IGH", is_cell == "true")
nrow(vdj1[(duplicated(vdj1$barcode) | duplicated(vdj1$barcode, fromLast = TRUE)), ])
nrow(vdj1 %>% filter(c_gene == ""))
vdjdubs = vdj1[(duplicated(vdj1$barcode) | duplicated(vdj1$barcode, fromLast = TRUE)), ]
vdj1 = vdj1[!(duplicated(vdj1$barcode) | duplicated(vdj1$barcode, fromLast = TRUE)), ]
vdj1 = vdj1 [!(vdj1$c_gene==""),]
vdj1 = vdj1[,c("barcode","c_gene","reads","umis")]
names(vdj1)[names(vdj1) == "c_gene"] = "Heavy_chain"
row.names(vdj1) = vdj1[,1]
ueRA1 = AddMetaData(ueRA1, metadata = vdj1) 

vdj2 = read.csv("C:/master_project/RA2/RA2_filtered_contig_annotations.csv",sep=",") 
vdj2 = vdj2 %>% filter(productive == "true", chain == "IGH", is_cell == "true")
nrow(vdj2[(duplicated(vdj2$barcode) | duplicated(vdj2$barcode, fromLast = TRUE)), ])
nrow(vdj2 %>% filter(c_gene == ""))
vdjdubs = vdj2[(duplicated(vdj2$barcode) | duplicated(vdj2$barcode, fromLast = TRUE)), ]
vdj2 = vdj2[!(duplicated(vdj2$barcode) | duplicated(vdj2$barcode, fromLast = TRUE)), ]
vdj2 = vdj2 [!(vdj2$c_gene==""),]
vdj2 = vdj2[,c("barcode","c_gene","reads","umis")]
names(vdj2)[names(vdj2) == "c_gene"] = "Heavy_chain"
row.names(vdj2) = vdj2[,1]
ueRA2 = AddMetaData(ueRA2, metadata = vdj2) 

vdj3 = read.csv("C:/master_project/RA3/RA3_filtered_contig_annotations.csv",sep=",") 
vdj3 = vdj3 %>% filter(productive == "true", chain == "IGH", is_cell == "true")
nrow(vdj3[(duplicated(vdj3$barcode) | duplicated(vdj3$barcode, fromLast = TRUE)), ])
nrow(vdj3 %>% filter(c_gene == ""))
vdjdubs = vdj3[(duplicated(vdj3$barcode) | duplicated(vdj3$barcode, fromLast = TRUE)), ]
vdj3 = vdj3[!(duplicated(vdj3$barcode) | duplicated(vdj3$barcode, fromLast = TRUE)), ]
vdj3 = vdj3 [!(vdj3$c_gene==""),]
vdj3 = vdj3[,c("barcode","c_gene","reads","umis")]
names(vdj3)[names(vdj3) == "c_gene"] = "Heavy_chain"
row.names(vdj3) = vdj3[,1]
ueRA3 = AddMetaData(ueRA3, metadata = vdj3) 

#Perform integration

sample_list = list(ueRA1,ueRA2,ueRA3)

features <- SelectIntegrationFeatures(object.list = sample_list,nfeatures = 3000)

sample_list <- PrepSCTIntegration(object.list = sample_list, anchor.features = features)

ueRA.anchors <- FindIntegrationAnchors(object.list = sample_list, normalization.method = "SCT",
                                       anchor.features = features)

ueRA.combined.sct <- IntegrateData(anchorset = ueRA.anchors, normalization.method = "SCT")

#Dimensionality reduction of integrated object

ueRA.combined.sct <- RunPCA(ueRA.combined.sct, verbose = FALSE)
ueRA.combined.sct <- RunUMAP(ueRA.combined.sct, reduction = "pca", dims = 1:30)

#Re-cluster integrated object

set.seed(1001)


#Finding appropriate resolution

ueRA = FindNeighbors(ueRA, dims = 1:30)

library(clustree)
res = seq.int(0.1, 0.5, 0.05) 
for (i in res){
  ueRA = FindClusters(object = ueRA, resolution = i)
}

#Supplementary figure 1A

clustree(ueRA@meta.data, prefix = "integrated_snn_res.") 

#Unsupervised clustering

ueRA = FindClusters(ueRA, resolution = 0.15)

#UMAP dimensionality reduction

ueRA <- RunUMAP(ueRA, reduction = "pca", dims = 1:30,
                n.components = 3L)

#Differential gene expression alaysis 

ueRA.marks = subset(FindAllMarkers(ueRA, only.pos = FALSE, min.pct = 0.25),
                    subset = p_val_adj < 0.05)
pos.marks = subset(ueRA.marks, subset = avg_log2FC > 0.0) 
neg.marks = subset(ueRA.marks, subset = avg_log2FC < 0.0)

#Saving the differentially expressed genes and the integrated Seurat object

write.xlsx(pos.marks, "C:/master_project/ueRA1-3-T3Dr_DEGs.xlsx", sheetName = "Upregulated", col.names = TRUE, row.names = FALSE)
write.xlsx(neg.marks, "C:/master_project/ueRA1-3-T3Dr_DEGs.xlsx", sheetName = "Downregulated", col.names = TRUE, row.names = FALSE, append = TRUE)

saveRDS(ueRA, file = "C:/master_project/ueRA1-3-T3Dr_base")

