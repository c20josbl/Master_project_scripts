## Master's project ueRA1-3-T3Dr script

#Loading the required packages

library(Seurat)
library(xlsx)
library(ggplot2)
library(gridExtra)
library(plotrix)
library(dplyr)
library(Nebulosa)

#Reading in the integrated ueRA Seurat object

ueRA = readRDS("C:/master_project/ueRA1-3-T3Dr_base")

#Figure 1A and 2A (ueRA, right)

DimPlot(ueRA,pt.size=1,label=TRUE,label.size = 7)+
  NoLegend()

#Supplementary figure 1C

DimPlot(ueRA,pt.size=1,label=T,dims = c(2,3))

#Figure 1A (cells per cluster)

table(Idents(ueRA))
sum(table(Idents(ueRA)))

df=data.frame(cluster=c(0:6),count=c(5033,2765,2439,
                                     607,291,269,131))
df$fraction=df$count/sum(df$count)
df$ymax=cumsum(df$fraction)
df$ymin=c(0,head(df$ymax,n=-1))
df$labelPosition=(df$ymax+df$ymin)/2
df$label=paste0(df$count)
df$Cluster=factor(df$cluster)


ggplot(df,aes(ymax=ymax,ymin=ymin,xmax=4,
              xmin=0,fill=Cluster))+
  geom_rect()+coord_polar(theta="y")+xlim(c(0,4))+theme_void()+
  NoLegend()

#Figure 1B

table(ueRA$orig.ident)

table(Idents(ueRA),ueRA$orig.ident)
Patient=rep(c("ueRA1","ueRA2","ueRA3"),7)
cluster=c(rep("0",3),rep("1",3),rep("2",3),
          rep("3",3),rep("4",3),rep("5",3),
          rep("6",3))
prop=c((1787/5033),(1500/5033),(1746/5033),(1046/2765),(956/2765),
       (763/2765),(771/2439),(1066/2439),
       (602/2439),(231/607),(258/607),(118/607),
       (104/291),(101/291),(86/291),(89/269),(81/269),(99/269),
       (35/131),(55/131),(41/131))
Igp<-data.frame(Patient,cluster,prop)

ggplot(Igp,aes(fill=Patient,y=prop,
               x=cluster))+
  geom_bar(position="fill",stat="identity")+
  xlab("Cluster")+ylab("Proportion")+
  theme_bw()+coord_flip()+
  scale_fill_brewer(palette="Pastel2")+
  theme(legend.title=element_text(size=15,face="bold"),
        legend.text = element_text(size=15),
        axis.title = element_text(size=15,face="bold"),
        axis.text = element_text(size=15))


#Sorting cells based on heavy chain expression

igm = WhichCells(object = ueRA, expression = Heavy_chain == "IGHM")

igg = c(WhichCells(object = ueRA, expression = Heavy_chain == "IGHG1"), 
        WhichCells(object = ueRA, expression = Heavy_chain == "IGHG2"), 
        WhichCells(object = ueRA, expression = Heavy_chain == "IGHG3"),
        WhichCells(object = ueRA, expression = Heavy_chain == "IGHG4")) 

iga = c(WhichCells(object = ueRA, expression = Heavy_chain == "IGHA1"), 
        WhichCells(object = ueRA, expression = Heavy_chain == "IGHA2"))

igd = WhichCells(object = ueRA, expression = Heavy_chain == "IGHD")

ige = WhichCells(object = ueRA, expression = Heavy_chain == "IGHE")

#Proportions of heavy chain per cluster

IgM=CellsByIdentities(ueRA,cells=igm)
str(IgM)
IgG=CellsByIdentities(ueRA,cells=igg)
str(IgG)
IgA=CellsByIdentities(ueRA,cells=iga)
str(IgA)
IgE=CellsByIdentities(ueRA,cells=ige)
str(IgE)
IgD=CellsByIdentities(ueRA,cells=igd)
str(IgD)

#for cluster 0
total0=3381+340+659+32
propM0=(3381/total0)*100
propG0=(340/total0)*100
propA0=(659/total0)*100
propE0=0
propD0=(32/total0)*100


#for cluster 1
total1=144+1071+401+2+15
propM1=(144/total1)*100
propG1=(1071/total1)*100
propA1=(401/total1)*100
propE1=(2/total1)*100
propD1=(15/total1)*100


#for cluster 2
total2=151+854+684+18
propM2=(151/total2)*100
propG2=(854/total2)*100
propA2=(684/total2)*100
propE2=0
propD2=(18/total1)*100


#for cluster 3
total3=279+124+65+9
propM3=(279/total3)*100
propG3=(124/total3)*100
propA3=(65/total3)*100
propE3=0
propD3=(9/total1)*100


#for cluster 4
total4=105+29+39+4
propM4=(105/total4)*100
propG4=(29/total4)*100
propA4=(39/total4)*100
propE4=0
propD4=(4/total1)*100


#for cluster 5
total5=137+27+52+4
propM5=(137/total4)*100
propG5=(27/total4)*100
propA5=(52/total4)*100
propE5=0
propD5=(4/total1)*100


#for cluster 6
total6=107+9+1+5
propM6=(107/total4)*100
propG6=(9/total4)*100
propA6=(1/total4)*100
propE6=0
propD6=(5/total1)*100

#Figure 1E

Ig=rep(c("IgM","IgA","IgG","IgE","IgD"),7)
cluster=c(rep("0",5),rep("1",5),rep("2",5),
          rep("3",5),rep("4",5),rep("5",5),
          rep("6",5))
prop=c(propM0,propA0,propG0,propE0,propD0,propM1,propA1,propG1,
       propE1,propD1,propM2,propA2,
       propG2,propE2,propD2,propM3,propA3,propG3,
       propE3,propD3,propM4,propA4,propG4,propE4,propD4,
       propM5,propA5,propG5,propE5,propD5,
       propM6,propA6,propG6,propE6,propD6)
Igp<-data.frame(Ig,cluster,prop)

ggplot(Igp,aes(fill=Ig,y=prop,
               x=cluster))+geom_bar(position="fill",
                                    stat="identity")+xlab("Cluster")+ylab("Proportion")+
  theme_bw()+coord_flip()+
  scale_fill_brewer(palette="Pastel2")+
  theme(legend.title=element_text(size=15,face="bold"),
        legend.text = element_text(size=15),
        axis.title = element_text(size=15,face="bold"),
        axis.text = element_text(size=15))


#Figure 1D

cd27 = WhichCells(ueRA, expression = CD27 > 0) 
DimPlot(ueRA,cells.highlight = cd27,cols.highlight = "steelblue",
        pt.size = 1,label=T,label.size = 7)+NoLegend()


#Figure 1C

library(clustree)
library(ggtree)

tree = BuildClusterTree(ueRA)
myphytree = Tool(tree, slot = "BuildClusterTree")
ggtree(myphytree)+geom_tiplab(size=10)+theme_tree()

#Supplementary figure 1B

up_degs = read.xlsx("C:/master_project/ueRA1-3-T3Dr_DEGs.xlsx", sheetName = "Upregulated")

DoHeatmap(ueRA, features=up_degs$gene)+theme(axis.text.y = element_blank())

clust0 = subset(up_degs, cluster == "0") 
clust1 = subset(up_degs, cluster == "1") 
clust2 = subset(up_degs, cluster == "2") 
clust3 = subset(up_degs, cluster == "3")
clust4 = subset(up_degs, cluster == "4") 
clust5 = subset(up_degs, cluster == "5") 
clust6 = subset(up_degs, cluster == "6") 

#Figure 1A (differentially expressed genes per cluster)

df=data.frame(cluster=c(0:6),count=c(303,296,410,
                                     453,288,224,329))
df$fraction=df$count/sum(df$count)
df$ymax=cumsum(df$fraction)
df$ymin=c(0,head(df$ymax,n=-1))
df$labelPosition=(df$ymax+df$ymin)/2
df$label=paste0("\n DEGs:",df$count)
df$Cluster=factor(df$cluster)


ggplot(df,aes(ymax=ymax,ymin=ymin,xmax=4,
              xmin=0,fill=Cluster))+
  geom_rect()+coord_polar(theta="y")+xlim(c(0,4))+theme_void()+
  theme(legend.title = element_text(size=15,face="bold"),
        legend.text=element_text(size=15))


#Supplementary figure 2A

DotPlot(ueRA,features=c("CXCR3","MTOR","TNFSF11"),assay="RNA")

#Figure 2B

genes=c("TCL1A", "APBB2",
        "MIR181A1HG","IKZF2","PLPP5","SELL","ZBTB16",
        "PCDH9","IL4R","IL21R","FCER2","CD38","CD27",
        "PARP15","COBLL1",
        "TCF4","CD1C",
        "RPS18","LINC01857","TEX9","COCH","MUC16","TOX","ITGAX","TBX21")

ueRA$groups = factor(Idents(ueRA))
ueRA$groups = factor(ueRA$groups, levels = c("6","0","1","2","3","4","5"))

DotPlot(ueRA, features = genes,group.by = "groups")+
  RotatedAxis()+coord_flip()+xlab("Genes")+
  ylab("Cluster")


#Figure 2A (ueRA, bottom panels)

plot_density(ueRA,features = c("TNFRSF13C","TNFRSF13B"),reduction = "umap",
             joint = FALSE,combine = TRUE, size = 1)

#Figure 2C

p1=VlnPlot(ueRA,features="IL6")+xlab("Cluster")+NoLegend()
p2=VlnPlot(ueRA,features="TNF")+xlab("Cluster")+NoLegend()
grid.arrange(p1,p2,ncol=2)


#3D UMAP (used for supplementary figure 1D)

TS_plot_3D <- function (seurat_object, dims = 1:30) {
  require(plotly)
  require(scales)
  temp1 <- RunUMAP(object = seurat_object, dims = dims, n.components = 3,
                   verbose = F)
  df <- data.frame(umap1 = temp1@reductions$umap@cell.embeddings[,
                                                                 1], umap2 = temp1@reductions$umap@cell.embeddings[, 2],
                   umap3 = temp1@reductions$umap@cell.embeddings[, 3], cell = Idents(temp1))
  plot_ly(df, x = ~umap1, y = ~umap2, z = ~umap3, color = ~cell,
          colors = hue_pal()(length(levels(df$cell)))) %>% add_markers() %>%
    layout(scene = list(xaxis = list(title = "umap 1"), yaxis = list(title = "umap 2"),
                        zaxis = list(title = "umap 3")))
}

TS_plot_3D(ueRA)
