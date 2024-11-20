library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(ggsci)
library(dplyr)
library(stringr)
rm(list=ls())

################Singlecell data##################
ovarysc_reference <- readRDS("ovary_reference.rds")
DimPlot(ovarysc_reference,label = T)+scale_color_d3(palette = "category20",alpha=0.4)
My_levels <- c("Stroma","Stroma","Stroma","Immunue","Granulosa","Granulosa","Theca","Endothelial",
               "Granulosa","Endothelial","Immunue","Granulosa","Perivascular","Granulosa","Granulosa",
               "Endothelial","Perivascular","Immunue")
names(My_levels) <- levels(Idents(ovarysc_reference))
ovarysc_reference <- RenameIdents(ovarysc_reference, My_levels)
Idents(ovarysc_reference)=factor(Idents(ovarysc_reference),levels=c("Granulosa","Theca","Stroma","Perivascular",
                                                                    "Endothelial","Immunue"))
mycolor=c("#C64F4F","#66CDAA","#87CEEB","#CD853F","#F7B6D2FF","#FFD700")
DimPlot(ovarysc_reference,label = T,cols=mycolor,label.size = 4,label.box = T)
gene=c("AMH","COL1A1","DCN","CXCR4","CD53","RGS5","TAGLN","VWF","CLDN5","CD34")
for(i in gene){
  pdf(paste0(i,".umap.pdf"),height=5,width=5)
  p=FeaturePlot(ovarysc_reference,features = i)
  print(p)
  dev.off()
}
