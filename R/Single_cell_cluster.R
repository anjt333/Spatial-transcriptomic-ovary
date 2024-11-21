Single_cell_cluster.Rlibrary(Seurat)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(ggsci)
library(dplyr)
library(stringr)
rm(list=ls())
ovarysc_reference <- readRDS("ovary_reference.rds")#2019NC.reference
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

################mapping##################
testAB.integrated=readRDS("allnew.11.rds")
cortex=testAB.integrated
anchors <- FindTransferAnchors(reference = ovarysc_reference, query = cortex)
predictions.assay <- TransferData(anchorset = anchors, refdata = Idents(ovarysc_reference), prediction.assay = TRUE)
cortex[["predictions"]] <- predictions.assay
DefaultAssay(cortex) <- "predictions"
FeaturePlot(cortex, features = c("0","1","2","3","4","5","6","7","8",
                                 "9","10","11","12","13","14","15","16","17"))

gene=c("0","2","4","5","6","10","12","13","15")
for(i in gene){
  pdf(paste0(i,".umap.pdf"),height=5,width=5)
  p=FeaturePlot(cortex,features = i)
  print(p)
  dev.off()
}
