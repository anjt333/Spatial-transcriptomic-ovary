library(Seurat)
library(ggplot2)
rm(list=ls())
memory.limit(2000000)
slideA=Load10X_Spatial("E:/spatial transcriptome/test.yanjie/yanjie.Asample",
                       filename = "filtered_feature_bc_matrix.h5",
                       assay = "Spatial", slice = "A_2020001_2", filter.matrix = TRUE, to.upper = FALSE)
slideB=Load10X_Spatial("E:/spatial transcriptome/test.yanjie/yanjie.Bsample",
                       filename = "filtered_feature_bc_matrix.h5",
                       assay = "Spatial", slice = "B_2020001_2", filter.matrix = TRUE, to.upper = FALSE)
slideC=Load10X_Spatial("E:/spatial transcriptome/test.yanjie/yanjie.Csample",
                       filename = "filtered_feature_bc_matrix.h5",
                       assay = "Spatial", slice = "C_S20006_2", filter.matrix = TRUE, to.upper = FALSE)
slideD=Load10X_Spatial("E:/spatial transcriptome/test.yanjie/yanjie.Dsample",
                       filename = "filtered_feature_bc_matrix.h5",
                       assay = "Spatial", slice = "D_S20006_2", filter.matrix = TRUE, to.upper = FALSE)
slideE=Load10X_Spatial("E:/spatial transcriptome/test.yanjie/yanjie.Esample",
                       filename = "filtered_feature_bc_matrix.h5",
                       assay = "Spatial", slice = "E_S20006_2", filter.matrix = TRUE, to.upper = FALSE)
slideF=Load10X_Spatial("E:/spatial transcriptome/test.yanjie/yanjie.Fsample",
                       filename = "filtered_feature_bc_matrix.h5",
                       assay = "Spatial", slice = "F_2020001_4", filter.matrix = TRUE, to.upper = FALSE)
slideG=Load10X_Spatial("E:/spatial transcriptome/test.yanjie/yanjie.Gsample",
                       filename = "filtered_feature_bc_matrix.h5",
                       assay = "Spatial", slice = "G_2020001_5", filter.matrix = TRUE, to.upper = FALSE)
slideI=Load10X_Spatial("E:/spatial transcriptome/test.yanjie/yanjie.Isample",
                       filename = "filtered_feature_bc_matrix.h5",
                       assay = "Spatial", slice = "I_S21004_10", filter.matrix = TRUE, to.upper = FALSE)
slideJ=Load10X_Spatial("E:/spatial transcriptome/test.yanjie/yanjie.Jsample",
                       filename = "filtered_feature_bc_matrix.h5",
                       assay = "Spatial", slice = "J_S21006_1", filter.matrix = TRUE, to.upper = FALSE)
slideK=Load10X_Spatial("E:/spatial transcriptome/test.yanjie/yanjie.Ksample",
                       filename = "filtered_feature_bc_matrix.h5",
                       assay = "Spatial", slice = "K_S21006_1", filter.matrix = TRUE, to.upper = FALSE)
slideL=Load10X_Spatial("E:/spatial transcriptome/test.yanjie/yanjie.Lsample",
                       filename = "filtered_feature_bc_matrix.h5",
                       assay = "Spatial", slice = "L_S21004_13", filter.matrix = TRUE, to.upper = FALSE)
slideM=Load10X_Spatial("E:/spatial transcriptome/test.yanjie/yanjie.Msample",
                       filename = "filtered_feature_bc_matrix.h5",
                       assay = "Spatial", slice = "M_S21004_1", filter.matrix = TRUE, to.upper = FALSE)
slideN=Load10X_Spatial("E:/spatial transcriptome/test.yanjie/yanjie.Nsample",
                       filename = "filtered_feature_bc_matrix.h5",
                       assay = "Spatial", slice = "N_S21004_1", filter.matrix = TRUE, to.upper = FALSE)
slideO=Load10X_Spatial("E:/spatial transcriptome/test.yanjie/yanjie.Osample",
                       filename = "filtered_feature_bc_matrix.h5",
                       assay = "Spatial", slice = "O_S21004_1", filter.matrix = TRUE, to.upper = FALSE)
slideS=Load10X_Spatial("E:/spatial transcriptome/spatialslide.batch2/POF.sample/S_2022003_2",
                       filename = "filtered_feature_bc_matrix.h5",
                       assay = "Spatial", slice = "S_2022003_2", filter.matrix = TRUE, to.upper = FALSE)
slidelist=c(slideA,slideB,slideC,slideD,slideE,slideF,slideG,slideI,slideJ,slideK,slideL,slideM,slideN,slideO,slideS)
for (i in 1:length(slidelist)){
  slidelist[[i]] <- SCTransform(slidelist[[i]], assay = "Spatial", verbose = FALSE)
  slidelist[[i]]=NormalizeData(slidelist[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  #slidelist[[i]]=FindVariableFeatures(slidelist[[i]], selection.method = "vst", nfeatures = 2000)
}

testAB.anchors <- FindIntegrationAnchors(object.list = slidelist, dims = 1:30)
rm(list=ls()[! ls() %in% c("testAB.anchors")])
AtoS <- IntegrateData(anchorset = testAB.anchors, dims = 1:30)
DefaultAssay(AtoS) <- "integrated"
# Run the standard workflow for visualization and clustering
AtoS <- ScaleData(AtoS, features = rownames(testAB.integrated))
AtoS <- RunPCA(AtoS, verbose = FALSE)
AtoS <- FindNeighbors(AtoS, dims = 1:30)
AtoS <- FindClusters(AtoS, resolution = 0.5)
AtoS <- RunUMAP(AtoS, dims = 1:30)
DimPlot(AtoS,label = T)

library(stringr)
Cluster=subset(AtoS,idents=c("stroma0","stroma1","stroma3"))
split_b<-str_split(rownames(Cluster@meta.data),"_")
Idents(Cluster)<-sapply(split_b,"[",2)
Idents(Cluster)
#My_levels <- c("slideA","slideB","slideC","slideD","slideE","slideF","slideG","slideI","slideJ","slideK","slideL","slideM","slideN","slideO","slideS")
My_levels <- c("normal","normal","normal","normal","normal","normal","normal","normal",
               "normal","normal","normal","normal","normal","normal","slideS")
names(My_levels) <- levels(Idents(Cluster))
Cluster<- RenameIdents(Cluster, My_levels)
Idents(Cluster)
DefaultAssay(Cluster)<-"Spatial"
Cluster <- NormalizeData(Cluster,normalization.method = "LogNormalize", scale.factor = 10000)
Cluster <- ScaleData(Cluster, features = rownames(Cluster))
VlnPlot(Cluster,features = c("WFIKKN2","WISP2","GAPDH"),pt.size=0)
VlnPlot(Cluster, features = "WFIKKN2", pt.size = 0,cols=c("#B1F4CF","#A3A0EF")) + NoLegend()+
  scale_fill_manual(values=c("#B1F4CF","#A3A0EF"))+ggtitle("WFIKKN2")+theme(plot.title = element_text(hjust = 0.5))#生成WFIKKN2.vlnplot.normalvsPOF.pdf
SpatialFeaturePlot(AtoS, images="S_2022003_2",features = "WFIKKN2",crop=F,pt.size.factor=1.5,alpha = c(0, 1))+
  scale_fill_gradientn(colors = rev(brewer.pal(11,"Spectral")),limits=c(-1,2))
