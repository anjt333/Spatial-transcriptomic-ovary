rm(list=ls())
memory.limit(2000000)
library(Seurat)
library(ggplot2)
library(dplyr)

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
slidelist=c(slideA,slideB,slideC,slideD,slideE,slideF,slideG,slideI,slideJ,slideK,slideL,slideM,slideN,slideO)
for (i in 1:length(slidelist)){
  slidelist[[i]] <- SCTransform(slidelist[[i]], assay = "Spatial", verbose = FALSE)
  slidelist[[i]]=NormalizeData(slidelist[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  #slidelist[[i]]=FindVariableFeatures(slidelist[[i]], selection.method = "vst", nfeatures = 2000)
}

testAB.anchors <- FindIntegrationAnchors(object.list = slidelist, dims = 1:30)
testAB.integrated <- IntegrateData(anchorset = testAB.anchors, dims = 1:30)
DefaultAssay(testAB.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
testAB.integrated <- ScaleData(testAB.integrated, features = rownames(testAB.integrated))
testAB.integrated <- RunPCA(testAB.integrated, verbose = FALSE)
testAB.integrated <- FindNeighbors(testAB.integrated, dims = 1:30)
testAB.integrated <- FindClusters(testAB.integrated, resolution = 0.3)
testAB.integrated <- RunUMAP(testAB.integrated, dims = 1:30)
SpatialDimPlot(testAB.integrated,label = T,label.size = 2)
DimPlot(testAB.integrated,label = T)
My_levels <- c("stroma0","stroma1","immune","etheca","stroma3","injured","stroma1",
               "muscle","endo","oocyte","gran","itheca")
names(My_levels) <- levels(Idents(testAB.integrated))
testAB.integrated <- RenameIdents(testAB.integrated, My_levels)
Idents(testAB.integrated)=factor(Idents(testAB.integrated),
                                   levels=c("gran","itheca","etheca",
                                            "stroma1","stroma0","stroma3","oocyte",
                                            "muscle","endo","immune","injured"))
