rm(list=ls())
memory.limit(2000000)
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
slideM=Load10X_Spatial("E:/spatial transcriptome/test.yanjie/yanjie.Msample",
                       filename = "filtered_feature_bc_matrix.h5",
                       assay = "Spatial", slice = "M_S21004_1", filter.matrix = TRUE, to.upper = FALSE)
slideO=Load10X_Spatial("E:/spatial transcriptome/test.yanjie/yanjie.Osample",
                       filename = "filtered_feature_bc_matrix.h5",
                       assay = "Spatial", slice = "O_S21004_1", filter.matrix = TRUE, to.upper = FALSE)

Mbarcode=read.csv("M.cluster.an.csv",header=T)
Obarcode=read.csv("O.cluster.an.csv",header=T)
slideM=subset(slideM,cells=Mbarcode$Barcode)
Idents(slideM)=Mbarcode$M
data.input1  <- slideM@assays$Spatial@data
which(data.input1["GAPDH",]==0)
data.input1=data.input1[,-580]
data.input1 <- sweep(data.input1, 2, data.input1["GAPDH",], `/`)
colnames(data.input1)=paste0(Mbarcode$Barcode[-580],"_M")
identity1 = data.frame(group =Idents(slideM)[-580],row.names = colnames(data.input1))
slideO=subset(slideO,cells=Obarcode$Barcode)
Idents(slideO)=Obarcode$O
data.input2  <- slideO@assays$Spatial@data
which(data.input2["GAPDH",]==0)
data.input2 <- sweep(data.input2, 2, data.input2["GAPDH",], `/`)
colnames(data.input2)=paste0(Obarcode$Barcode,"_O")
identity2 = data.frame(group =Idents(slideO),row.names = colnames(data.input2))
data.input=cbind(data.input1,data.input2)
identity=rbind(identity1,identity2)

unique(identity$group)
identity$group=factor(identity$group,levels=c("oocyte.AF1","gran2.AF1","cgran.AF1","mgran.AF1",
                                              "itheca.AF1","etheca.AF1",
                                              "gran2.AF2","cgran.AF2","mgran.AF2",
                                              "itheca.AF2","etheca.AF2",
                                              "oocyte.AF3","gran2.AF3","cgran.AF3","mgran.AF3",
                                              "itheca.AF3","etheca.AF3",
                                              "AF4","cluster9","stroma0"))
cellchat <- createCellChat(data.input)


###############Network#######################
cellchat
summary(cellchat)
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
rm(list=ls()[! ls() %in% c("cellchat")])
cellchat <- setIdent(cellchat, ident.use = "labels")
levels(cellchat@idents) 
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- CellChatDB.human 
colnames(CellChatDB$interaction)
CellChatDB$interaction[1:4,1:4]
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use
unique(CellChatDB$interaction$annotation)
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)  
mycomputeCommunProb <-edit(computeCommunProb)
environment(mycomputeCommunProb) <- environment(computeCommunProb)
cellchat <- mycomputeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat@netP$pathways
head(cellchat@LR$LRsig)
cellchat@netP$pathways
levels(cellchat@idents) 

pathways.show <- cellchat@netP$pathways
netVisual_aggregate(cellchat, signaling = c("AMH"))
netVisual_aggregate(cellchat, signaling = c("AMH"),top=0.1)
setwd("C:/Users/Dell/Desktop")
pdf("pathwayall.top0.1.pdf")
for(i in 1:length(pathways.show)){
  netVisual_aggregate(cellchat, signaling = pathways.show[i],top=0.1)
}
dev.off()

netVisual_heatmap(cellchat, signaling = c("AMH"), color.heatmap = "Reds")
pdf("pathwayheatmap.pdf")
for(i in 1:length(pathways.show)){
  p=netVisual_heatmap(cellchat, signaling = pathways.show[i],color.heatmap = "Reds")
  print(p)
}
dev.off()

netVisual_heatmap(cellchat, signaling = pathways.show) 
plotGeneExpression(cellchat,signaling = pathways.show)
plotGeneExpression(cellchat,signaling = "BMP")
netVisual_bubble(cellchat,signaling = pathways.show, remove.isolate = FALSE)

cellchat <- netAnalysis_computeCentrality(cellchat,slot.name = "netP")
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show)
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",height = 20)
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",height = 20)

Mco=read.table("M.cluster.an.coordinate.txt",header=T,row.names = 1)
Oco=read.table("O.cluster.an.coordinate.txt",header=T,row.names = 1)
Mco$y_cent=Mco$x
Mco$x_cent=Mco$y-30000
Oco$y_cent=-Oco$x
Oco$x_cent=Oco$y-30000
Oytemp=-Oco$x
Oxtemp=Oco$y-30000
Oco$x_cent=-0.5*Oxtemp+0.866*Oytemp+5000
Oco$y_cent=-0.5*Oytemp-0.866*Oxtemp
coo=rbind(Mco,Oco)
cellchat@images=list()
cellchat@images$coordinates=coo[,3:4]

cellchat <- netAnalysis_computeCentrality(cellchat,slot.name = "netP")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",height = 20,color.heatmap="RdPu")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",height = 20,color.heatmap="GnBu")
#########AMH###################
netVisual_heatmap(cellchat, signaling = c("AMH"), color.heatmap = "Reds")#生成AMH.heatmap.pdf
netVisual_aggregate(cellchat, signaling = "AMH",top=0.03,remove.isolate = FALSE)#生成AMH.network.pdf
netVisual_aggregate(cellchat, signaling = "AMH", layout = "spatial",edge.width.max = 2,top=0.03,
                    vertex.size.max = 1, alpha.image = 0.6, vertex.label.cex =4, remove.isolate = FALSE,point.size = 3)#生成AMH.spatial.pdf
netVisual_bubble(cellchat,signaling = "AMH",
                 sources.use = c("mgran.AF3","cgran.AF3","gran2.AF3","oocyte.AF3","cgran.AF2",
                                 "mgran.AF1","cgran.AF1","gran2.AF1","oocyte.AF1","AF4"),
                 targets.use = c("mgran.AF3","AF4"),sort.by.target=TRUE,angle.x = 45,
                 remove.isolate = FALSE)
plotGeneExpression(cellchat,signaling = "AMH")

#########OPIOID###################
netVisual_heatmap(cellchat, signaling = c("OPIOID"), color.heatmap = "Reds")
netVisual_aggregate(cellchat, signaling = "OPIOID",top=0.1,remove.isolate = FALSE)
netVisual_aggregate(cellchat, signaling = "OPIOID", layout = "spatial",edge.width.max = 2,top=0.1,
                    vertex.size.max = 1, alpha.image = 0.6, vertex.label.cex =4, remove.isolate = FALSE,point.size = 3)
netVisual_bubble(cellchat,signaling = "OPIOID",
                 sources.use = c("gran2.AF1","itheca.AF1","etheca.AF1","gran2.AF2","itheca.AF2","etheca.AF2",
                                 "gran2.AF3","cgran.AF3","mgran.AF3","itheca.AF3","etheca.AF3",
                                 "AF4","cluster9","stroma0"),
                 targets.use = c("oocyte.AF1"),sort.by.target=TRUE,angle.x = 45,
                 remove.isolate = FALSE,thresh=0.05)
plotGeneExpression(cellchat,signaling = "OPIOID")


#########GAS###################
netVisual_heatmap(cellchat, signaling = c("GAS"), color.heatmap = "Reds")
netVisual_aggregate(cellchat, signaling = "GAS",top=0.03,remove.isolate = FALSE)
netVisual_aggregate(cellchat, signaling = "GAS", layout = "spatial",edge.width.max = 2,top=0.03,
                    vertex.size.max = 1, alpha.image = 0.6, vertex.label.cex =4, remove.isolate = FALSE,point.size = 3)
netVisual_bubble(cellchat,signaling = "GAS",
                 sources.use = c("gran2.AF1","cgran.AF1","mgran.AF1","itheca.AF1","etheca.AF1",
                                 "etheca.AF3","AF4","cluster9","stroma0"),
                 targets.use = c("AF4","cluster9","stroma0"),sort.by.target=TRUE,angle.x = 45,
                 remove.isolate = FALSE)
plotGeneExpression(cellchat,signaling = "GAS")
