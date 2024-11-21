gran=readRDS("gran.rds")
gran <- ScaleData(gran, features = rownames(gran))
gran <- RunPCA(gran,verbose = FALSE,features = VariableFeatures(object = gran))
gran <- FindNeighbors(gran, dims = 1:20)
gran <- RunUMAP(gran, dims = 1:20)
DimPlot(gran,label = F,pt.size=2,cols=c("#e0c45c","#85203b","#cb3b3b"))#生成gran.umap.pdf

slidelist=c("A_2020001_2","B_2020001_2","C_S20006_2","D_S20006_2","E_S20006_2","F_2020001_4","G_2020001_5",
            "I_S21004_10","J_S21006_1","K_S21006_1","L_S21004_13","M_S21004_1","N_S21004_1","O_S21004_1")

for(i in slidelist){
  pdf(paste0(i,".gran.pdf"),height=6,width=8)
  p=SpatialDimPlot(gran,images=i,cols=c("#e0c45c","#85203b","#cb3b3b"),crop=F,pt.size.factor=1)+
    scale_fill_manual(values=c("#e0c45c","#85203b","#cb3b3b"))+ggtitle(i)+theme(plot.title = element_text(hjust = 0.5))
  print(p)
  dev.off()
}
#########gran.markers##################
pbmc.markers <- FindAllMarkers(gran, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
cgran.marker=pbmc.markers[pbmc.markers$cluster=="cgran",]
cgran.marker=cgran.marker[,c(7,1:5)]
gran2.marker=pbmc.markers[pbmc.markers$cluster=="gran2",]
gran2.marker=gran2.marker[,c(7,1:5)]
mgran.marker=pbmc.markers[pbmc.markers$cluster=="mgran",]
mgran.marker=mgran.marker[,c(7,1:5)]
DoHeatmap(gran, features = pbmc.markers$gene,group.colors=c("#e0c45c","#85203b","#cb3b3b"))+
  scale_fill_gradientn(colors = c("blue","white","red"))


gene="GSTA1"
VlnPlot(gran, features = gene, pt.size = 0,cols=c("#e0c45c","#85203b","#cb3b3b")) + NoLegend()+
  scale_fill_manual(values=c("#e0c45c","#85203b","#cb3b3b"))+ggtitle(gene)+theme(plot.title = element_text(hjust = 0.5))+
  geom_boxplot(width=.2,col="black",fill="white")

for(i in slidelist){
  pdf(paste0(i,".GSTA1.pdf"),height=6,width=6)
  p=SpatialFeaturePlot(gran, images=i,features = "GSTA1",crop=F,pt.size.factor=1.2)+
    scale_fill_gradientn(colors = rev(brewer.pal(11,"Spectral")),limits=c(1,4))+
    ggtitle(i)+theme(plot.title = element_text(hjust = 0.5))
  print(p)
  dev.off()
}
for(i in slidelist){
  pdf(paste0(i,".AMH.pdf"),height=6,width=6)
  p=SpatialFeaturePlot(gran, images=i,features = "AMH",crop=F,pt.size.factor=1.2)+
    scale_fill_gradientn(colors = rev(brewer.pal(11,"Spectral")),limits=c(1,3.5))+
    ggtitle(i)+theme(plot.title = element_text(hjust = 0.5))
  print(p)
  dev.off()
}
for(i in slidelist){
  pdf(paste0(i,".COL1A2.pdf"),height=6,width=6)
  p=SpatialFeaturePlot(gran, images=i,features = "COL1A2",crop=F,pt.size.factor=1.2)+
    scale_fill_gradientn(colors = rev(brewer.pal(11,"Spectral")),limits=c(0,3))+
    ggtitle(i)+theme(plot.title = element_text(hjust = 0.5))
  print(p)
  dev.off()
}
###########Cyto###############
library(CytoTRACE)
raw.data <- as.matrix(gran@assays$Spatial@counts)
x=as.character(Idents(gran))
names(x)=Cells(gran)
c = sweep(raw.data, 2, raw.data["GAPDH",], `/`)
rm(list=ls()[! ls() %in% c("c","x")])
setwd("C:/Users/Dell/Desktop")
results <- CytoTRACE(c)
plotCytoTRACE(results, phenotype = x)
gran[["cyto"]]=1-results$CytoTRACE
VlnPlot(gran, features = "cyto", pt.size = 0,cols=c("#e0c45c","#85203b","#cb3b3b")) + NoLegend()+
  scale_fill_manual(values=c("#e0c45c","#85203b","#cb3b3b"))+ggtitle("cyto")+theme(plot.title = element_text(hjust = 0.5))+
  geom_boxplot(width=.2,col="black",fill="white")
SpatialFeaturePlot(gran,features="cyto")
###########Stem marker##########
PGC=read.table("E:/spatial transcriptome/PGC/merge.female.pgc.somatic.txt",header=T,row.names = 1)
cluster=read.table("E:/spatial transcriptome/PGC/final.cluster.txt",header=T,row.names = NULL)
soma.cluster= cluster[which(substr(cluster$Cluster,1,11) == "Female_Soma"),]
soma.cluster= soma.cluster[which(soma.cluster$Cell%in%colnames(PGC)),]
PGC.soma=PGC[,soma.cluster$Cell]
PGC.soma.gene=PGC.soma[which(rowSums(PGC.soma>1)>=2), ]
Pluripotency.markers=read.table("E:/spatial transcriptome/PGC/Pluripotency markers.2013.txt",sep="\t")
Pluripotency.markers=Pluripotency.markers[seq(5,401,3),]
Pluripotency.markers=unique(Pluripotency.markers)
Pluripotency.genes=read.table("E:/spatial transcriptome/PGC/Pluripotency genes.2015.txt",sep="\t",header=T)
Pluripotency.genes=unique(Pluripotency.genes$Pluripotency.genes)
Pluripotency.markers=Pluripotency.markers[which(Pluripotency.markers %in% rownames(PGC.soma.gene))]
Pluripotency.genes=Pluripotency.genes[which(Pluripotency.genes %in% rownames(PGC.soma.gene))]
Pluripotencylist=c(Pluripotency.markers,Pluripotency.genes)
Pluripotencylist=unique(Pluripotencylist)

raw.data <- as.matrix(gran@assays$Spatial@counts)
raw.data <- raw.data[which(rownames(raw.data) %in% c(Pluripotencylist,"GAPDH")),]
mdata = sweep(raw.data, 2, raw.data["GAPDH",], `/`)
mdata = mdata[-which(rownames(mdata)=="GAPDH"),]
median=rep(0,ncol(mdata))
for (i in 1:ncol(mdata)){
  median[i]=median(mdata[,i])
}
gran[["stem"]]=median
VlnPlot(gran, features = "stem",cols=c("#e0c45c","#85203b","#cb3b3b")) + NoLegend()+
  scale_fill_manual(values=c("#e0c45c","#85203b","#cb3b3b"))+ggtitle("stem")+theme(plot.title = element_text(hjust = 0.5))
