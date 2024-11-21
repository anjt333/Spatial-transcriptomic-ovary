##################膜细胞部分###########################
theca=readRDS("theca.rds")
theca <- ScaleData(theca, features = rownames(theca))
theca <- RunPCA(theca,verbose = FALSE,features = VariableFeatures(object = theca))
theca <- FindNeighbors(theca, dims = 1:20)
DimPlot(theca,cols=c("#98DF8AFF","#009933"),label=T,label.size = 5,pt.size=1)

slidelist=c("A_2020001_2","B_2020001_2","C_S20006_2","D_S20006_2","E_S20006_2","F_2020001_4","G_2020001_5",
            "I_S21004_10","J_S21006_1","K_S21006_1","L_S21004_13","M_S21004_1","N_S21004_1","O_S21004_1")
for(i in slidelist){
  pdf(paste0(i,".theca.pdf"),height=6,width=8)
  p=SpatialDimPlot(theca,images=i,cols=c("#98DF8AFF","#009933"),crop=F,pt.size.factor=1)+
    scale_fill_manual(values=c("#98DF8AFF","#009933"))+ggtitle(i)+theme(plot.title = element_text(hjust = 0.5))
  print(p)
  dev.off()
}
#########theca.markers##################
pbmc.markers <- FindAllMarkers(theca, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
DoHeatmap(theca, features = pbmc.markers$gene,group.colors=c("#98DF8AFF","#009933"))+
  scale_fill_gradientn(colors = c("blue","white","red"))
itheca.marker=pbmc.markers[pbmc.markers$cluster=="itheca",]
itheca.marker=itheca.marker[,c(7,1:5)]
etheca.marker=pbmc.markers[pbmc.markers$cluster=="etheca",]
etheca.marker=etheca.marker[,c(7,1:5)]

#########GO plot (by David database)################
GO=read.table("etheca.selected.GO.david.txt",header=T,sep="\t")
GO$Genes<-gsub(", ","/",GO$Genes)
split_Term<-str_split(GO$Term,"~")
GO$Term1<-sapply(split_Term,"[",1)
GO$Term2<-sapply(split_Term,"[",2)
GO$GeneRatio=as.numeric(GO$Count)/as.numeric(GO$List.Total)
GO$BGRatio <- as.numeric(GO$Pop.Hits)/as.numeric(GO$Pop.Total)
GO.dot=data.frame(cbind(Term=GO$Term1,Description=GO$Term2,GeneRatio=GO$GeneRatio,BGRatio=GO$BGRatio,
                        pvalue=GO$PValue,FDR=GO$FDR,GeneID=GO$Genes,Count=GO$Count))
GO.dot$pvalue=as.numeric(GO.dot$pvalue)
GO.dot$FDR=as.numeric(GO.dot$FDR)
GO.dot$Count=as.numeric(GO.dot$Count)
GO.dot$GeneRatio=round(as.numeric(GO.dot$GeneRatio),2)
GO.dot <- GO.dot[order(GO.dot$GeneRatio,-GO.dot$pvalue),]
GO.dot$Description <- factor(GO.dot$Description,levels=GO.dot$Description)
ggplot(GO.dot,aes(x = GeneRatio,y = Description))+
  geom_point(aes(color = pvalue,
                 size = Count))+
  scale_color_gradient(low = "#9370DB", high = "#FF00FF")+ggtitle("Enriched ontology terms of etheca feature genes")+
  theme_test()

##############markergene#######################
VlnPlot(theca, features = "ACTA2", pt.size = 0,cols=c("#98DF8AFF","#009933")) + NoLegend()+
  scale_fill_manual(values=c("#98DF8AFF","#009933"))+ggtitle("ACTA2")+theme(plot.title = element_text(hjust = 0.5))+
  geom_boxplot(width=.2,col="black",fill="white")#生成ACTA2.vlnplot.pdf
SpatialFeaturePlot(theca, features = "SCARB1",crop=F,pt.size.factor=1)
SpatialFeaturePlot(theca, features = "ACTA2",crop=F,pt.size.factor=1)

for(i in slidelist){
  pdf(paste0(i,".ACTA2.pdf"),height=6,width=6)
  p=SpatialFeaturePlot(theca, images=i,features = "ACTA2",crop=F,pt.size.factor=1.2)+
    scale_fill_gradientn(colors = rev(brewer.pal(11,"Spectral")),limits=c(-1,4))+
    ggtitle(i)+theme(plot.title = element_text(hjust = 0.5))
  print(p)
  dev.off()
}
