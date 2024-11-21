testAB.integrated=readRDS("allnew.11.rds")
stroma=subset(testAB.integrated,idents=c("stroma1","stroma0","stroma3"))
stroma <- ScaleData(stroma, features = rownames(stroma))
stroma <- RunPCA(stroma,verbose = FALSE,features = VariableFeatures(object = stroma))
stroma <- FindNeighbors(stroma, dims = 1:20)
DimPlot(stroma,cols=c("#17BECFFF","#AEC7E8FF","#FF9896FF"),label=T,label.size = 5,pt.size=1)

slidelist=c("A_2020001_2","B_2020001_2","C_S20006_2","D_S20006_2","E_S20006_2","F_2020001_4","G_2020001_5",
            "I_S21004_10","J_S21006_1","K_S21006_1","L_S21004_13","M_S21004_1","N_S21004_1","O_S21004_1")
for(i in slidelist){
  pdf(paste0(i,".stroma.pdf"),height=6,width=8)
  p=SpatialDimPlot(stroma,images=i,cols=c("#17BECFFF","#AEC7E8FF","#FF9896FF"),crop=F,pt.size.factor=1)+
    scale_fill_manual(values=c("#17BECFFF","#AEC7E8FF","#FF9896FF"))+ggtitle(i)+theme(plot.title = element_text(hjust = 0.5))
  print(p)
  dev.off()
}

#########stroma.markers##################
pbmc.markers <- FindAllMarkers(stroma, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
stroma1.marker=pbmc.markers[pbmc.markers$cluster=="stroma1",]
stroma1.marker=stroma1.marker[,c(7,1:5)]
stroma0.marker=pbmc.markers[pbmc.markers$cluster=="stroma0",]
stroma0.marker=stroma0.marker[,c(7,1:5)]
stroma3.marker=pbmc.markers[pbmc.markers$cluster=="stroma3",]
stroma3.marker=stroma3.marker[,c(7,1:5)]

DoHeatmap(stroma, features = c(stroma1.marker$gene,stroma3.marker$gene),group.colors=c("#17BECFFF","#AEC7E8FF","#FF9896FF"))+
  scale_fill_gradientn(colors = c("blue","white","red"))

stroma0.marker <- FindMarkers(stroma, min.pct = 0.25, logfc.threshold = 0.25,ident.1 = "stroma0")
stroma0.marker2 <- FindMarkers(stroma, min.pct = 0.25, logfc.threshold = 0.05,ident.1 = "stroma0")
stroma0.marker2$gene=rownames(stroma0.marker2)
stroma0.marker2$change = ifelse(stroma0.marker2$p_val_adj < 0.05 & abs(stroma0.marker2$avg_log2FC) >=0.25, 
                        ifelse(stroma0.marker2$avg_log2FC>=0.25 ,'Up','Down'),
                        'Stable')
ggplot(stroma0.marker2,aes(avg_log2FC, -log10(p_val_adj),colour=change))+geom_point()+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")+
  geom_vline(xintercept = c(-0.25,0.25), linetype = "dashed", color = "black")+
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  scale_size_continuous(range = c(0,1))+
  theme_test()#生成stroma0.volcano.pdf

#########GO enrichment (by David database)################
GO=read.table("stroma3.selected.GO.david.txt",header=T,sep="\t")
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
  scale_color_gradient(low = "#9370DB", high = "#FF00FF")+ggtitle("Enriched ontology terms of stroma3 feature genes")+
  theme_test()#生成stroma3.GO.pdf


SpatialFeaturePlot(stroma, images=i,features = "WISP2")
slidelist=c("A_2020001_2","B_2020001_2","C_S20006_2","D_S20006_2","E_S20006_2","F_2020001_4","G_2020001_5",
            "I_S21004_10","J_S21006_1","K_S21006_1","L_S21004_13","M_S21004_1","N_S21004_1","O_S21004_1")

for(i in slidelist){
  pdf(paste0(i,".WFIKKN2.pdf"),height=6,width=6)
  p=SpatialFeaturePlot(stroma, images=i,features = "WFIKKN2",crop=F,pt.size.factor=1.2)+
    scale_fill_gradientn(colors = rev(brewer.pal(11,"Spectral")),limits=c(-1,2))+
    ggtitle(i)+theme(plot.title = element_text(hjust = 0.5))
  print(p)
  dev.off()
}
for(i in slidelist){
  pdf(paste0(i,".WISP2.pdf"),height=6,width=6)
  p=SpatialFeaturePlot(stroma, images=i,features = "WISP2",crop=F,pt.size.factor=1.2)+
    scale_fill_gradientn(colors = rev(brewer.pal(11,"Spectral")),limits=c(-1,2.5))+
    ggtitle(i)+theme(plot.title = element_text(hjust = 0.5))
  print(p)
  dev.off()
}
SpatialFeaturePlot(stroma, features = "WFIKKN2",crop=F,pt.size.factor=1)
