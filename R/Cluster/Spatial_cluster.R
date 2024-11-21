testAB.integrated=readRDS("allnew.11.rds")
mycolor=c("#C64F4F","#009933","#98DF8AFF","#17BECFFF","#AEC7E8FF","#FF9896FF",
          "#C0C0C0","#CD853F","#F7B6D2FF","#FFD700","#DB7093")
DimPlot(testAB.integrated,cols=mycolor,label=T,label.size = 5,pt.size=0.8)

slidelist=c("A_2020001_2","B_2020001_2","C_S20006_2","D_S20006_2","E_S20006_2","F_2020001_4","G_2020001_5",
            "I_S21004_10","J_S21006_1","K_S21006_1","L_S21004_13","M_S21004_1","N_S21004_1","O_S21004_1")
SpatialDimPlot(subset(testAB.integrated,idents="stroma0"),images="M_S21004_1",cols=mycolor[5],crop=F)+
    scale_fill_manual(values=mycolor[5])+theme(plot.title = element_text(hjust = 0.5))
SpatialDimPlot(subset(testAB.integrated,idents="injured"),images="E_S20006_2",cols=mycolor[11],crop=F,pt.size=1.2)+
  scale_fill_manual(values=mycolor[11])+theme(plot.title = element_text(hjust = 0.5))
for(i in slidelist){
  pdf(paste0(i,".cluster.pdf"),height=6,width=8)
  p=SpatialDimPlot(testAB.integrated,images=i,cols=mycolor,crop=F,pt.size=1.2)+
    scale_fill_manual(values=mycolor)+ggtitle(i)+theme(plot.title = element_text(hjust = 0.5))
  print(p)
  dev.off()
}
for(i in slidelist){
  pdf(paste0(i,".HE.pdf"),height=6,width=8)
  p=SpatialDimPlot(testAB.integrated,images=i,cols=mycolor,alpha = 0,crop=F)+
    scale_fill_manual(values=mycolor)+ggtitle(i)+theme(plot.title = element_text(hjust = 0.5))
  print(p)
  dev.off()
}
gene=c("AMH","STAR","CYP17A1","COL1A1","DCN","ZP3","TAGLN","MYH11","CCL21")
for(i in gene){
  pdf(paste0(i,".umap.pdf"),height=5,width=5)
  p=FeaturePlot(testAB.integrated,features = i,cols=c("grey","lightgrey","red"))
  print(p)
  dev.off()
}

split_b<-str_split(rownames(testAB.integrated@meta.data),"_")
Idents(testAB.integrated)<-sapply(split_b,"[",2)
My_levels <- c("slideA","slideB","slideC","slideD","slideE","slideF","slideG","slideI",
               "slideJ","slideK","slideL","slideM","slideN","slideO")
names(My_levels) <- levels(Idents(testAB.integrated))
testAB.integrated <- RenameIdents(testAB.integrated, My_levels)
DimPlot(testAB.integrated)

VlnPlot(testAB.integrated, features = "nFeature_Spatial", pt.size = 0) + NoLegend()+ggtitle("Number of detected genes")+theme(plot.title = element_text(hjust = 0.5))+
  geom_boxplot(width=.2,col="black",fill="white")

mean(testAB.integrated[["nFeature_Spatial"]]$nFeature_Spatial,na.rm=T)
sum(colSums(subset(testAB.integrated,idents="slideO")@assays$Spatial[])>=1)

##################celltype propotion###################
split_b<-str_split(rownames(testAB.integrated@meta.data),"_")
testAB.integrated$sampleID<-sapply(split_b,"[",2)
table(testAB.integrated$sampleID)
table(Idents(testAB.integrated),testAB.integrated$sampleID)
count=table(Idents(testAB.integrated),testAB.integrated$sampleID)
count=as.data.frame.array(count[,c(1,7:14,2:6)])
colnames(count)=c("A","B","C","D","E","F","G","I","J","K","L","M","N","O")
data=data.frame("G"=count$G,"L"=count$L,"I"=count$I,"A"=count$A,"B"=count$B,
                "C"=count$C,"D"=count$D,"J"=count$J,"K"=count$K,
                "F"=count$F,"M"=count$M,"N"=count$N,"O"=count$O,"E"=count$E)
rownames(data)=rownames(count)
data_rownames <- rownames(data)
data_colnames <- colnames(data)
data$cluster <- data_rownames
data_m <- melt(data, id.vars=c("cluster"))
data_m$cluster <- factor(data_m$cluster, levels=c("gran","itheca","etheca",
                                                  "stroma1","stroma0","stroma3","oocyte",
                                                  "muscle","endo","immune","injured"))
data_m

ggplot(data_m, aes(x=value, y=variable)) +
  geom_bar(stat="identity", position="fill", aes(fill=cluster)) +
  scale_x_continuous(labels = scales::percent)+
  scale_fill_manual(values=mycolor)+
  theme_test()+xlab("Percent")+ylab("Slide")

##############marker###################
pbmc.markers <- FindAllMarkers(testAB.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
DoHeatmap(testAB.integrated, features = pbmc.markers$gene,group.colors=mycolor)+
  scale_fill_gradientn(colors = c("blue","white","red"))

#################detected gene######################
VlnPlot(testAB.integrated, features = "nFeature_Spatial", pt.size = 0,cols=mycolor) + NoLegend()+
  scale_fill_manual(values=mycolor)+ggtitle("Number of detected genes")+theme(plot.title = element_text(hjust = 0.5))+
  geom_boxplot(width=.2,col="black",fill="white")

data=data.frame(testAB.integrated[["nFeature_Spatial"]])
data$nFeature_Spatial
ggplot(data = data,aes(x=nFeature_Spatial))+geom_histogram(bins = 30,fill="#6495ED", color="#e9ecef")+
  theme_test()+labs(x="Number of detected genes",y="Number of spots")

slidelist=c("A_2020001_2","B_2020001_2","C_S20006_2","D_S20006_2","E_S20006_2","F_2020001_4","G_2020001_5",
            "I_S21004_10","J_S21006_1","K_S21006_1","L_S21004_13","M_S21004_1","N_S21004_1","O_S21004_1")
for(i in slidelist){
  pdf(paste0(i,".detected.pdf"),height=6,width=6)
  p=SpatialFeaturePlot(testAB.integrated,images=i,features = "nFeature_Spatial",image.alpha = 0,crop=F,pt.size=1.2)+
    scale_fill_gradientn(colors = c("#7EEBA6","#FAFF81","#E64229","#970F0F"),limits=c(0,10000))+
    ggtitle(i)+theme(plot.title = element_text(hjust = 0.5))
  print(p)
  dev.off()
}
