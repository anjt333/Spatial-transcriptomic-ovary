rm(list=ls())
testAB.integrated=readRDS(allnew.11.rds")
slidelist=c("A","B","C","D","E","F","G","I","J","K","L","M","N","O")
library(stringr)
barc=data.frame()
for(s in 1:14){
  data=paste0(slidelist[s],".csv")
  if(file.exists(data)){
    temp=read.table(data,header=T,sep=",")
    colnames(temp)=c("Barcode","Follicle")
    temp$Follicle=str_to_title(temp$Follicle)
    temp$Barcode=paste0(temp$Barcode,"_",s)
    barc=rbind(barc,temp)
  }
}
unique(barc$Follicle)
barc$Follicle=str_sub(barc$Follicle,1,nchar(barc$Follicle)-2)
subsetf=subset(testAB.integrated,cells=barc$Barcode)
barc$Cluster=Idents(subsetf)
table(barc$Cluster,barc$Follicle)
count=table(barc$Cluster,barc$Follicle)
data=as.data.frame.array(count[,c(5,4,1,3,2)])
data_rownames <- rownames(data)
data_colnames <- colnames(data)
data$cluster <- data_rownames
library(reshape2)
data_m <- melt(data, id.vars=c("cluster"))
data_m=data_m[data_m$variable!="Health"]
data_m$cluster <- factor(data_m$cluster, levels=c("gran","itheca","etheca",
                                                  "stroma1","stroma0","stroma3","oocyte",
                                                  "muscle","endo","immune","injured"))
data_m
library(ggplot2)
mycolor=c("#C64F4F","#009933","#98DF8AFF","#17BECFFF","#AEC7E8FF","#FF9896FF",
          "#C0C0C0","#CD853F","#F7B6D2FF","#FFD700","#DB7093")
ggplot(data_m, aes(x=value, y=variable)) +
  geom_bar(stat="identity", position="fill", aes(fill=cluster)) +
  scale_x_continuous(labels = scales::percent)+
  scale_fill_manual(values=mycolor)+
  theme_test()+xlab("Percent")+ylab("Slide")

library(dplyr)
library(ggplot2)
barc$Follicle=factor(barc$Follicle,levels=c("Healthy","Atretic Antral","Atretic White","Atretic Pink"))
Idents(subsetf)

x="endo"
sub=subset(subsetf,cells=barc[barc$Cluster==x,1])
DefaultAssay(subsetf)="SCT"
sub=ScaleData(sub, features = rownames(sub))
Idents(sub)=barc[barc$Cluster==x,2]
levels(Idents(sub))
datax=cbind(data.frame(apply(subset(sub,idents="Healthy")@assays$SCT@data,1,median)),
            data.frame(apply(subset(sub,idents="Atretic Antral")@assays$SCT@data,1,median)),
            data.frame(apply(subset(sub,idents="Atretic White")@assays$SCT@data,1,median)),
            data.frame(apply(subset(sub,idents="Atretic Pink")@assays$SCT@data,1,median)))
colnames(datax)=levels(Idents(sub))
o=datax[c("ATG3","NBAS","ENY2"),]
o$gene=factor(rownames(o),levels=rownames(o))
o=melt(o,measure.vars = levels(Idents(sub)),variable.name = "Stage")
o$value=as.numeric(o$value)
ggplot(data = o,aes(x=Stage,y=value,group=gene))+geom_line()+
  ylab("Expression")+theme_classic()
