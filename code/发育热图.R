library(reshape2)
library("clusterProfiler")
library("S4Vectors")
library("org.Hs.eg.db")
gene<-read.csv('**.txt',sep = "\t",header = F)
g1<-as.character(unique(gene$V1))
tran2<-bitr(g1, fromType="SYMBOL", toType=c("ENSEMBL"), OrgDb="org.Hs.eg.db")

######################################
load("E:/课题/ZB_eyedisease/发育/development.retina.RData")
library("limma")
library("ggplot2")
library(ggalt)
library(reshape2)
library(pheatmap)
#############################################################

fgene<-final2[which(rownames(final2)%in%tran2$ENSEMBL),]
fgene1<-fgene
fgene1<-((2-(-2))*((fgene1-min(fgene1))/(max(fgene1)-min(fgene1))))+(-2)



fple_up<-fgene1[which(rownames(fgene1)%in%rownames(Tsta)[which(Tsta$EE>0)]),]
fple_do<-fgene1[which(rownames(fgene1)%in%rownames(Tsta)[which(Tsta$EE<0)]),]




#save(file="E:/课题/ZB_eyedisease/发育/development.retina.RData", fgene, final2,Tsta,
# siggene,nonple,fple,heaemap_ple,heaemap_nonple)
fple<-as.data.frame(fgene1)

allgene<-rbind(fple)
allgene$direct<-c()
for (i in 1:nrow(allgene)) {
  allgene$direct[i]<-Tsta$EE[which(rownames(Tsta)%in%rownames(allgene)[i])]
}

for (i in 1:nrow(allgene)) {
  if(length(which(tran2$ENSEMBL%in%rownames(allgene)[i]))>0){
    rownames(allgene)[i]<-as.character(tran2$SYMBOL[which(tran2$ENSEMBL%in%rownames(allgene)[i])])
  }else{
    rownames(allgene)[i]<-as.character(tran2$SYMBOL[which(tran2$ENSEMBL%in%rownames(allgene)[i])])
  }
}
allgene<-allgene[order(allgene$direct),]

heaemap_allgene<-pheatmap(allgene,color = colorRampPalette(colors = c("#3952A3","white","#E7161B"))(100),
                          cluster_rows = F,cluster_cols = F,clustering_method = 'average')


