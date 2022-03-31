ensg2gene <- read.table("./ensg2gene")
gene_biotype <-read.table("./Homo_sapiens.GRCh38.95.genebiotype",
                          sep = "\t",
                          as.is = T,
                          col.names = c("gene", "biotype")
)
pcg <- gene_biotype[which(gene_biotype$biotype=="protein_coding"),]
pcg <- merge(pcg,ensg2gene,by.x="gene",by.y="V1")

load('./GSE98370_gene_expression.Rdata')
#expres<-read.table("./GSE98370_normalisedcounts.txt",header = T,sep ="\t")

sample<-read.table("./ACAT/sample_info.txt",sep = "\t",header = F)
sample_info<-data.frame("n"=colnames(expres),"s"=sample$V7,"t"=sample$V10)
sample_info$t[c(30,31,32)]<-"Ad"
################## Integrate sample at the same time points (mean) #######################
colnames(fpkm)<-sample_info$t
final11<-c()
for(i in unique(sample_info$t)){
  if(length(which(colnames(fpkm)%in%i))==1){
    final11<-cbind(final11,fpkm[,which(colnames(fpkm)%in%i)])
  }else{
    final11<-cbind(final11,apply(fpkm[,which(colnames(fpkm)%in%i)], 1, mean))
  }}
dim(final11)
colnames(final11)<-unique(sample_info$t)
boxplot(final11)
############### Exclude low expressed genes ######################
lowfpkm<-c()
for (i in 1:length(final11[,1])) {
  if(max(final11[i,])<1){
    lowfpkm<-c(lowfpkm,i)
  }
}
expr<-final11[-lowfpkm,]  
boxplot(expr)
dim(expr)

###################### Normalization ################
final2<-log2(expr+1)
boxplot(final2)
fpkm_pro<-final2[which(rownames(final2)%in%pcg$gene),]
final2<-scale(fpkm_pro,center = T,scale = F)
dim(final2)
boxplot(final2)

#################### Selecting genes with significant expression variation during retinal development #################
TT<-c()
TT2<-c()
EE<-c()
PP<-c()
FP<-c()
x<-1:16
for (j in 1:nrow(final2)) {
  t_sta<-lm(formula = final2[j,] ~ x)
  T_value<-summary(t_sta)$coef[2,3]
  F_value<-anova(t_sta)$`F value`[1]
  TT<-c(TT,T_value)
  TT2<-c(TT2,F_value)
  EE<-c(EE,summary(t_sta)$coef[2,1])
  PP<-c(PP,summary(t_sta)$coef[2,4])
}
TT1<-cbind(TT,EE,PP)
rownames(TT1)<-rownames(final2)
TT1<-as.data.frame(TT1)
TT1$fdr<-p.adjust(TT1$PP,method = "fdr")
siggene<-TT1[which(TT1$fdr<0.05),]

################## Save ##############
save(file="**/development.retina.RData", final2,Tsta,siggene)

#heatmap
library(reshape2)
library("clusterProfiler")
library("S4Vectors")
library("org.Hs.eg.db")
gene<-read.csv('**.txt',sep = "\t",header = F)
g1<-as.character(unique(gene$V1))
tran2<-bitr(g1, fromType="SYMBOL", toType=c("ENSEMBL"), OrgDb="org.Hs.eg.db")

######################################
load("./development.retina.RData")
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


#save(file="./development.retina.RData", fgene, final2,Tsta,
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

