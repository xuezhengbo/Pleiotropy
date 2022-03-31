library(ggplot2)
library(dplyr)
library(forcats)
library(cowplot)
library(stringr)
ensg2gene <- read.table("/share/pub/chenfk/Data/eye_database/multi/ensg2gene")
gene_biotype <-
  read.table(
    "/share/pub/chenfk/Data/eye_database/multi/Homo_sapiens.GRCh38.95.genebiotype",
    sep = "\t",
    as.is = T,
    col.names = c("gene", "biotype")
  )
pcg <- gene_biotype[which(gene_biotype$biotype=="protein_coding"),]
pcg <- merge(pcg,ensg2gene,by.x="gene",by.y="V1")
#################��ѡ����######################
gene<-read.csv("/share/pub/chenfk/project/ZB_eyedisease/gene_list/gene_list_new.txt",sep = "\t")
g1<-as.character(unique(gene$Peiotropy))
g1<-g1[which(g1%in%pcg$V2)]
g2<-as.character(unique(gene$Non.peiotropy))
g2<-g2[-which(g2%in%"")]
g2<-g2[which(g2%in%pcg$V2)]
############################################
expregene_organoid<-read.table("/share/pub/chenfk/project/ZB_eyedisease/single/GSE104827/cluster_data/organoid/big_group/expres0.25.txt",sep = " ",header = F)
expregene_organoid<-expregene_organoid[which(expregene_organoid$V2%in%pcg$V2),]
length(unique(expregene_organoid$V2))
ple_gene<-g1[which(g1%in%expregene_organoid$V2)]
nonple_gene<-g2[which(g2%in%expregene_organoid$V2)]

write.table(ple_gene,"/share/pub/chenfk/project/ZB_eyedisease/single/GSE104827/cluster_data/fig7/ple_gene.txt",sep = "\t",col.names = F,row.names = F,quote = F)
write.table(nonple_gene,"/share/pub/chenfk/project/ZB_eyedisease/single/GSE104827/cluster_data/fig7/nonple_gene.txt",sep = "\t",col.names = F,row.names = F,quote = F)



org<-read.table("/share/pub/chenfk/project/ZB_eyedisease/single/GSE104827/cluster_data/organoid/exprMatrix_singlegene.tsv",sep = " ",header = T,quote = "")
org_expre<-org[which(rownames(org)%in%c(g1,g2)),]
write.table(org_expre,"/share/pub/chenfk/project/ZB_eyedisease/single/GSE104827/cluster_data/fig7/organoid_candidate.txt",sep = "\t",col.names = F,row.names = F,quote = F)

per<-read.table("/share/pub/chenfk/project/ZB_eyedisease/single/GSE104827/cluster_data/peripheral/exprMatrix_singlegene.tsv",sep = " ",header = T,quote = "")
per_expre<-org[which(rownames(per)%in%c(g1,g2)),]
write.table(per_expre,"/share/pub/chenfk/project/ZB_eyedisease/single/GSE104827/cluster_data/fig7/peripheral_candidate.txt",sep = "\t",col.names = F,row.names = F,quote = F)

fove<-read.table("/share/pub/chenfk/project/ZB_eyedisease/single/GSE104827/cluster_data/foveal/exprMatrix_singlegene.tsv",sep = " ",header = T,quote = "")
fove_expre<-org[which(rownames(fove)%in%c(g1,g2)),]
write.table(fove_expre,"/share/pub/chenfk/project/ZB_eyedisease/single/GSE104827/cluster_data/fig7/foveal_candidate.txt",sep = "\t",col.names = F,row.names = F,quote = F)





##############################################################


org<-read.table("/share/pub/chenfk/project/ZB_eyedisease/single/GSE104827/cluster_data/fig7/organoid_candidate.txt",sep = "\t",header = T,quote = "")
org<-org[which(org[,1]%in%c(ple_gene,nonple_gene)),]
colnames(org)<-gsub("X","",colnames(org))
per<-read.table("/share/pub/chenfk/project/ZB_eyedisease/single/GSE104827/cluster_data/fig7/peripheral_candidate.txt",sep = "\t",header = T,quote = "")
per<-per[which(per[,1]%in%c(ple_gene,nonple_gene)),]
colnames(per)<-gsub("X","",colnames(per))
fove<-read.table("/share/pub/chenfk/project/ZB_eyedisease/single/GSE104827/cluster_data/fig7/foveal_candidate.txt",sep = "\t",header = T,quote = "")
fove<-fove[which(fove[,1]%in%c(ple_gene,nonple_gene)),]
colnames(fove)<-gsub("X","",colnames(fove))
rownames(org)<-org$gene
rownames(per)<-per$gene
rownames(fove)<-fove$gene
fove<-fove[,-1]
org<-org[,-1]
per<-per[,-1]
peripheral_meta <- read.table("/share/pub/chenfk/project/ZB_eyedisease/single/GSE104827/cluster_data/peripheral/meta.tsv",sep = "\t",header = T)
organoid_meta <- read.table("/share/pub/chenfk/project/ZB_eyedisease/single/GSE104827/cluster_data/organoid/meta.tsv",sep = "\t",header = T)
foveal_meta <- read.table("/share/pub/chenfk/project/ZB_eyedisease/single/GSE104827/cluster_data/foveal/meta.tsv",sep = "\t",header = T)
peripheral_meta$cell_type<-as.character(peripheral_meta$cell_type)
organoid_meta$cell_type<-as.character(organoid_meta$cell_type)
foveal_meta$cell_type<-as.character(foveal_meta$cell_type)
foveal<-foveal_meta[which(foveal_meta$cell_type_group%in%"ganglion"),]
#peripheral_umap<-read.table("/share/pub/chenfk/project/ZB_eyedisease/single/GSE104827/cluster_data/peripheral/UMAP.coords.tsv")
for (i in 1:nrow(peripheral_meta)) {
  if(str_detect(peripheral_meta$cell_type[i],"cone")){
    peripheral_meta$cell_type[i]<-"C"
  }else if(str_detect(peripheral_meta$cell_type[i],"rod")){
    peripheral_meta$cell_type[i]<-"R"
  }else if(str_detect(peripheral_meta$cell_type[i],"TCell")){
    peripheral_meta$cell_type[i]<-"TCE"
  }else if(str_detect(peripheral_meta$cell_type[i],"BCell")){
    peripheral_meta$cell_type[i]<-"BCE"
  }else if(str_detect(peripheral_meta$cell_type[i],"Ast")){
    peripheral_meta$cell_type[i]<-"AST"
  }else if(str_detect(peripheral_meta$cell_type[i],"_")){
    t<-unlist(strsplit(as.character(peripheral_meta$cell_type[i]),"_"))
    if(str_detect(t[1],"hBC")){
      peripheral_meta$cell_type[i]<-"HBC"
    }else if(str_detect(t[1],"dBC")){
      peripheral_meta$cell_type[i]<-"DBC"
    }else{
      peripheral_meta$cell_type[i]<-t[1]
    }
  }else{
    
  }
}


for (i in 1:nrow(organoid_meta)) {
  if(str_detect(organoid_meta$cell_type[i],"cone")){
    organoid_meta$cell_type[i]<-"C"
  }else if(str_detect(organoid_meta$cell_type[i],"rod")){
    organoid_meta$cell_type[i]<-"R"
  }else if(str_detect(organoid_meta$cell_type[i],"TCell")){
    organoid_meta$cell_type[i]<-"TCE"
  }else if(str_detect(organoid_meta$cell_type[i],"BCell")){
    organoid_meta$cell_type[i]<-"BCE"
  }else if(str_detect(organoid_meta$cell_type[i],"Ast")){
    organoid_meta$cell_type[i]<-"AST"
  }else if(str_detect(organoid_meta$cell_type[i],"_")){
    t<-unlist(strsplit(as.character(organoid_meta$cell_type[i]),"_"))
    if(str_detect(t[1],"hBC")){
      organoid_meta$cell_type[i]<-"HBC"
    }else if(str_detect(t[1],"dBC")){
      organoid_meta$cell_type[i]<-"DBC"
    }else{
      organoid_meta$cell_type[i]<-t[1]
    }
  }else{
    
  }
}

for (i in 1:nrow(foveal_meta)) {
  if(str_detect(foveal_meta$cell_type[i],"cone")){
    foveal_meta$cell_type[i]<-"C"
  }else if(str_detect(foveal_meta$cell_type[i],"rod")){
    foveal_meta$cell_type[i]<-"R"
  }else if(str_detect(foveal_meta$cell_type[i],"TCell")){
    foveal_meta$cell_type[i]<-"TCE"
  }else if(str_detect(foveal_meta$cell_type[i],"BCell")){
    foveal_meta$cell_type[i]<-"BCE"
  }else if(str_detect(foveal_meta$cell_type[i],"Ast")){
    foveal_meta$cell_type[i]<-"AST"
  }else if(str_detect(foveal_meta$cell_type[i],"_")){
    t<-unlist(strsplit(as.character(foveal_meta$cell_type[i]),"_"))
    if(str_detect(t[1],"hBC")){
      foveal_meta$cell_type[i]<-"HBC"
    }else if(str_detect(t[1],"dBC")){
      foveal_meta$cell_type[i]<-"DBC"
    }else{
      foveal_meta$cell_type[i]<-t[1]
    }
  }else{
    
  }
}


colnames(org)<-organoid_meta$cell_type
colnames(per)<-peripheral_meta$cell_type
colnames(fove)<-foveal_meta$cell_type
cell_type<-unique(c(organoid_meta$cell_type,peripheral_meta$cell_type,foveal_meta$cell_type))
cell_type<-cell_type[-c(3,4,8)]
p<-c()
o<-c()
f<-c()
for (i in cell_type) {
  print(i)
  if(length(which(colnames(per)%in%i))>0){
    t<-per[,which(colnames(per)%in%i)]
    p<-cbind(p,apply(t, 1, mean))
    colnames(p)[length(p[1,])]<-i
  }
  if(length(which(colnames(org)%in%i))>0){
    t<-org[,which(colnames(org)%in%i)]
    o<-cbind(o,apply(t, 1, mean))
    colnames(o)[length(o[1,])]<-i
  }
  if(length(which(colnames(fove)%in%i))>0){
    t<-fove[,which(colnames(fove)%in%i)]
    f<-cbind(f,apply(t, 1, mean))
    colnames(f)[length(f[1,])]<-i
  }
}
p_sort<-c("RPE","MC","GC","C","R","CM","AST","NK","FB","AC","PER","uG","HBC","MO","MAST","HC","DBC","RBC","TCE","END")
f_sort<-c("RPE","MC","GC","C","R","BCE","AST","NK","FB","AC","PER","uG","HBC","MO","HC","DBC","RBC","TCE","END")
o_sort<-c("RPE","MC","C","R","AST","AC","HBC","HC","DBC","RBC")

allexpre<-cbind(f[,pmatch(f_sort,colnames(f))],p[,pmatch(p_sort,colnames(p))],o[,pmatch(o_sort,colnames(o))])
tmp<-c()
for (i in 1:nrow(allexpre)) {
  if(length(which(allexpre[i,]>0))==49){
    tmp<-c(tmp,i)
  }
}
library(pheatmap)
allexpre1<-allexpre[-tmp,]
allexpre1<-(allexpre-min(allexpre))/(max(allexpre)-min(allexpre))
allexpre1_ple<-allexpre1[which(rownames(allexpre1)%in%g1),]
allexpre1_nonple<-allexpre1[which(rownames(allexpre1)%in%g2),]
cv_ple<-c()
for (i in 1:length(allexpre1_ple[,1])) {
  cv<-sd(allexpre1_ple[i,])/mean(allexpre1_ple[i,])
  cv_ple<-c(cv_ple,cv)
}
cv_ple<-data.frame(CV=cv_ple,Symbol=rownames(allexpre1_ple))
cv_ple$CV<-as.character(cv_ple$CV)
cv_ple$CV<-as.numeric(cv_ple$CV)

cv_ple<-cv_ple[order(cv_ple$CV,decreasing=T),]
cv_ple<-cv_ple[which(cv_ple$CV>0.5),]


cv_nonple<-c()
for (i in 1:length(allexpre1_nonple[,1])) {
  cv<-sd(allexpre1_nonple[i,])/mean(allexpre1_nonple[i,])
  cv_nonple<-c(cv_nonple,cv)
}
cv_nonple<-data.frame(CV=cv_nonple,Symbol=rownames(allexpre1_nonple))
cv_nonple$CV<-as.character(cv_nonple$CV)
cv_nonple$CV<-as.numeric(cv_nonple$CV)

cv_nonple<-cv_nonple[order(cv_nonple$CV,decreasing=T),]
cv_nonple<-cv_nonple[which(cv_nonple$CV>0.5),]
allexpre1_ple<-allexpre1_ple[which(rownames(allexpre1_ple)%in%cv_ple$Symbol),]
allexpre1_nonple<-allexpre1_nonple[which(rownames(allexpre1_nonple)%in%cv_nonple$Symbol),]

p<-pheatmap(allexpre1_ple,cluster_rows = T,cluster_cols = F,clustering_method = "single")
n<-pheatmap(allexpre1_nonple,cluster_rows = T,cluster_cols = F,clustering_method = "single")
allexpre1_ple<-allexpre1_ple[p$tree_row$order,][-c(18:23),]
allexpre1_nonple<-allexpre1_nonple[n$tree_row$order,][c(1:6,11,12),]

p<-pheatmap(allexpre1_ple,cluster_rows = T,cluster_cols = F,clustering_method = "single")
n<-pheatmap(allexpre1_nonple,cluster_rows = T,cluster_cols = F,clustering_method = "single")
allexpre1<-rbind(allexpre1_ple[p$tree_row$order,],allexpre1_nonple[n$tree_row$order,])

annotation_col = data.frame("Celltype"=colnames(allexpre1) )
colnames(allexpre1)<-rownames(annotation_col)
rownames(annotation_col) = colnames(allexpre1)

ann_colors = list(Celltype=c("RPE"="#8200FE","MC"="#0C0155","GC"="#D45183","C"="#989C69","R"="#000000",
                  "BCE"="#C9DC80","AST"="#2901FE", "NK"="#FAFD44", "FB"="#9F45F2","AC"="#C71C21",
                  "PER"="#7E029C","uG"="#9B95F6","HBC"="#F8FD3B","MO"="#98995C","HC"="#C77D1F",
                  "DBC"="#558D37", "RBC"="#1B3506", "TCE"="#D1DB5D","END"="#4B0079","CM"="#CF01A1",
                  "MAST"="#BEC25F"))
pheatmap(allexpre1,cluster_rows = F,cluster_cols = F,annotation_col = annotation_col,annotation_colors = ann_colors,
         gaps_col = c(19,39),gaps_row = c(20),color = colorRampPalette(colors = c("white","#8E1E20"))(100),
         show_colnames = F,annotation_names_col = F,clustering_method = "complete")
setwd("/share/pub/chenfk/project/ZB_eyedisease/single/GSE104827/cluster_data/fig7")
ggsave('./fig7_redcolor.pdf', t, width = 10,height = 8)
