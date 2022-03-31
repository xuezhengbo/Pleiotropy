#GTEx tissue-specific enrichment
library(ggplot2)
ensg2gene <- read.table("./ensg2gene")
gene_biotype <-read.table("./Homo_sapiens.GRCh38.95.genebiotype",
                          sep = "\t",
                          as.is = T,
                          col.names = c("gene", "biotype")
)
pcg <- gene_biotype[which(gene_biotype$biotype=="protein_coding"),]
pcg <- merge(pcg,ensg2gene,by.x="gene",by.y="V1")
###########################
gene<-read.csv('*.txt',sep = "\t",header = T)
g1<-as.character(unique(gene$symbol))
g1<-g1[which(g1%in%pcg$V2)]
g2<-as.character(unique(gene$Non.pleiotropy))
g2<-g2[-which(g2%in%"")]
g2<-g2[which(g2%in%pcg$V2)]
#####################################
ff1<-read.table("./Gtex/gtex55tissue_protein_new.txt",header = T,sep = "\t")
##################################################
allgene<-as.character(unique(rownames(ff1)))

gtex_t<-read.table("./Gtex/T_sta_info.txt",header = T,sep = "\t")
tissue<-read.table("./Gtex/tissue_name.txt",header = F,sep = "\t")
colnames(gtex_t)<-tissue$V1
tissue<-as.character(tissue$V1)
final1<-c()
for (i in 1:ncol(gtex_t)) {
  t<-gtex_t[,i]
  t<-t(matrix(as.numeric(unlist(strsplit(as.character(t),";"))),nrow = 3))
  
  t<-as.data.frame(t)
  
  t$fdr<-p.adjust(t$V3,method = "fdr")
  rownames(t)<-rownames(gtex_t)
  t1<-t[which(t$V2>3&t$fdr<0.05),]
  t1<-rownames(t1)
  
  dd<-c()
  
  gene_overlap<-which(t1%in%c(g1))
  gene_overlap_other<-which(allgene[which(!(allgene%in%t1))]%in%c(g1))
  nmoudle_nonoverlap<-length(t1)-length(gene_overlap)
  othermoudle_over<-length(gene_overlap_other)
  compare<-matrix(c(length(gene_overlap),othermoudle_over,nmoudle_nonoverlap,(length(allgene)-othermoudle_over-length(t1))),nr = 2)
  
  result<-fisher.test(compare ,conf.level = 0.95)
  logP<--(log10(result$p.value))
  pvalue<-result$p.value  
  if(length(t1[gene_overlap])==0){
    dd<-cbind(tissue[i],pvalue,logP,length(gene_overlap),othermoudle_over,nmoudle_nonoverlap,(length(allgene)-othermoudle_over-length(t1)),"")
  }else{
    dd<-cbind(tissue[i],pvalue,logP,length(gene_overlap),othermoudle_over,nmoudle_nonoverlap,(length(allgene)-othermoudle_over-length(t1)),paste(t1[gene_overlap],collapse =  ";"))
  }
  
  final1<-rbind(final1,dd)
}

final1<-as.data.frame(final1)
colnames(final1)[1]<-"Tissue"
final1$pvalue<-as.numeric(as.character(final1$pvalue))
final1$logP<-as.numeric(as.character(final1$logP))
final1<-final1[order(final1$pvalue),]
final1$p.adj<-p.adjust(final1$pvalue,method = "fdr")
final1<-final1[order(final1$p.adj),]
#num<-num+1
final1[1:10,1:3]


write.table(final1,"./Top50_gene_enrich.txt",sep = "\t",col.names = T,row.names = F,quote = F)
gtex<-read.table("./Top50_gene_enrich.txt",sep = "\t",header = T)
lable<-read.table("./Gtex/tissue_system.txt",header=T,sep="\t")
type<-read.table("./Gtex/tissue_type.txt",header=F,sep="\t")
newgtex<-c()
for (i in 1:55) {
  newgtex<-rbind(newgtex,gtex[which(gtex$Tissue%in%type$V1[i]),])
}
newgtex$Type<-type$V2

color<-cbind(as.character(unique(type$V2)),c("#F75E5E","#FFB366","#FFFF66","#B3FF66",
                                             "#66FF66","#66FFB3","#66FFFF","#66B3FF","#6666FF","#B366FF"))
newgtex$Color<-c()


for (i in 1:55) {
  newgtex$Color[i]<-color[which(color[,1]%in%newgtex$Type[i]),2]
}
newgtex$Tissue<-as.character(newgtex$Tissue)
newgtex1<-c()
for (i in 1:55) {
  newgtex1<-rbind(newgtex1,newgtex[which(newgtex$Tissue%in%newgtex$Tissue[56-i]),])
}
newgtex1$Tis<-"Retina"
newgtex1$logP<-as.numeric(as.character(newgtex1$logP))
p1<-ggplot(newgtex1,aes(x=Tissue,y=logP,fill=Type))+geom_bar(stat = "identity")+theme(
  panel.background = element_rect(fill = "transparent",colour = "black"),
  plot.background = element_rect(fill = "transparent",colour = NA),
  #panel.grid =element_blank(),
  panel.grid.major=element_line(colour="#EBEBEB"),
  axis.text= element_text(size=10, color="black", vjust = 0.5, hjust = 0.5),
  axis.text.x =element_text(size=10, color="black", hjust =1) )+
  xlab("Tissue")+ylab("-log10 p-value")+ scale_x_discrete(limits=factor(newgtex1[,1]))+
  theme(axis.title = element_text(size = 10, color = "black", vjust = 0.5, hjust = 0.5))
p1+coord_flip()+
  labs(title="Tissue specific Expression Enrichment")+
  theme(plot.title = element_text(hjust = 0.5))+ 
  scale_fill_manual("Tissue_Category", values = c("Cardiovascular" = "#F75E5E", 
                                                  "Digestive" = "#FFB366", 
                                                  "Exo-/Endocrine" = "#FFFF66",
                                                  "Hemic and Immune" = "#B3FF66", 
                                                  "Integumentary" = "#66FF66", 
                                                  "Musculoskeletal" = "#66FFB3",
                                                  "Nervous" = "#66FFFF", 
                                                  "Respiratory" = "#66B3FF", 
                                                  "Retina" = "#6666FF",
                                                  "Urogenital" = "#B366FF"))




##################################
ple=overlap_retinagene[which(overlap_retinagene%in%g1)]
nonple=overlap_retinagene[which(overlap_retinagene%in%g2)]
gtex_all<-ff1[which(rownames(ff1)%in%c(ple,nonple)),]
colnames(gtex_all)<-tissue
gtex_all1<-((2-(-2))*((gtex_all-min(gtex_all))/(max(gtex_all)-min(gtex_all))))+(-2)

gtex_ple1<-gtex_all1[which(rownames(gtex_all1)%in%ple),]
gtex_nonple1<-gtex_all1[which(rownames(gtex_all1)%in%nonple),]


t<-pheatmap(gtex_ple1,color = colorRampPalette(colors = c("blue","white","red"))(100),
            cluster_rows = T,cluster_cols = T,treeheight_row = 0,treeheight_col = 0,clustering_method = 'centroid')
t1<-pheatmap(gtex_nonple1,color = colorRampPalette(colors = c("blue","white","red"))(100),
             cluster_rows = T,cluster_cols = T,treeheight_row = 0,treeheight_col = 0,clustering_method = 'centroid')




mmp<-gtex_ple1[t$tree_row$order,]
mmn<-gtex_nonple1[t1$tree_row$order,]
gtex_all<-rbind(mmp,mmn)
ann_colors = list(
  Gene = c("Pleiotropy" = "#EC870E", "Non-Pleiotropy"="#205AA7")
)
annotation_row = data.frame(Gene=as.factor(c(rep("Pleiotropy",6),rep("Non-Pleiotropy",6))) )
rownames(annotation_row) = rownames(gtex_all)
pheatmap(gtex_all,color = colorRampPalette(colors = c("blue","white","#EF2A29"))(100),annotation_row = annotation_row,
         annotation_colors = ann_colors,cluster_rows = F,cluster_cols = F,treeheight_row = 0,treeheight_col = 0,clustering_method = 'centroid')



