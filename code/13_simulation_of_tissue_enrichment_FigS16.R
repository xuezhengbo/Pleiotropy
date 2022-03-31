
ensg2gene <- read.table("./ensg2gene")
gene_biotype <-read.table("./Homo_sapiens.GRCh38.95.genebiotype",
                          sep = "\t",
                          as.is = T,
             col.names = c("gene", "biotype")
)
pcg <- gene_biotype[which(gene_biotype$biotype=="protein_coding"),]
pcg <- merge(pcg,ensg2gene,by.x="gene",by.y="V1")
###########################
gene<-read.csv('./GTEX/gene_list_1025.txt',sep = "\t",header = T)
g1<-as.character(unique(gene$Pleiotropy))
g1<-g1[which(g1%in%pcg$V2)]
gene1<-gene[which(gene$Pleiotropy%in%g1),]
g2<-as.character(unique(gene$Non.pleiotropy))
g2<-g2[-which(g2%in%"")]
g2<-g2[which(g2%in%pcg$V2)]
gene2<-gene[which(gene$Non.pleiotropy%in%g2),]
#####################################
ff1<-read.table("./Gtex/gtex55tissue_protein_new.txt",header = T,sep = "\t")
##################################################
allgene<-as.character(unique(rownames(ff1)))

gtex_t<-read.table("./Gtex/T_sta_info.txt",header = T,sep = "\t")
tissue<-read.table("./Gtex/tissue_name.txt",header = F,sep = "\t")
colnames(gtex_t)<-tissue$V1
tissue<-as.character(tissue$V1)
final1_random<-c()

for (j in 1:1000) {
  print(j)
  final1<-c()
for (i in 1:1) {
  t<-gtex_t[,i]
  t<-t(matrix(as.numeric(unlist(strsplit(as.character(t),";"))),nrow = 3))
  
  t<-as.data.frame(t)
  gg<-rownames(ff1)[sample(x =1:17795 ,size = length(c(g1,g2)))]
  t$fdr<-p.adjust(t$V3,method = "fdr")
  rownames(t)<-rownames(gtex_t)
  t1<-t[which(t$V2>3&t$fdr<0.05),]
  t1<-rownames(t1)
  
  dd<-c()
  
  gene_overlap<-which(t1%in%c(gg))
  gene_overlap_other<-which(allgene[which(!(allgene%in%t1))]%in%c(gg))
  nmoudle_nonoverlap<-length(t1)-length(gene_overlap)
  othermoudle_over<-length(gene_overlap_other)
  compare<-matrix(c(length(gene_overlap),othermoudle_over,nmoudle_nonoverlap,(length(allgene)-othermoudle_over-length(t1))),nr = 2)
  
  result<-fisher.test(compare ,conf.level = 0.95)
  logP<--(log10(result$p.value))
  pvalue<-result$p.value  
  OR<-result$estimate
  if(length(t1[gene_overlap])==0){
    dd<-cbind(tissue[i],pvalue,logP,OR,length(gene_overlap),othermoudle_over,nmoudle_nonoverlap,(length(allgene)-othermoudle_over-length(t1)),"")
  }else{
    dd<-cbind(tissue[i],pvalue,logP,OR,length(gene_overlap),othermoudle_over,nmoudle_nonoverlap,(length(allgene)-othermoudle_over-length(t1)),paste(t1[gene_overlap],collapse =  ";"))
  }
  final1<-rbind(final1,dd)
}

final1<-as.data.frame(final1)
colnames(final1)[1]<-"Tissue"
final1$pvalue<-as.numeric(as.character(final1$pvalue))
final1$logP<-as.numeric(as.character(final1$logP))
final1<-final1[order(final1$pvalue),]
final1$p.adj<-p.adjust(final1$pvalue,method = "fdr")
cc<-cbind(final1$pvalue,as.character(final1$OR),final1$p.adj)
final1_random<-rbind(final1_random,cc)
}


final1_random<-as.data.frame(final1_random)

final1_random$V2<-as.numeric(as.character(final1_random$V2))

hist(final1_random$V2)
E<-hist(final1_random$V2,freq = T,lwd = 3,
        col="#CF242A",
        cex.axis = 1.5,
        cex.lab = 2,
        font.lab = 1,
        font.axis = 2,border = F ,breaks =8,
        xlab = "Enrichment score (OR)",main = "")
abline(v=3.45014616209649,col="#5BABED", lty=2,lwd=3)
