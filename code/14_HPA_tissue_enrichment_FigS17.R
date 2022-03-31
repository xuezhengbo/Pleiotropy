firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
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
gene<-read.csv("****",sep = "\t")
g1<-as.character(unique(gene$Peiotropy))
g1<-g1[which(g1%in%pcg$V2)]
g2<-as.character(unique(gene$Non.peiotropy))
g2<-g2[-which(g2%in%"")]
g2<-g2[which(g2%in%pcg$V2)]

HPA<-read.csv("/share/pub/chenfk/Data/HPA/tissue_enhanced_genes1.txt",sep = "\t",header = T)
HPA$RNA.tissue.specific<-as.character(HPA$RNA.tissue.specific)
HPA$RNA.tissue.specific<-firstup(HPA$RNA.tissue.specific)
HPA$RNA.tissue.specific<-gsub(" ","_",HPA$RNA.tissue.specific)
HPA$RNA.tissue.specific<-gsub(",","",HPA$RNA.tissue.specific)
HPA<-HPA[which(HPA$Gene%in%pcg$V2),]
tissue<-as.character(unique(HPA$RNA.tissue.specific))


allgene<-as.character(unique(pcg$V2))

ti<-c()
dd<-c()
final1<-c()
for (i in 1:length(tissue)) {
  line<-HPA$Gene[which(HPA$RNA.tissue.specific%in%tissue[i])]
  gene_overlap<-which(line%in%c(g1,g2))
  gene_overlap_other<-which(allgene[which(!(allgene%in%line))]%in%c(g1,g2))
  other_tissue_num<-length(allgene[which(!(allgene%in%line))])
  nmoudle_nonoverlap<-length(line)-length(gene_overlap)
  othermoudle_over<-length(gene_overlap_other)
  compare<-matrix(c(length(gene_overlap),othermoudle_over,length(line)-length(gene_overlap),length(allgene[which(!(allgene%in%line))])-othermoudle_over),nr = 2)
  
  result<-fisher.test(compare ,conf.level = 0.95)
  logP<--(log10(result$p.value))
  pvalue<-result$p.value  
  if(length(line[gene_overlap])==0){
    dd<-cbind(tissue[i],pvalue,logP,length(gene_overlap),length(line)-length(gene_overlap),othermoudle_over,length(allgene[which(!(allgene%in%line))])-othermoudle_over,"")
  }else{
    dd<-cbind(tissue[i],pvalue,logP,length(gene_overlap),length(line)-length(gene_overlap),othermoudle_over,length(allgene[which(!(allgene%in%line))])-othermoudle_over,paste(line[gene_overlap],collapse =  ";"))
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

write.table(final1,"**",sep = "\t",col.names = T,row.names = F,quote = F)
