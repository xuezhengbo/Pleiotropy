#GTEx tissue-specific enrichment
#get tissue-specific genes in GTEx tissues
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


  for (i in 1:length(tissue)) {
    t<-gtex_t[,i]
    t<-t(matrix(as.numeric(unlist(strsplit(as.character(t),";"))),nrow = 3))
    t<-as.data.frame(t)
    t$fdr<-p.adjust(t$V3,method = "fdr")
    rownames(t)<-rownames(gtex_t)
    t1<-t[which(t$V2>3&t$fdr<0.05),]
    t1<-rownames(t1)
    loca<-paste("./Gtex/GTEx_tissue_gene/",tissue[i],".txt")
    write.table(t1,loca,sep = "\t",row.names = F,col.names = F,quote = F)
  }
