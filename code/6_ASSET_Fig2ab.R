library(ASSET)
library(data.table)

#read gwas summary, extract BETA and SE
gwas_AMD <- fread("/share/pub/xuezb/biobank/GWAS/AMD/data_68304/result/AMD_10pc.txt")
gwas_DR <- fread("/share/pub/xuezb/biobank/GWAS/DR/data_64032/result/DR_10pc.txt")
gwas_GLC <- fread("/share/pub/xuezb/biobank/GWAS/GLC/data_70359/result/GLC_10pc.txt")
gwas_RD <- fread("/share/pub/xuezb/biobank/GWAS/RD/data_65837/result/RD_10pc.txt")
gwas_myopia <- fread("/share/pub/xuezb/biobank/GWAS/myopia/data_66133/result/myopia_10pc.txt")
gwas_AMD <- gwas_AMD[,c(1,9,10)]
gwas_DR <- gwas_DR[,c(1,9,10)]
gwas_GLC <- gwas_GLC[,c(1,9,10)]
gwas_RD <- gwas_RD[,c(1,9,10)]
gwas_myopia <- gwas_myopia[,c(1,9,10)]
names(gwas_AMD) <- c("SNP", "AMD.BETA", "AMD.SE")
names(gwas_DR) <- c("SNP", "DR.BETA", "DR.SE")
names(gwas_GLC) <- c("SNP", "GLC.BETA", "GLC.SE")
names(gwas_RD) <- c("SNP", "RD.BETA", "RD.SE")
names(gwas_myopia) <- c("SNP", "Myopia.BETA", "Myopia.SE")

gwas_merge <- merge(merge(merge(merge(gwas_AMD, gwas_DR), gwas_GLC), gwas_RD), gwas_myopia)

snps <- as.vector(as.matrix(gwas_merge[, "SNP"]))
traits.lab <- c("AMD","DR","GLC","RD","Myopia")
beta.hat <- as.numeric(as.matrix(gwas_merge[, c(2,4,6,8,10)]))
sigma.hat <- as.numeric(as.matrix(gwas_merge[, c(3,5,7,9,11)]))

N11 <- as.matrix(fread("N11.txt", drop = 1))
N10 <- as.matrix(fread("N10.txt", drop = 1))
N00 <- as.matrix(fread("N00.txt", drop = 1))
rownames(N11) <- traits.lab
rownames(N10) <- traits.lab
rownames(N00) <- traits.lab
cor <- list(N11=N11, N00=N00, N10=N10)
ncase <- diag(N11)
ncntl <- diag(N00)

#asset
res <- h.traits(snps, traits.lab, beta.hat, sigma.hat, ncase, ncntl, cor=cor, cor.numr=FALSE, side = 2,
                meta=TRUE, zmax.args=NULL, meth.pval="DLM")
summary <- h.summary(res)
save.image("asset.RData")
