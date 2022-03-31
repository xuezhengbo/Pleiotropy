library(qqman)
library(data.table)

args <- commandArgs(trailingOnly=TRUE)
path <- args[1]
gwasResults <- fread(path, header=T, sep="\t")
         
#calculate lambda
p_value <- gwasResults$P
z <- qnorm(p_value/ 2)
lambda <- round(median(z^2, na.rm = TRUE) / qchisq(0.5,1), 3)
print(lambda)
case <- as.numeric(args[2])
control <- as.numeric(args[3])
lambda1000 <- 1 + (lambda - 1) * (1 / case + 1 / control) * 500
print(lambda1000)

color_set <- c("#4169E1","#87CEFA")
#manhattan
png(args[4])
manhattan(gwasResults,chr = "CHR", bp = "POS", p = "P", snp = "SNP", col = color_set, ylim = c(0,35), cex.axis = 1, cex.lab = 1)
while (!is.null(dev.list())) dev.off()

#QQplot
png(args[5])
qq(gwasResults$P, cex.axis = 1, cex.lab = 1)
while (!is.null(dev.list())) dev.off()


