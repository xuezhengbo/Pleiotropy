library(data.table)
library(GenomicSEM)

#munge the summary statistics
munge(files=c("AMD_10pc_info_noMHC.txt", "DR_10pc_info_noMHC.txt", "GLC_10pc_info_noMHC.txt", "RD_10pc_info_noMHC.txt", "myopia_10pc_info_noMHC.txt"), hm3 = "w_hm3.noMHC.snplist",trait.names=c("AMD", "DR", "GLC", "RD", "Myopia"), N=c(66387,62229,68390,64011,64268), info.filter = 0.9, maf.filter = 0.01)

#run multivariable LDSC to create the S and V matrices
traits <- c("AMD.sumstats.gz", "DR.sumstats.gz", "GLC.sumstats.gz", "RD.sumstats.gz", "Myopia.sumstats.gz")
sample.prev <- c(0.0885, 0.0265, 0.1151, 0.0539, 0.4356)
population.prev <- c(0.123, 0.0176, 0.0251, 0.0002, 0.306)
ld <- "/share/pub/xuezb/software/ldsc/eur_w_ld_chr/"
wld <- "/share/pub/xuezb/software/ldsc/eur_w_ld_chr/"
trait.names <- c("AMD", "DR", "GLC", "RD", "Myopia")
anthro <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)

#CommonFactor model
CommonFactor_DWLS <- commonfactor(covstruc = anthro, estimation="DWLS")
CommonFactor_DWLS
CommonFactor_ML <- commonfactor(covstruc = anthro, estimation="ML")
CommonFactor_ML

################################# EFA and CFA #######################################
#smooth the S matrix for EFA using the nearPD function in the Matrix package. 
require(Matrix)
Ssmooth <- as.matrix((nearPD(anthro$S, corr = FALSE))$mat)

#run EFA with promax rotation and 2 factors using the factanal function in the stats package
require(stats)
EFA <- factanal(covmat = Ssmooth, factors = 2, rotation = "promax")
EFA$loadings

#Specify the Genomic confirmatory factor model
CFAofEFA <- 'F1 =~ NA*AMD + DR + GLC
             F2 =~ NA*RD + Myopia
F1~~F2
AMD ~~ a*AMD
a > .001
RD ~~ b*RD
b > .001'
#run the model
Anthro_DWLS <- usermodel(anthro, estimation = "DWLS", model = CFAofEFA, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
Anthro_ML <- usermodel(anthro, estimation = "ML", model = CFAofEFA, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
#print the results
Anthro_DWLS
Anthro_ML

