setwd("H:/FUMA_ASSET/eQTL_plot")
library(locuscomparer)
library(data.table)
library(dplyr)

eqtl <- fread("H:/FUMA_ASSET_with_not_sig_eqtl/eqtl.txt")
credible_snp <- fread("H:/FUMA_ASSET/Credible_SNP_func.txt")

eqtl_credible_snp <- merge(eqtl[,c(1,3,6,13)], credible_snp[,c(1,9,11)])

#FGF5
FGF5_eqtl_Esophagus_Muscularis <- eqtl_credible_snp %>% 
  filter(symbol == "FGF5" & tissue == "Esophagus_Muscularis")
FGF5_eqtl_Cells_Cultured_fibroblasts <- eqtl_credible_snp %>% 
  filter(symbol == "FGF5" & tissue == "Cells_Cultured_fibroblasts")
fwrite(FGF5_eqtl_Esophagus_Muscularis[,c(5,6)], "FGF5_Esophagus_Muscularis_gwas.tsv", sep = "\t")
fwrite(FGF5_eqtl_Esophagus_Muscularis[,c(5,3)], "FGF5_Esophagus_Muscularis_eqtl.tsv", sep = "\t")
fwrite(FGF5_eqtl_Cells_Cultured_fibroblasts[,c(5,6)], "FGF5_Cells_Cultured_fibroblasts_gwas.tsv", sep = "\t")
fwrite(FGF5_eqtl_Cells_Cultured_fibroblasts[,c(5,3)], "FGF5_Cells_Cultured_fibroblasts_eqtl.tsv", sep = "\t")
#FGF5 Esophagus_Muscularis
gwas_fn = 'FGF5_Esophagus_Muscularis_gwas.tsv'
eqtl_fn = 'FGF5_Esophagus_Muscularis_eqtl.tsv'
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Meta-GWAS', title2 = 'Esophagus Muscularis eQTL',
             marker_col1 = "SNP", pval_col1 = "P",marker_col2 = "SNP", pval_col2 = "p",
             snp = "rs7678123", legend_position = "topleft", combine = F)
#FGF5 Cells_Cultured_fibroblasts
gwas_fn = 'FGF5_Cells_Cultured_fibroblasts_gwas.tsv'
eqtl_fn = 'FGF5_Cells_Cultured_fibroblasts_eqtl.tsv'
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Meta-GWAS', title2 = 'Cells Cultured Fibroblasts eQTL',
             marker_col1 = "SNP", pval_col1 = "P",marker_col2 = "SNP", pval_col2 = "p",
             snp = "rs2055178", legend_position = "topright", combine = F)

#PRKG2
PRKG2_eqtl_BIOS <- eqtl_credible_snp %>% 
  filter(symbol == "PRKG2" & tissue == "BIOS_eQTL_geneLevel")
PRKG2_eqtl_Lung <- eqtl_credible_snp %>% 
  filter(symbol == "PRKG2" & tissue == "Lung")
PRKG2_eqtl_Skin_Suprapubic <- eqtl_credible_snp %>% 
  filter(symbol == "PRKG2" & tissue == "Skin_Not_Sun_Exposed_Suprapubic")
PRKG2_eqtl_Skin_Lower_leg <- eqtl_credible_snp %>% 
  filter(symbol == "PRKG2" & tissue == "Skin_Sun_Exposed_Lower_leg")
fwrite(PRKG2_eqtl_BIOS[,c(5,6)], "PRKG2_BIOS_gwas.tsv", sep = "\t")
fwrite(PRKG2_eqtl_BIOS[,c(5,3)], "PRKG2_BIOS_eqtl.tsv", sep = "\t")
fwrite(PRKG2_eqtl_Lung[,c(5,6)], "PRKG2_Lung_gwas.tsv", sep = "\t")
fwrite(PRKG2_eqtl_Lung[,c(5,3)], "PRKG2_Lung_eqtl.tsv", sep = "\t")
fwrite(PRKG2_eqtl_Skin_Suprapubic[,c(5,6)], "PRKG2_Skin_Suprapubic_gwas.tsv", sep = "\t")
fwrite(PRKG2_eqtl_Skin_Suprapubic[,c(5,3)], "PRKG2_Skin_Suprapubic_eqtl.tsv", sep = "\t")
fwrite(PRKG2_eqtl_Skin_Lower_leg[,c(5,6)], "PRKG2_Skin_Lower_leg_gwas.tsv", sep = "\t")
fwrite(PRKG2_eqtl_Skin_Lower_leg[,c(5,3)], "PRKG2_Skin_Lower_leg_eqtl.tsv", sep = "\t")
#PRKG2 BIOS
gwas_fn = 'PRKG2_BIOS_gwas.tsv'
eqtl_fn = 'PRKG2_BIOS_eqtl.tsv'
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Meta-GWAS', title2 = 'BIOS blood eQTL',
             marker_col1 = "SNP", pval_col1 = "P",marker_col2 = "SNP", pval_col2 = "p",
             snp = NULL, legend_position = "bottomright", combine = F)
#PRKG2 Lung
gwas_fn = 'PRKG2_Lung_gwas.tsv'
eqtl_fn = 'PRKG2_Lung_eqtl.tsv'
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Meta-GWAS', title2 = 'Lung eQTL',
             marker_col1 = "SNP", pval_col1 = "P",marker_col2 = "SNP", pval_col2 = "p",
             snp = "rs2055178", legend_position = "bottomright", combine = F)
#PRKG2 Skin_Not_Sun_Exposed_Suprapubic
gwas_fn = 'PRKG2_Skin_Suprapubic_gwas.tsv'
eqtl_fn = 'PRKG2_Skin_Suprapubic_eqtl.tsv'
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Meta-GWAS', title2 = 'Skin_Suprapubic eQTL',
             marker_col1 = "SNP", pval_col1 = "P",marker_col2 = "SNP", pval_col2 = "p",
             snp = "rs2055178", legend_position = "bottomright", combine = F)
#PRKG2 Skin_Lower_leg
gwas_fn = 'PRKG2_Skin_Lower_leg_gwas.tsv'
eqtl_fn = 'PRKG2_Skin_Lower_leg_eqtl.tsv'
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Meta-GWAS', title2 = 'Skin_Lower_leg eQTL',
             marker_col1 = "SNP", pval_col1 = "P",marker_col2 = "SNP", pval_col2 = "p",
             snp = "rs2055178", legend_position = "bottomright", combine = F)

#BMP3
BMP3_eqtl_Adipose_Subcutaneous <- eqtl_credible_snp %>% 
  filter(symbol == "BMP3" & tissue == "Adipose_Subcutaneous")
BMP3_eqtl_Skin_Suprapubic <- eqtl_credible_snp %>% 
  filter(symbol == "BMP3" & tissue == "Skin_Not_Sun_Exposed_Suprapubic")
BMP3_eqtl_Skin_Lower_leg <- eqtl_credible_snp %>% 
  filter(symbol == "BMP3" & tissue == "Skin_Sun_Exposed_Lower_leg")
fwrite(BMP3_eqtl_Adipose_Subcutaneous[,c(5,6)], "BMP3_Adipose_Subcutaneous_gwas.tsv", sep = "\t")
fwrite(BMP3_eqtl_Adipose_Subcutaneous[,c(5,3)], "BMP3_Adipose_Subcutaneous_eqtl.tsv", sep = "\t")
fwrite(BMP3_eqtl_Skin_Suprapubic[,c(5,6)], "BMP3_Skin_Suprapubic_gwas.tsv", sep = "\t")
fwrite(BMP3_eqtl_Skin_Suprapubic[,c(5,3)], "BMP3_Skin_Suprapubic_eqtl.tsv", sep = "\t")
fwrite(BMP3_eqtl_Skin_Lower_leg[,c(5,6)], "BMP3_Skin_Lower_leg_gwas.tsv", sep = "\t")
fwrite(BMP3_eqtl_Skin_Lower_leg[,c(5,3)], "BMP3_Skin_Lower_leg_eqtl.tsv", sep = "\t")
#BMP3 Adipose_Subcutaneous
gwas_fn = 'BMP3_Adipose_Subcutaneous_gwas.tsv'
eqtl_fn = 'BMP3_Adipose_Subcutaneous_eqtl.tsv'
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Meta-GWAS', title2 = 'Adipose_Subcutaneous eQTL',
             marker_col1 = "SNP", pval_col1 = "P",marker_col2 = "SNP", pval_col2 = "p",
             snp = "rs7678123", legend_position = "topright", combine = F)
#BMP3 Skin_Not_Sun_Exposed_Suprapubic
gwas_fn = 'BMP3_Skin_Suprapubic_gwas.tsv'
eqtl_fn = 'BMP3_Skin_Suprapubic_eqtl.tsv'
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Meta-GWAS', title2 = 'Skin_Suprapubic eQTL',
             marker_col1 = "SNP", pval_col1 = "P",marker_col2 = "SNP", pval_col2 = "p",
             snp = "rs7678123", legend_position = "topright", combine = F)
#BMP3 Skin_Lower_leg
gwas_fn = 'BMP3_Skin_Lower_leg_gwas.tsv'
eqtl_fn = 'BMP3_Skin_Lower_leg_eqtl.tsv'
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Meta-GWAS', title2 = 'Skin_Lower_leg eQTL',
             marker_col1 = "SNP", pval_col1 = "P",marker_col2 = "SNP", pval_col2 = "p",
             snp = "rs7678123", legend_position = "topleft", combine = F)


#DPYSL4
DPYSL4_eqtl_Adipose_Subcutaneous <- eqtl_credible_snp %>% 
  filter(symbol == "DPYSL4" & tissue == "Adipose_Subcutaneous")
DPYSL4_eqtl_Artery_Tibial <- eqtl_credible_snp %>% 
  filter(symbol == "DPYSL4" & tissue == "Artery_Tibial")
DPYSL4_eqtl_Nerve_Tibial <- eqtl_credible_snp %>% 
  filter(symbol == "DPYSL4" & tissue == "Nerve_Tibial")
DPYSL4_eqtl_Skin_Lower_leg <- eqtl_credible_snp %>% 
  filter(symbol == "DPYSL4" & tissue == "Skin_Sun_Exposed_Lower_leg")
fwrite(DPYSL4_eqtl_Adipose_Subcutaneous[,c(5,6)], "DPYSL4_Adipose_Subcutaneous_gwas.tsv", sep = "\t")
fwrite(DPYSL4_eqtl_Adipose_Subcutaneous[,c(5,3)], "DPYSL4_Adipose_Subcutaneous_eqtl.tsv", sep = "\t")
fwrite(DPYSL4_eqtl_Artery_Tibial[,c(5,6)], "DPYSL4_Artery_Tibial_gwas.tsv", sep = "\t")
fwrite(DPYSL4_eqtl_Artery_Tibial[,c(5,3)], "DPYSL4_Artery_Tibial_eqtl.tsv", sep = "\t")
fwrite(DPYSL4_eqtl_Nerve_Tibial[,c(5,6)], "DPYSL4_Nerve_Tibial_gwas.tsv", sep = "\t")
fwrite(DPYSL4_eqtl_Nerve_Tibial[,c(5,3)], "DPYSL4_Nerve_Tibial_eqtl.tsv", sep = "\t")
fwrite(DPYSL4_eqtl_Skin_Lower_leg[,c(5,6)], "DPYSL4_Skin_Lower_leg_gwas.tsv", sep = "\t")
fwrite(DPYSL4_eqtl_Skin_Lower_leg[,c(5,3)], "DPYSL4_Skin_Lower_leg_eqtl.tsv", sep = "\t")
#DPYSL4 Adipose_Subcutaneous
gwas_fn = 'DPYSL4_Adipose_Subcutaneous_gwas.tsv'
eqtl_fn = 'DPYSL4_Adipose_Subcutaneous_eqtl.tsv'
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Meta-GWAS', title2 = 'Adipose_Subcutaneous eQTL',
             marker_col1 = "SNP", pval_col1 = "P",marker_col2 = "SNP", pval_col2 = "p",
             snp = "rs7908634", legend_position = "bottomright", combine = F)
#DPYSL4 Artery_Tibial
gwas_fn = 'DPYSL4_Artery_Tibial_gwas.tsv'
eqtl_fn = 'DPYSL4_Artery_Tibial_eqtl.tsv'
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Meta-GWAS', title2 = 'Artery_Tibial eQTL',
             marker_col1 = "SNP", pval_col1 = "P",marker_col2 = "SNP", pval_col2 = "p",
             snp = "rs7908634", legend_position = "bottomright", combine = F)
#DPYSL4 Nerve_Tibial
gwas_fn = 'DPYSL4_Nerve_Tibial_gwas.tsv'
eqtl_fn = 'DPYSL4_Nerve_Tibial_eqtl.tsv'
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Meta-GWAS', title2 = 'Nerve_Tibial eQTL',
             marker_col1 = "SNP", pval_col1 = "P",marker_col2 = "SNP", pval_col2 = "p",
             snp = "rs7908634", legend_position = "topright", combine = F)
#DPYSL4 Skin_Lower_leg
gwas_fn = 'DPYSL4_Skin_Lower_leg_gwas.tsv'
eqtl_fn = 'DPYSL4_Skin_Lower_leg_eqtl.tsv'
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Meta-GWAS', title2 = 'Skin_Lower_leg eQTL',
             marker_col1 = "SNP", pval_col1 = "P",marker_col2 = "SNP", pval_col2 = "p",
             snp = "rs7908634", legend_position = "topright", combine = F)

#STK32C
STK32C_eqtl_Adipose_Subcutaneous <- eqtl_credible_snp %>% 
  filter(symbol == "STK32C" & tissue == "Adipose_Subcutaneous")
STK32C_eqtl_Whole_Blood <- eqtl_credible_snp %>% 
  filter(symbol == "STK32C" & tissue == "Whole_Blood")
STK32C_eqtl_Artery_Aorta <- eqtl_credible_snp %>% 
  filter(symbol == "STK32C" & tissue == "Artery_Aorta")
STK32C_eqtl_Breast_Mammary_Tissue <- eqtl_credible_snp %>% 
  filter(symbol == "STK32C" & tissue == "Breast_Mammary_Tissue")
STK32C_eqtl_Muscle_Skeletal <- eqtl_credible_snp %>% 
  filter(symbol == "STK32C" & tissue == "Muscle_Skeletal")
STK32C_eqtl_Nerve_Tibial <- eqtl_credible_snp %>% 
  filter(symbol == "STK32C" & tissue == "Nerve_Tibial")
fwrite(STK32C_eqtl_Adipose_Subcutaneous[,c(5,6)], "STK32C_Adipose_Subcutaneous_gwas.tsv", sep = "\t")
fwrite(STK32C_eqtl_Adipose_Subcutaneous[,c(5,3)], "STK32C_Adipose_Subcutaneous_eqtl.tsv", sep = "\t")
fwrite(STK32C_eqtl_Whole_Blood[,c(5,6)], "STK32C_Whole_Blood_gwas.tsv", sep = "\t")
fwrite(STK32C_eqtl_Whole_Blood[,c(5,3)], "STK32C_Whole_Blood_eqtl.tsv", sep = "\t")
fwrite(STK32C_eqtl_Artery_Aorta[,c(5,6)], "STK32C_Artery_Aorta_gwas.tsv", sep = "\t")
fwrite(STK32C_eqtl_Artery_Aorta[,c(5,3)], "STK32C_Artery_Aorta_eqtl.tsv", sep = "\t")
fwrite(STK32C_eqtl_Breast_Mammary_Tissue[,c(5,6)], "STK32C_Breast_Mammary_Tissue_gwas.tsv", sep = "\t")
fwrite(STK32C_eqtl_Breast_Mammary_Tissue[,c(5,3)], "STK32C_Breast_Mammary_Tissue_eqtl.tsv", sep = "\t")
fwrite(STK32C_eqtl_Muscle_Skeletal[,c(5,6)], "STK32C_Muscle_Skeletal_gwas.tsv", sep = "\t")
fwrite(STK32C_eqtl_Muscle_Skeletal[,c(5,3)], "STK32C_Muscle_Skeletal_eqtl.tsv", sep = "\t")
fwrite(STK32C_eqtl_Nerve_Tibial[,c(5,6)], "STK32C_Nerve_Tibial_gwas.tsv", sep = "\t")
fwrite(STK32C_eqtl_Nerve_Tibial[,c(5,3)], "STK32C_Nerve_Tibial_eqtl.tsv", sep = "\t")
#STK32C Adipose_Subcutaneous
gwas_fn = 'STK32C_Adipose_Subcutaneous_gwas.tsv'
eqtl_fn = 'STK32C_Adipose_Subcutaneous_eqtl.tsv'
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Meta-GWAS', title2 = 'Adipose_Subcutaneous eQTL',
             marker_col1 = "SNP", pval_col1 = "P",marker_col2 = "SNP", pval_col2 = "p",
             snp = "rs7908634", legend_position = "bottomright", combine = F)
#STK32C Artery_Aorta
gwas_fn = 'STK32C_Artery_Aorta_gwas.tsv'
eqtl_fn = 'STK32C_Artery_Aorta_eqtl.tsv'
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Meta-GWAS', title2 = 'Artery_Aorta eQTL',
             marker_col1 = "SNP", pval_col1 = "P",marker_col2 = "SNP", pval_col2 = "p",
             snp = "rs7908634", legend_position = "bottomright", combine = F)
#STK32C Whole_Blood
gwas_fn = 'STK32C_Whole_Blood_gwas.tsv'
eqtl_fn = 'STK32C_Whole_Blood_eqtl.tsv'
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Meta-GWAS', title2 = 'Whole_Blood eQTL',
             marker_col1 = "SNP", pval_col1 = "P",marker_col2 = "SNP", pval_col2 = "p",
             snp = "rs7908634", legend_position = "bottomright", combine = F)
#STK32C Breast_Mammary_Tissue
gwas_fn = 'STK32C_Breast_Mammary_Tissue_gwas.tsv'
eqtl_fn = 'STK32C_Breast_Mammary_Tissue_eqtl.tsv'
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Meta-GWAS', title2 = 'Breast_Mammary_Tissue eQTL',
             marker_col1 = "SNP", pval_col1 = "P",marker_col2 = "SNP", pval_col2 = "p",
             snp = "rs7908634", legend_position = "bottomright", combine = F)
#STK32C Muscle_Skeletal
gwas_fn = 'STK32C_Muscle_Skeletal_gwas.tsv'
eqtl_fn = 'STK32C_Muscle_Skeletal_eqtl.tsv'
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Meta-GWAS', title2 = 'Muscle_Skeletal eQTL',
             marker_col1 = "SNP", pval_col1 = "P",marker_col2 = "SNP", pval_col2 = "p",
             snp = "rs7908634", legend_position = "bottomright", combine = F)
#STK32C Nerve_Tibial
gwas_fn = 'STK32C_Nerve_Tibial_gwas.tsv'
eqtl_fn = 'STK32C_Nerve_Tibial_eqtl.tsv'
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Meta-GWAS', title2 = 'Nerve_Tibial eQTL',
             marker_col1 = "SNP", pval_col1 = "P",marker_col2 = "SNP", pval_col2 = "p",
             snp = "rs7908634", legend_position = "bottomright", combine = F)

#LRRC27
LRRC27_eqtl_Brain_Cortex <- eqtl_credible_snp %>% 
  filter(symbol == "LRRC27" & tissue == "Brain_Cortex")
LRRC27_eqtl_Whole_Blood <- eqtl_credible_snp %>% 
  filter(symbol == "LRRC27" & tissue == "Whole_Blood")
LRRC27_eqtl_Brain_Frontal_Cortex_BA9 <- eqtl_credible_snp %>% 
  filter(symbol == "LRRC27" & tissue == "Brain_Frontal_Cortex_BA9")
LRRC27_eqtl_Cells_Cultured_fibroblasts <- eqtl_credible_snp %>% 
  filter(symbol == "LRRC27" & tissue == "Cells_Cultured_fibroblasts")
LRRC27_eqtl_Brain_Caudate_basal_ganglia <- eqtl_credible_snp %>% 
  filter(symbol == "LRRC27" & tissue == "Brain_Caudate_basal_ganglia")
LRRC27_eqtl_Nerve_Tibial <- eqtl_credible_snp %>% 
  filter(symbol == "LRRC27" & tissue == "Nerve_Tibial")
fwrite(LRRC27_eqtl_Brain_Cortex[,c(5,6)], "LRRC27_Brain_Cortex_gwas.tsv", sep = "\t")
fwrite(LRRC27_eqtl_Brain_Cortex[,c(5,3)], "LRRC27_Brain_Cortex_eqtl.tsv", sep = "\t")
fwrite(LRRC27_eqtl_Whole_Blood[,c(5,6)], "LRRC27_Whole_Blood_gwas.tsv", sep = "\t")
fwrite(LRRC27_eqtl_Whole_Blood[,c(5,3)], "LRRC27_Whole_Blood_eqtl.tsv", sep = "\t")
fwrite(LRRC27_eqtl_Cells_Cultured_fibroblasts[,c(5,6)], "LRRC27_Cells_Cultured_fibroblasts_gwas.tsv", sep = "\t")
fwrite(LRRC27_eqtl_Cells_Cultured_fibroblasts[,c(5,3)], "LRRC27_Cells_Cultured_fibroblasts_eqtl.tsv", sep = "\t")
fwrite(LRRC27_eqtl_Nerve_Tibial[,c(5,6)], "LRRC27_Nerve_Tibial_gwas.tsv", sep = "\t")
fwrite(LRRC27_eqtl_Nerve_Tibial[,c(5,3)], "LRRC27_Nerve_Tibial_eqtl.tsv", sep = "\t")
#LRRC27 Brain_Cortex
gwas_fn = 'LRRC27_Brain_Cortex_gwas.tsv'
eqtl_fn = 'LRRC27_Brain_Cortex_eqtl.tsv'
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Meta-GWAS', title2 = 'Brain_Cortex eQTL',
             marker_col1 = "SNP", pval_col1 = "P",marker_col2 = "SNP", pval_col2 = "p",
             snp = "rs7908634", legend_position = "topright", combine = F)
#LRRC27 Cells_Cultured_fibroblasts
gwas_fn = 'LRRC27_Cells_Cultured_fibroblasts_gwas.tsv'
eqtl_fn = 'LRRC27_Cells_Cultured_fibroblasts_eqtl.tsv'
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Meta-GWAS', title2 = 'Cells_Cultured_fibroblasts eQTL',
             marker_col1 = "SNP", pval_col1 = "P",marker_col2 = "SNP", pval_col2 = "p",
             snp = "rs7908634", legend_position = "bottomright", combine = F)
#LRRC27 Nerve_Tibial
gwas_fn = 'LRRC27_Nerve_Tibial_gwas.tsv'
eqtl_fn = 'LRRC27_Nerve_Tibial_eqtl.tsv'
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Meta-GWAS', title2 = 'Nerve_Tibial eQTL',
             marker_col1 = "SNP", pval_col1 = "P",marker_col2 = "SNP", pval_col2 = "p",
             snp = "rs7908634", legend_position = "bottomright", combine = F)

#PWWP2B
PWWP2B_eqtl_Adipose_Subcutaneous <- eqtl_credible_snp %>% 
  filter(symbol == "PWWP2B" & tissue == "Adipose_Subcutaneous")
PWWP2B_eqtl_Adipose_Visceral_Omentum <- eqtl_credible_snp %>% 
  filter(symbol == "PWWP2B" & tissue == "Adipose_Visceral_Omentum")
PWWP2B_eqtl_Muscle_Skeletal <- eqtl_credible_snp %>% 
  filter(symbol == "PWWP2B" & tissue == "Muscle_Skeletal")
PWWP2B_eqtl_Pancreas <- eqtl_credible_snp %>% 
  filter(symbol == "PWWP2B" & tissue == "Pancreas")
fwrite(PWWP2B_eqtl_Adipose_Subcutaneous[,c(5,6)], "PWWP2B_Adipose_Subcutaneous_gwas.tsv", sep = "\t")
fwrite(PWWP2B_eqtl_Adipose_Subcutaneous[,c(5,3)], "PWWP2B_Adipose_Subcutaneous_eqtl.tsv", sep = "\t")
fwrite(PWWP2B_eqtl_Adipose_Visceral_Omentum[,c(5,6)], "PWWP2B_Adipose_Visceral_Omentum_gwas.tsv", sep = "\t")
fwrite(PWWP2B_eqtl_Adipose_Visceral_Omentum[,c(5,3)], "PWWP2B_Adipose_Visceral_Omentum_eqtl.tsv", sep = "\t")
fwrite(PWWP2B_eqtl_Muscle_Skeletal[,c(5,6)], "PWWP2B_Muscle_Skeletal_gwas.tsv", sep = "\t")
fwrite(PWWP2B_eqtl_Muscle_Skeletal[,c(5,3)], "PWWP2B_Muscle_Skeletal_eqtl.tsv", sep = "\t")
fwrite(PWWP2B_eqtl_Pancreas[,c(5,6)], "PWWP2B_Pancreas_gwas.tsv", sep = "\t")
fwrite(PWWP2B_eqtl_Pancreas[,c(5,3)], "PWWP2B_Pancreas_eqtl.tsv", sep = "\t")
#PWWP2B Adipose_Subcutaneous
gwas_fn = 'PWWP2B_Adipose_Subcutaneous_gwas.tsv'
eqtl_fn = 'PWWP2B_Adipose_Subcutaneous_eqtl.tsv'
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Meta-GWAS', title2 = 'Adipose_Subcutaneous eQTL',
             marker_col1 = "SNP", pval_col1 = "P",marker_col2 = "SNP", pval_col2 = "p",
             snp = "rs7908634", legend_position = "topright", combine = F)
#PWWP2B Adipose_Visceral_Omentum
gwas_fn = 'PWWP2B_Adipose_Visceral_Omentum_gwas.tsv'
eqtl_fn = 'PWWP2B_Adipose_Visceral_Omentum_eqtl.tsv'
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Meta-GWAS', title2 = 'Adipose_Visceral_Omentum eQTL',
             marker_col1 = "SNP", pval_col1 = "P",marker_col2 = "SNP", pval_col2 = "p",
             snp = NULL, legend_position = "topright", combine = F)
#PWWP2B Muscle_Skeletal
gwas_fn = 'PWWP2B_Muscle_Skeletal_gwas.tsv'
eqtl_fn = 'PWWP2B_Muscle_Skeletal_eqtl.tsv'
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Meta-GWAS', title2 = 'Muscle_Skeletal eQTL',
             marker_col1 = "SNP", pval_col1 = "P",marker_col2 = "SNP", pval_col2 = "p",
             snp = "rs7908634", legend_position = "bottomright", combine = F)
#PWWP2B Pancreas
gwas_fn = 'PWWP2B_Pancreas_gwas.tsv'
eqtl_fn = 'PWWP2B_Pancreas_eqtl.tsv'
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Meta-GWAS', title2 = 'Pancreas eQTL',
             marker_col1 = "SNP", pval_col1 = "P",marker_col2 = "SNP", pval_col2 = "p",
             snp = "rs7908634", legend_position = "topright", combine = F)

#ME3
ME3_eqtl_EyeGEx <- eqtl_credible_snp %>% 
  filter(symbol == "ME3" & tissue == "EyeGEx")
ME3_eqtl_Brain_Spinal_cord <- eqtl_credible_snp %>% 
  filter(symbol == "ME3" & tissue == "Brain_Spinal_cord_cervical_c-1")
ME3_eqtl_Nerve_Tibial <- eqtl_credible_snp %>% 
  filter(symbol == "ME3" & tissue == "Nerve_Tibial")
ME3_eqtl_Cells_Cultured_fibroblasts <- eqtl_credible_snp %>% 
  filter(symbol == "ME3" & tissue == "Cells_Cultured_fibroblasts")
fwrite(ME3_eqtl_EyeGEx[,c(5,6)], "ME3_EyeGEx_gwas.tsv", sep = "\t")
fwrite(ME3_eqtl_EyeGEx[,c(5,3)], "ME3_EyeGEx_eqtl.tsv", sep = "\t")
fwrite(ME3_eqtl_Brain_Spinal_cord[,c(5,6)], "ME3_Brain_Spinal_cord_gwas.tsv", sep = "\t")
fwrite(ME3_eqtl_Brain_Spinal_cord[,c(5,3)], "ME3_Brain_Spinal_cord_eqtl.tsv", sep = "\t")
fwrite(ME3_eqtl_Nerve_Tibial[,c(5,6)], "ME3_Nerve_Tibial_gwas.tsv", sep = "\t")
fwrite(ME3_eqtl_Nerve_Tibial[,c(5,3)], "ME3_Nerve_Tibial_eqtl.tsv", sep = "\t")
fwrite(ME3_eqtl_Cells_Cultured_fibroblasts[,c(5,6)], "ME3_Cells_Cultured_fibroblasts_gwas.tsv", sep = "\t")
fwrite(ME3_eqtl_Cells_Cultured_fibroblasts[,c(5,3)], "ME3_Cells_Cultured_fibroblasts_eqtl.tsv", sep = "\t")
#ME3 EyeGEx
gwas_fn = 'ME3_EyeGEx_gwas.tsv'
eqtl_fn = 'ME3_EyeGEx_eqtl.tsv'
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Meta-GWAS', title2 = 'Retina eQTL',
             marker_col1 = "SNP", pval_col1 = "P",marker_col2 = "SNP", pval_col2 = "p",
             snp = "rs9667489", legend_position = "topright", combine = F)
#ME3 Brain_Spinal_cord
gwas_fn = 'ME3_Brain_Spinal_cord_gwas.tsv'
eqtl_fn = 'ME3_Brain_Spinal_cord_eqtl.tsv'
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Meta-GWAS', title2 = 'Brain_Spinal_cord_cervical_c-1 eQTL',
             marker_col1 = "SNP", pval_col1 = "P",marker_col2 = "SNP", pval_col2 = "p",
             snp = "rs9667489", legend_position = "bottomright", combine = F)
#ME3 Nerve_Tibial
gwas_fn = 'ME3_Nerve_Tibial_gwas.tsv'
eqtl_fn = 'ME3_Nerve_Tibial_eqtl.tsv'
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Meta-GWAS', title2 = 'Nerve_Tibial eQTL',
             marker_col1 = "SNP", pval_col1 = "P",marker_col2 = "SNP", pval_col2 = "p",
             snp = "rs9667489", legend_position = "topright", combine = F)
#ME3 Cells_Cultured_fibroblasts
gwas_fn = 'ME3_Cells_Cultured_fibroblasts_gwas.tsv'
eqtl_fn = 'ME3_Cells_Cultured_fibroblasts_eqtl.tsv'
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Meta-GWAS', title2 = 'Cells_Cultured_fibroblasts eQTL',
             marker_col1 = "SNP", pval_col1 = "P",marker_col2 = "SNP", pval_col2 = "p",
             snp = "rs9667489", legend_position = "bottomright", combine = F)

#PRSS23
PRSS23_eqtl_Adipose_Visceral_Omentum <- eqtl_credible_snp %>% 
  filter(symbol == "PRSS23" & tissue == "Adipose_Visceral_Omentum")
fwrite(PRSS23_eqtl_Adipose_Visceral_Omentum[,c(5,6)], "PRSS23_Adipose_Visceral_Omentum_gwas.tsv", sep = "\t")
fwrite(PRSS23_eqtl_Adipose_Visceral_Omentum[,c(5,3)], "PRSS23_Adipose_Visceral_Omentum_eqtl.tsv", sep = "\t")
PRSS23_eqtl_Artery_Aorta <- eqtl_credible_snp %>% 
  filter(symbol == "PRSS23" & tissue == "Artery_Aorta")
fwrite(PRSS23_eqtl_Artery_Aorta[,c(5,6)], "PRSS23_Artery_Aorta_gwas.tsv", sep = "\t")
fwrite(PRSS23_eqtl_Artery_Aorta[,c(5,3)], "PRSS23_Artery_Aorta_eqtl.tsv", sep = "\t")
PRSS23_eqtl_Brain_Caudate_basal_ganglia <- eqtl_credible_snp %>% 
  filter(symbol == "PRSS23" & tissue == "Brain_Caudate_basal_ganglia")
fwrite(PRSS23_eqtl_Brain_Caudate_basal_ganglia[,c(5,6)], "PRSS23_Brain_Caudate_basal_ganglia_gwas.tsv", sep = "\t")
fwrite(PRSS23_eqtl_Brain_Caudate_basal_ganglia[,c(5,3)], "PRSS23_Brain_Caudate_basal_ganglia_eqtl.tsv", sep = "\t")

#PRSS23 Adipose_Visceral_Omentum
gwas_fn = 'PRSS23_Adipose_Visceral_Omentum_gwas.tsv'
eqtl_fn = 'PRSS23_Adipose_Visceral_Omentum_eqtl.tsv'
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Meta-GWAS', title2 = 'Adipose_Visceral_Omentum eQTL',
             marker_col1 = "SNP", pval_col1 = "P",marker_col2 = "SNP", pval_col2 = "p",
             snp = "rs9667489", legend_position = "bottomright", combine = F)
#PRSS23 Artery_Aorta
gwas_fn = 'PRSS23_Artery_Aorta_gwas.tsv'
eqtl_fn = 'PRSS23_Artery_Aorta_eqtl.tsv'
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Meta-GWAS', title2 = 'Artery_Aorta eQTL',
             marker_col1 = "SNP", pval_col1 = "P",marker_col2 = "SNP", pval_col2 = "p",
             snp = "rs9667489", legend_position = "topleft", combine = F)
#PRSS23 Brain_Caudate_basal_ganglia
gwas_fn = 'PRSS23_Brain_Caudate_basal_ganglia_gwas.tsv'
eqtl_fn = 'PRSS23_Brain_Caudate_basal_ganglia_eqtl.tsv'
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Meta-GWAS', title2 = 'Brain_Caudate_basal_ganglia eQTL',
             marker_col1 = "SNP", pval_col1 = "P",marker_col2 = "SNP", pval_col2 = "p",
             snp = "rs9667489", legend_position = "topleft", combine = F)

#CCDC81
CCDC81_eqtl_EyeGEx <- eqtl_credible_snp %>% 
  filter(symbol == "CCDC81" & tissue == "EyeGEx")
fwrite(CCDC81_eqtl_EyeGEx[,c(5,6)], "CCDC81_EyeGEx_gwas.tsv", sep = "\t")
fwrite(CCDC81_eqtl_EyeGEx[,c(5,3)], "CCDC81_EyeGEx_eqtl.tsv", sep = "\t")
#CCDC81 EyeGEx
gwas_fn = 'CCDC81_EyeGEx_gwas.tsv'
eqtl_fn = 'CCDC81_EyeGEx_eqtl.tsv'
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Meta-GWAS', title2 = 'Retina eQTL',
             marker_col1 = "SNP", pval_col1 = "P",marker_col2 = "SNP", pval_col2 = "p",
             snp = "rs9667489", legend_position = "bottomright", combine = F)
