#METASOFT
#input file contains the Beta and SE of the SNPs in individual GWASs (AMD, DR, GLC, RD, and myopia)
java -jar /share/pub/xuezb/software/Metasoft/Metasoft.jar \
-input input_correct_allele.txt \
-pvalue_table /share/pub/xuezb/software/Metasoft/HanEskinPvalueTable.txt \
-binary_effects \
-mvalue \
-output asset_meta_correct_allele.out


#visualize by ForestPMPlot
python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
1:164194417_TA_T \
- \
1:164194417_TA_T.pdf

python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
rs41393947 \
EFEMP1 \
rs41393947.pdf

python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
2:146913702_CCTCT_C \
ARHGAP15 \
2:146913702_CCTCT_C.pdf

python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
rs62181740 \
DYNC1I2,METTL8,AC068039.1,SLC25A12,DCAF17,METAP1D,HAT1,DLX1,DLX2,ITGA6,PDK1 \
rs62181740.pdf

python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
rs7678123 \
FGF5,C4orf22,BMP3,PRKG2,PAQR3 \
rs7678123.pdf

python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
rs13118211 \
RASGEF1B \
rs13118211.pdf

python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
rs10036789 \
ZNF366,MAP1B,MRPS27,PTCD2,TNPO1 \
rs10036789.pdf

python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
rs7744813 \
KCNQ5 \
rs7744813.pdf

python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
rs138650617 \
BICC1,UBE2D1,TFAM,PHYHIPL,CCDC6 \
rs138650617.pdf

python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
rs10887262 \
LRIT1,RGR \
rs10887262.pdf

python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
rs12570944 \
LRRC27,STK32C,PWWP2B \
rs12570944.pdf

python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
rs9667489 \
ME3 \
rs9667489.pdf

python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
rs2738265 \
BMP4 \
rs2738265.pdf

python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
rs112115087 \
C14orf39,DHRS7,HIF1A,PCNXL4,SIX1,SIX4,SIX6,SLC38A6,TMEM30B,TRMT5 \
rs112115087.pdf

python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
rs12950511 \
ABI3,GNGT2,PHOSPHO1,ZNF652,RP11-81K2.1 \
rs12950511.pdf

python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
rs9330814 \
WNT7B \
rs9330814.pdf
python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
rs10529326 \
ALDH9A1,TMCO1,UCK2,MGST3 \
rs10529326.pdf

python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
rs41393947 \
EFEMP1 \
rs41393947.pdf

python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
2:146913702_CCTCT_C \
ARHGAP15 \
2:146913702_CCTCT_C.pdf

python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
rs2573222 \
PRSS56,CHRND,CHRNG,TIGD1,EIF4E2 \
rs2573222.pdf

python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
rs773674648 \
CADM2 \
rs773674648.pdf

python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
rs12193446 \
LAMA2,PTPRK,L3MBTL3 \
rs12193446.pdf

python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
rs72644322 \
ZMAT4 \
rs72644322.pdf

python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
rs10089517 \
TOX \
rs10089517.pdf

python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
9:22053956_TA_T \
RP11-145E5.5,CDKN2A,CDKN2B,MTAP \
9:22053956_TA_T.pdf

python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
rs2472493 \
ABCA1,FSD1L \
rs2472493.pdf

python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
rs10824539 \
KCNMA1 \
rs10824539.pdf

python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
rs1367024 \
FRMPD2,AGAP10,MAPK8,NPY4R \
rs1367024.pdf

python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
rs36090025 \
TCF7L2\
rs36090025.pdf

python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
rs3750846 \
ARMS2,HTRA1,BTBD16,PLEKHA1 \
rs3750846.pdf

python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
rs11606250 \
LRRC4C \
rs11606250.pdf

python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
rs5442 \
GNB3,CDCA3,USP5,SPSB2,LRRC23,ATN1,C12orf57,PTPN6,LPAR5,COPS7A,LAG3,ENO2,C1R,C1RL,RBP5,CLSTN3 \
rs5442.pdf

python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
rs3138142 \
ITGA7,RP11-644F5.10,BLOC1S1,RDH5,CD63,SARNP,RP11-762I7.5,ORMDL2,DNAJC14,PMEL,SUOX,RPS26,GDF11,WIBG,DGKA,RP11-603J24.9,PA2G4,RPL41,ZC3H10,ESYT1 \
rs3138142.pdf

python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
rs1347190 \
ZIC2,PCCA \
rs1347190.pdf

python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
rs524952 \
GOLGA8A,GJD2,ACTC1 \
rs524952.pdf

python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
rs11631829 \
RASGRF1 \
rs11631829.pdf

python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
rs7184522 \
RBFOX1 \
rs7184522.pdf

python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
rs113941606 \
SHISA6,DNAH9 \
rs113941606.pdf

python pmplot.py \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/input.txt \
/share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/metasoft/asset_meta_correct_allele.out \
study.name.txt \
study.order.txt \
rs7405453 \
C17orf70,NPLOC4,PDE6G,OXLD1,CCDC137,ARL16,HGS,MRPL12,SLC25A10,ACTG1,FSCN2,CHMP6,AC127496.1,AATK,ENTHD2,C17orf89,SLC38A10,RP11-1055B8.6,RP11-1055B8.7,ARHGDIA,ALYREF,ANAPC11,SIRT7,NOTUM,ASPSCR1,RAC3,DCXR,CSNK1D,C17orf62 \
rs7405453.pdf
