########### AMD #############################
python ./ldsc/munge_sumstats.py \
--sumstats ./AMD_10pc_info_noMHC.txt \
--frq Freq1 \
--info INFO \
--N-cas 5873 \
--N-con 60514 \
--ignore Beta \
--chunksize 500000 \
--merge-alleles ./ldsc/w_hm3.snplist \
--out AMD_info

python ./ldsc/ldsc.py \
--h2 AMD_info.sumstats.gz \
--samp-prev 0.0885 \
--pop-prev 0.123 \
--ref-ld-chr ./ldsc/eur_w_ld_chr/ \
--w-ld-chr  ./ldsc/eur_w_ld_chr/ \
--out ./h2/AMD_info

############### DR ############################
python ./ldsc/munge_sumstats.py \
--sumstats ./DR_10pc_info_noMHC.txt \
--frq Freq1 \
--info INFO \
--N-cas 1652 \
--N-con 60577 \
--ignore Beta \
--chunksize 500000 \
--merge-alleles ./ldsc/w_hm3.snplist \
--out DR_info

python ./ldsc/ldsc.py \
--h2 DR_info.sumstats.gz \
--samp-prev 0.0265 \
--pop-prev 0.0176 \
--ref-ld-chr ./ldsc/eur_w_ld_chr/ \
--w-ld-chr  ./ldsc/eur_w_ld_chr/ \
--out ./h2/DR_info

############ GLC ##############################
python ./ldsc/munge_sumstats.py \
--sumstats ./GLC_10pc_info_noMHC.txt \
--frq Freq1 \
--info INFO \
--N-cas 7873 \
--N-con 60517 \
--ignore Beta \
--chunksize 500000 \
--merge-alleles ./ldsc/w_hm3.snplist \
--out GLC_info

python ./ldsc/ldsc.py \
--h2 GLC_info.sumstats.gz \
--samp-prev 0.1151 \
--pop-prev 0.0251 \
--ref-ld-chr ./ldsc/eur_w_ld_chr/ \
--w-ld-chr  ./ldsc/eur_w_ld_chr/ \
--out ./h2/GLC_info

########### RD #################################
python ./ldsc/munge_sumstats.py \
--sumstats RD_10pc_info_noMHC.txt \
--frq Freq1 \
--info INFO \
--N-cas 3449 \
--N-con 60562 \
--ignore Beta \
--chunksize 500000 \
--merge-alleles ./ldsc/w_hm3.snplist \
--out RD_info

python ./ldsc/ldsc.py \
--h2 RD_info.sumstats.gz \
--samp-prev 0.0539 \
--pop-prev 0.0002 \
--ref-ld-chr ./ldsc/eur_w_ld_chr/ \
--w-ld-chr  ./ldsc/eur_w_ld_chr/ \
--out ./h2/RD_info

############### Myopia #####################
python /ldsc/munge_sumstats.py \
--sumstats ./myopia_10pc_info_noMHC.txt \
--frq Freq1 \
--info INFO \
--N-cas 27993 \
--N-con 36275 \
--chunksize 500000 \
--merge-alleles /ldsc/w_hm3.snplist \
--out Myopia_info

python /ldsc/ldsc.py \
--h2 Myopia_info.sumstats.gz \
--samp-prev 0.4356 \
--pop-prev 0.306 \
--ref-ld-chr /ldsc/eur_w_ld_chr/ \
--w-ld-chr  /ldsc/eur_w_ld_chr/ \
--out ./h2/Myopia_info
###############################################

########### heritability of ASSET result###########
python ./ldsc/munge_sumstats.py \
--sumstats ./Meta_info_noMHC.txt \
--N-cas 43877 \
--N-con 44373 \
--frq freq \
--info INFO \
--chunksize 500000 \
--merge-alleles ./ldsc/w_hm3.snplist \
--out ASSET_info

python ./ldsc/ldsc.py \
--h2 ASSET_info.sumstats.gz \
--samp-prev 0.497 \
--pop-prev 0.497 \
--ref-ld-chr ./ldsc/eur_w_ld_chr/ \
--w-ld-chr  ./ldsc/eur_w_ld_chr/ \
--out ./h2/ASSET_info

############## partition heritability of ASSET ###################
python ./ldsc/ldsc.py \
--h2 ASSET_info.sumstats.gz \
--ref-ld-chr ./ldsc/baseline_v1.2/baseline. \
--w-ld-chr ./ldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--overlap-annot \
--frqfile-chr ./ldsc/1000G_Phase3_frq/1000G.EUR.QC. \
--out ./h2/ASSET_info_baseline




