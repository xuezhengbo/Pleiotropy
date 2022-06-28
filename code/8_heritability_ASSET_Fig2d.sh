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
