##### genetic correlations across five eye diseases #########
python /share/pub/xuezb/software/GNOVA/gnova.py /share/pub/xuezb/biobank/GWAS/GNOVA/AMD_info.sumstats.gz /share/pub/xuezb/biobank/GWAS/GNOVA/DR_info.sumstats.gz \
--N1 66387 \
--N2 62229 \
--bfile /share/pub/xuezb/1000Genomes/Genotype/1000G_EUR_maf001/1000G_eur_chr@_maf001 \
--out /share/pub/xuezb/biobank/GWAS/GNOVA/AMD_DR_1000G_maf001.txt
python /share/pub/xuezb/software/GNOVA/gnova.py /share/pub/xuezb/biobank/GWAS/GNOVA/AMD_info.sumstats.gz /share/pub/xuezb/biobank/GWAS/GNOVA/GLC_info.sumstats.gz \
--N1 66387 \
--N2 68390 \
--bfile /share/pub/xuezb/1000Genomes/Genotype/1000G_EUR_maf001/1000G_eur_chr@_maf001 \
--out /share/pub/xuezb/biobank/GWAS/GNOVA/AMD_GLC_1000G_maf001.txt
python /share/pub/xuezb/software/GNOVA/gnova.py /share/pub/xuezb/biobank/GWAS/GNOVA/AMD_info.sumstats.gz /share/pub/xuezb/biobank/GWAS/GNOVA/RD_info.sumstats.gz \
--N1 66387 \
--N2 64011 \
--bfile /share/pub/xuezb/1000Genomes/Genotype/1000G_EUR_maf001/1000G_eur_chr@_maf001 \
--out /share/pub/xuezb/biobank/GWAS/GNOVA/AMD_RD_1000G_maf001.txt
python /share/pub/xuezb/software/GNOVA/gnova.py /share/pub/xuezb/biobank/GWAS/GNOVA/AMD_info.sumstats.gz /share/pub/xuezb/biobank/GWAS/GNOVA/Myopia_info.sumstats.gz \
--N1 66387 \
--N2 64268 \
--bfile /share/pub/xuezb/1000Genomes/Genotype/1000G_EUR_maf001/1000G_eur_chr@_maf001 \
--out /share/pub/xuezb/biobank/GWAS/GNOVA/AMD_Myopia_1000G_maf001.txt
python /share/pub/xuezb/software/GNOVA/gnova.py /share/pub/xuezb/biobank/GWAS/GNOVA/DR_info.sumstats.gz /share/pub/xuezb/biobank/GWAS/GNOVA/GLC_info.sumstats.gz \
--N1 62229 \
--N2 68390 \
--bfile /share/pub/xuezb/1000Genomes/Genotype/1000G_EUR_maf001/1000G_eur_chr@_maf001 \
--out /share/pub/xuezb/biobank/GWAS/GNOVA/DR_GLC_1000G_maf001.txt
python /share/pub/xuezb/software/GNOVA/gnova.py /share/pub/xuezb/biobank/GWAS/GNOVA/DR_info.sumstats.gz /share/pub/xuezb/biobank/GWAS/GNOVA/RD_info.sumstats.gz \
--N1 62229 \
--N2 64011 \
--bfile /share/pub/xuezb/1000Genomes/Genotype/1000G_EUR_maf001/1000G_eur_chr@_maf001 \
--out /share/pub/xuezb/biobank/GWAS/GNOVA/DR_RD_1000G_maf001.txt
python /share/pub/xuezb/software/GNOVA/gnova.py /share/pub/xuezb/biobank/GWAS/GNOVA/DR_info.sumstats.gz /share/pub/xuezb/biobank/GWAS/GNOVA/Myopia_info.sumstats.gz \
--N1 62229 \
--N2 64268 \
--bfile /share/pub/xuezb/1000Genomes/Genotype/1000G_EUR_maf001/1000G_eur_chr@_maf001 \
--out /share/pub/xuezb/biobank/GWAS/GNOVA/DR_Myopia_1000G_maf001.txt
python /share/pub/xuezb/software/GNOVA/gnova.py /share/pub/xuezb/biobank/GWAS/GNOVA/GLC_info.sumstats.gz /share/pub/xuezb/biobank/GWAS/GNOVA/RD_info.sumstats.gz \
--N1 68390 \
--N2 64011 \
--bfile /share/pub/xuezb/1000Genomes/Genotype/1000G_EUR_maf001/1000G_eur_chr@_maf001 \
--out /share/pub/xuezb/biobank/GWAS/GNOVA/GLC_RD_1000G_maf001.txt
python /share/pub/xuezb/software/GNOVA/gnova.py /share/pub/xuezb/biobank/GWAS/GNOVA/GLC_info.sumstats.gz /share/pub/xuezb/biobank/GWAS/GNOVA/Myopia_info.sumstats.gz \
--N1 68390 \
--N2 64268 \
--bfile /share/pub/xuezb/1000Genomes/Genotype/1000G_EUR_maf001/1000G_eur_chr@_maf001 \
--out /share/pub/xuezb/biobank/GWAS/GNOVA/GLC_Myopia_1000G_maf001.txt
python /share/pub/xuezb/software/GNOVA/gnova.py /share/pub/xuezb/biobank/GWAS/GNOVA/RD_info.sumstats.gz /share/pub/xuezb/biobank/GWAS/GNOVA/Myopia_info.sumstats.gz \
--N1 64011 \
--N2 64268 \
--bfile /share/pub/xuezb/1000Genomes/Genotype/1000G_EUR_maf001/1000G_eur_chr@_maf001 \
--out /share/pub/xuezb/biobank/GWAS/GNOVA/RD_Myopia_1000G_maf001.txt

################ genetic correlations between eye diseases and risk traits ########################
#prepare sumstats.gz
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/alcohol_Karlsson_noMHC.txt --frq EAF_A1 --maf-min 0.05 --N 414343 --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out alcohol_Karlsson_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/alcohol_liu_noMHC.txt --a1 ALT --a2 REF --frq AF --maf-min 0.05 --N-col EFFECTIVE_N --ignore N --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out alcohol_liu_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/BMI_Watanabe_noMHC.txt --frq MAF --maf-min 0.05 --info INFO --N 385336 --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out BMI_Watanabe_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/BMI_Wood2016_noMHC.txt --frq Freq --maf-min 0.05 --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out BMI_Wood_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/DBP_withID_noMHC.txt --frq Freq1 --maf-min 0.05 --N-col N_effective --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out DBP_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/SBP_noMHC.txt --frq Freq1 --maf-min 0.05 --N-col N_effective --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out SBP_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/EduYears_Okbay_noMHC.txt --frq freq --maf-min 0.05 --N 405072 --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out EduYears_Okbay_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/HbA1c_noMHC.txt --frq freq --maf-min 0.05 --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out HbA1c_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/Height_Watanabe_noMHC.txt --frq MAF --maf-min 0.05 --info INFO --N 385748 --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out Height_Watanabe_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/highbloodpressure_ZhuZ_noMHC.txt --frq MAF --info INFO --maf-min 0.05 --N-cas 144793 --N-con 313761 --ignore N --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out HBP_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/Intelligence_Watanabe_noMHC.txt --frq MAF --maf-min 0.05 --info INFO --N 125935 --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out Intelligence_Watanabe_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/jointGwasMc_HDL_noMHC.txt --frq Freq1 --maf-min 0.05 --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out HDL_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/jointGwasMc_LDL_noMHC.txt --frq Freq1 --maf-min 0.05 --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out LDL_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/jointGwasMc_TC_noMHC.txt --frq Freq1 --maf-min 0.05 --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out TC_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/jointGwasMc_TG_noMHC.txt --frq Freq1 --maf-min 0.05 --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out TG_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/Klimentidis_AccAve_noMHC.txt --a2 ALLELE0 --frq A1FREQ --maf-min 0.05 --info INFO --N 91084 --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out Klimentidis_AccAve_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/Klimentidis_MVPA_noMHC.txt --a2 ALLELE0 --frq A1FREQ --maf-min 0.05 --info INFO --N 377234 --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out Klimentidis_MVPA_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/Klimentidis_SSOE_noMHC.txt --a2 ALLELE0 --frq A1FREQ --maf-min 0.05 --info INFO --N 350492 --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out Klimentidis_SSOE_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/Klimentidis_VPA_noMHC.txt --a2 ALLELE0 --frq A1FREQ --maf-min 0.05 --info INFO --N 261055 --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out Klimentidis_VPA_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/Lane_Insomnia_noMHC.txt --frq FREQ1 --maf-min 0.05 --info INFO --N-cas 345022 --N-con 108357 --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out Insomnia_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/Dashti_sleepduration_noMHC.txt --frq FREQ1 --maf-min 0.05 --info INFO --N 446118 --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out sleepduration_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/MAGIC_Manning_FastingGlucose_MainEffect.txt --frq maf --maf-min 0.05 --N 58074 --p BMIadjMainP --signed-sumstats BMIadjMainEffects,0 --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out FG_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/MAGIC_Manning_lnFastingInsulin_MainEffect.txt --frq maf --maf-min 0.05 --N 51750 --p BMIadjMainP --signed-sumstats BMIadjMainEffects,0 --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out FI_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/MAGNETIC_Glucose_noMHC.txt --frq freq --maf-min 0.05 --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out MAGNETIC_Glucose_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/mcv_noMHC.txt --frq EAF --maf-min 0.05 --info INFO --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out mcv_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/mpv_noMHC.txt --frq EAF --maf-min 0.05 --info INFO --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out mpv_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/plt_noMHC.txt --frq EAF --maf-min 0.05 --info INFO --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out plt_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/moderate_activity_Doherty_noMHC.txt --a2 ALLELE0 --frq A1FREQ --maf-min 0.05 --info INFO --p P_LINREG --N 91105 --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out moderate_activity_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/overall_activity_Doherty_noMHC.txt --a2 ALLELE0 --frq A1FREQ --maf-min 0.05 --info INFO --p P_LINREG --N 91105 --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out overall_activity_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/sedentary_Doherty_noMHC.txt --a2 ALLELE0 --frq A1FREQ --maf-min 0.05 --info INFO --p P_LINREG --N 91105 --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out sedentary_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/walking_Doherty_noMHC.txt --a2 ALLELE0 --frq A1FREQ --maf-min 0.05 --info INFO --p P_LINREG --N 91105 --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out walking_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/T2D_Wood2016_noMHC.txt --frq Freq --maf-min 0.05 --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out T2D_Wood_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/T2D_Xue2018_noMHC.txt --frq frq_A1 --maf-min 0.05 --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out T2D_Xue_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/tag.cpd_noMHC.txt --frq FRQ_A --ignore FRQ_U --maf-min 0.05 --info INFO --N 38181 --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out tag_cpd_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/tag.evrsmk_noMHC.txt --frq FRQ_A --ignore FRQ_U --maf-min 0.05 --info INFO --N-cas 39022 --N-con 30387 --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out tag_evrsmk_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/tag.former_noMHC.txt --frq FRQ_A --ignore FRQ_U --maf-min 0.05 --info INFO --N-cas 20619 --N-con 15226 --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out tag_former_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/tag.logonset_noMHC.txt --frq FRQ_A --ignore FRQ_U --maf-min 0.05 --info INFO --N 22438 --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out tag_logonset_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/TelevisionTime_Watanabe_noMHC.txt --frq MAF --maf-min 0.05 --info INFO --N 365236 --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out TelevisionTime_Watanabe_maf005
python /share/pub/xuezb/software/ldsc/munge_sumstats.py --sumstats /share/pub/xuezb/GWAS_other_traits/VitD_Jiang_noMHC.txt --chunksize 500000 --merge-alleles /share/pub/xuezb/software/ldsc/w_hm3.snplist --out VitD_Jiang

#estimate genetic correlation
python /share/pub/xuezb/software/ldsc/ldsc.py --rg AMD_info_maf005.sumstats.gz,SBP_maf005.sumstats.gz,DBP_maf005.sumstats.gz,HBP_maf005.sumstats.gz,MAGNETIC_Glucose_maf005.sumstats.gz,FG_maf005.sumstats.gz,FI_maf005.sumstats.gz,HbA1c_maf005.sumstats.gz,T2D_Wood_maf005.sumstats.gz,T2D_Xue_maf005.sumstats.gz,HDL_maf005.sumstats.gz,LDL_maf005.sumstats.gz,TC_maf005.sumstats.gz,TG_maf005.sumstats.gz,mcv_maf005.sumstats.gz,mpv_maf005.sumstats.gz,plt_maf005.sumstats.gz,VitD_Jiang.sumstats.gz,BMI_Watanabe_maf005.sumstats.gz,BMI_Wood_maf005.sumstats.gz,Height_Watanabe_maf005.sumstats.gz,Intelligence_Watanabe_maf005.sumstats.gz,Insomnia_maf005.sumstats.gz,sleepduration_maf005.sumstats.gz,alcohol_Karlsson_maf005.sumstats.gz,alcohol_liu_maf005.sumstats.gz,tag_cpd_maf005.sumstats.gz,tag_evrsmk_maf005.sumstats.gz,tag_former_maf005.sumstats.gz,tag_logonset_maf005.sumstats.gz,EduYears_Okbay_maf005.sumstats.gz,TelevisionTime_Watanabe_maf005.sumstats.gz,Klimentidis_AccAve_maf005.sumstats.gz,Klimentidis_MVPA_maf005.sumstats.gz,Klimentidis_SSOE_maf005.sumstats.gz,Klimentidis_VPA_maf005.sumstats.gz,moderate_activity_maf005.sumstats.gz,overall_activity_maf005.sumstats.gz,sedentary_maf005.sumstats.gz,walking_maf005.sumstats.gz --ref-ld-chr /share/pub/xuezb/software/ldsc/eur_w_ld_chr/ --w-ld-chr /share/pub/xuezb/software/ldsc/eur_w_ld_chr/ --samp-prev 0.0885,nan,nan,0.3158,nan,nan,nan,nan,0.0343,0.0954,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,0.761,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan --pop-prev 0.123,nan,nan,0.311,nan,nan,nan,nan,0.063,0.063,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,0.33,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan --out /share/pub/xuezb/biobank/GWAS/LDSC/rg/info_maf005/AMD_other_traits
python /share/pub/xuezb/software/ldsc/ldsc.py --rg DR_info_maf005.sumstats.gz,SBP_maf005.sumstats.gz,DBP_maf005.sumstats.gz,HBP_maf005.sumstats.gz,MAGNETIC_Glucose_maf005.sumstats.gz,FG_maf005.sumstats.gz,FI_maf005.sumstats.gz,HbA1c_maf005.sumstats.gz,T2D_Wood_maf005.sumstats.gz,T2D_Xue_maf005.sumstats.gz,HDL_maf005.sumstats.gz,LDL_maf005.sumstats.gz,TC_maf005.sumstats.gz,TG_maf005.sumstats.gz,mcv_maf005.sumstats.gz,mpv_maf005.sumstats.gz,plt_maf005.sumstats.gz,VitD_Jiang.sumstats.gz,BMI_Watanabe_maf005.sumstats.gz,BMI_Wood_maf005.sumstats.gz,Height_Watanabe_maf005.sumstats.gz,Intelligence_Watanabe_maf005.sumstats.gz,Insomnia_maf005.sumstats.gz,sleepduration_maf005.sumstats.gz,alcohol_Karlsson_maf005.sumstats.gz,alcohol_liu_maf005.sumstats.gz,tag_cpd_maf005.sumstats.gz,tag_evrsmk_maf005.sumstats.gz,tag_former_maf005.sumstats.gz,tag_logonset_maf005.sumstats.gz,EduYears_Okbay_maf005.sumstats.gz,TelevisionTime_Watanabe_maf005.sumstats.gz,Klimentidis_AccAve_maf005.sumstats.gz,Klimentidis_MVPA_maf005.sumstats.gz,Klimentidis_SSOE_maf005.sumstats.gz,Klimentidis_VPA_maf005.sumstats.gz,moderate_activity_maf005.sumstats.gz,overall_activity_maf005.sumstats.gz,sedentary_maf005.sumstats.gz,walking_maf005.sumstats.gz --ref-ld-chr /share/pub/xuezb/software/ldsc/eur_w_ld_chr/ --w-ld-chr /share/pub/xuezb/software/ldsc/eur_w_ld_chr/ --samp-prev 0.0265,nan,nan,0.3158,nan,nan,nan,nan,0.0343,0.0954,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,0.761,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan --pop-prev 0.0176,nan,nan,0.311,nan,nan,nan,nan,0.063,0.063,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,0.33,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan --out /share/pub/xuezb/biobank/GWAS/LDSC/rg/info_maf005/DR_other_traits
python /share/pub/xuezb/software/ldsc/ldsc.py --rg GLC_info_maf005.sumstats.gz,SBP_maf005.sumstats.gz,DBP_maf005.sumstats.gz,HBP_maf005.sumstats.gz,MAGNETIC_Glucose_maf005.sumstats.gz,FG_maf005.sumstats.gz,FI_maf005.sumstats.gz,HbA1c_maf005.sumstats.gz,T2D_Wood_maf005.sumstats.gz,T2D_Xue_maf005.sumstats.gz,HDL_maf005.sumstats.gz,LDL_maf005.sumstats.gz,TC_maf005.sumstats.gz,TG_maf005.sumstats.gz,mcv_maf005.sumstats.gz,mpv_maf005.sumstats.gz,plt_maf005.sumstats.gz,VitD_Jiang.sumstats.gz,BMI_Watanabe_maf005.sumstats.gz,BMI_Wood_maf005.sumstats.gz,Height_Watanabe_maf005.sumstats.gz,Intelligence_Watanabe_maf005.sumstats.gz,Insomnia_maf005.sumstats.gz,sleepduration_maf005.sumstats.gz,alcohol_Karlsson_maf005.sumstats.gz,alcohol_liu_maf005.sumstats.gz,tag_cpd_maf005.sumstats.gz,tag_evrsmk_maf005.sumstats.gz,tag_former_maf005.sumstats.gz,tag_logonset_maf005.sumstats.gz,EduYears_Okbay_maf005.sumstats.gz,TelevisionTime_Watanabe_maf005.sumstats.gz,Klimentidis_AccAve_maf005.sumstats.gz,Klimentidis_MVPA_maf005.sumstats.gz,Klimentidis_SSOE_maf005.sumstats.gz,Klimentidis_VPA_maf005.sumstats.gz,moderate_activity_maf005.sumstats.gz,overall_activity_maf005.sumstats.gz,sedentary_maf005.sumstats.gz,walking_maf005.sumstats.gz --ref-ld-chr /share/pub/xuezb/software/ldsc/eur_w_ld_chr/ --w-ld-chr /share/pub/xuezb/software/ldsc/eur_w_ld_chr/ --samp-prev 0.1151,nan,nan,0.3158,nan,nan,nan,nan,0.0343,0.0954,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,0.761,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan --pop-prev 0.0251,nan,nan,0.311,nan,nan,nan,nan,0.063,0.063,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,0.33,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan --out /share/pub/xuezb/biobank/GWAS/LDSC/rg/info_maf005/GLC_other_traits
python /share/pub/xuezb/software/ldsc/ldsc.py --rg RD_info_maf005.sumstats.gz,SBP_maf005.sumstats.gz,DBP_maf005.sumstats.gz,HBP_maf005.sumstats.gz,MAGNETIC_Glucose_maf005.sumstats.gz,FG_maf005.sumstats.gz,FI_maf005.sumstats.gz,HbA1c_maf005.sumstats.gz,T2D_Wood_maf005.sumstats.gz,T2D_Xue_maf005.sumstats.gz,HDL_maf005.sumstats.gz,LDL_maf005.sumstats.gz,TC_maf005.sumstats.gz,TG_maf005.sumstats.gz,mcv_maf005.sumstats.gz,mpv_maf005.sumstats.gz,plt_maf005.sumstats.gz,VitD_Jiang.sumstats.gz,BMI_Watanabe_maf005.sumstats.gz,BMI_Wood_maf005.sumstats.gz,Height_Watanabe_maf005.sumstats.gz,Intelligence_Watanabe_maf005.sumstats.gz,Insomnia_maf005.sumstats.gz,sleepduration_maf005.sumstats.gz,alcohol_Karlsson_maf005.sumstats.gz,alcohol_liu_maf005.sumstats.gz,tag_cpd_maf005.sumstats.gz,tag_evrsmk_maf005.sumstats.gz,tag_former_maf005.sumstats.gz,tag_logonset_maf005.sumstats.gz,EduYears_Okbay_maf005.sumstats.gz,TelevisionTime_Watanabe_maf005.sumstats.gz,Klimentidis_AccAve_maf005.sumstats.gz,Klimentidis_MVPA_maf005.sumstats.gz,Klimentidis_SSOE_maf005.sumstats.gz,Klimentidis_VPA_maf005.sumstats.gz,moderate_activity_maf005.sumstats.gz,overall_activity_maf005.sumstats.gz,sedentary_maf005.sumstats.gz,walking_maf005.sumstats.gz --ref-ld-chr /share/pub/xuezb/software/ldsc/eur_w_ld_chr/ --w-ld-chr /share/pub/xuezb/software/ldsc/eur_w_ld_chr/ --samp-prev 0.0539,nan,nan,0.3158,nan,nan,nan,nan,0.0343,0.0954,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,0.761,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan --pop-prev 0.0002,nan,nan,0.311,nan,nan,nan,nan,0.063,0.063,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,0.33,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan --out /share/pub/xuezb/biobank/GWAS/LDSC/rg/info_maf005/RD_other_traits
python /share/pub/xuezb/software/ldsc/ldsc.py --rg Myopia_info_maf005.sumstats.gz,SBP_maf005.sumstats.gz,DBP_maf005.sumstats.gz,HBP_maf005.sumstats.gz,MAGNETIC_Glucose_maf005.sumstats.gz,FG_maf005.sumstats.gz,FI_maf005.sumstats.gz,HbA1c_maf005.sumstats.gz,T2D_Wood_maf005.sumstats.gz,T2D_Xue_maf005.sumstats.gz,HDL_maf005.sumstats.gz,LDL_maf005.sumstats.gz,TC_maf005.sumstats.gz,TG_maf005.sumstats.gz,mcv_maf005.sumstats.gz,mpv_maf005.sumstats.gz,plt_maf005.sumstats.gz,VitD_Jiang.sumstats.gz,BMI_Watanabe_maf005.sumstats.gz,BMI_Wood_maf005.sumstats.gz,Height_Watanabe_maf005.sumstats.gz,Intelligence_Watanabe_maf005.sumstats.gz,Insomnia_maf005.sumstats.gz,sleepduration_maf005.sumstats.gz,alcohol_Karlsson_maf005.sumstats.gz,alcohol_liu_maf005.sumstats.gz,tag_cpd_maf005.sumstats.gz,tag_evrsmk_maf005.sumstats.gz,tag_former_maf005.sumstats.gz,tag_logonset_maf005.sumstats.gz,EduYears_Okbay_maf005.sumstats.gz,TelevisionTime_Watanabe_maf005.sumstats.gz,Klimentidis_AccAve_maf005.sumstats.gz,Klimentidis_MVPA_maf005.sumstats.gz,Klimentidis_SSOE_maf005.sumstats.gz,Klimentidis_VPA_maf005.sumstats.gz,moderate_activity_maf005.sumstats.gz,overall_activity_maf005.sumstats.gz,sedentary_maf005.sumstats.gz,walking_maf005.sumstats.gz --ref-ld-chr /share/pub/xuezb/software/ldsc/eur_w_ld_chr/ --w-ld-chr /share/pub/xuezb/software/ldsc/eur_w_ld_chr/ --samp-prev 0.4356,nan,nan,0.3158,nan,nan,nan,nan,0.0343,0.0954,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,0.761,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan --pop-prev 0.306,nan,nan,0.311,nan,nan,nan,nan,0.063,0.063,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,0.33,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan --out /share/pub/xuezb/biobank/GWAS/LDSC/rg/info_maf005/Myopia_other_traits






