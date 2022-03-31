#!/bin/bash
#PBS -N AMD_gwas
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=4
#PBS -j oe
cd /share/pub/xuezb/biobank/GWAS/AMD/data_68304
/share/pub/xuezb/software/plink2/plink2 --bfile merge_maf001_no_relative_v2 --covar covar_v2.txt --covar-variance-standardize --glm hide-covar --adjust --out ./result_v2/AMD_10pc --threads 4
