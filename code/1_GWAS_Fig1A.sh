#select variants with INFO > 0.4
for CHR in {1..22}
do
  value1="/.../ukb_mfi_chr"${CHR}"_v3.txt"
  value2="/.../chr"${CHR}"_info04.txt"
  awk '{ if ($8 >= 0.4) print $2 }' $value1 > $value2
done
awk '{ if ($8 >= 0.4) print $2 }' /.../ukb_mfi_chrX_v3.txt > /.../chrX_info04.txt

#extract samples and variants, and QC
for CHR in {1..22}
do
  value1="/.../ukb_imp_chr"${CHR}"_v3.bgen"
  value2="/.../ukb45270_imp_chr"${CHR}"_v3_s487296.sample"
  value3="/.../chr"${CHR}"_info04.txt"
  value4="chr"${CHR}
  plink2 --bgen $value1 ref-first --sample $value2 --extract $value3 --keep myopia_list.txt --geno 0.05 --hwe 1e-6 --maf 0.01 --make-bed --threads 4 --out $value4
done
plink2 --bgen /.../ukb_imp_chrX_v3.bgen ref-first --sample /.../ukb45270_imp_chrX_v3_s486645.sample --extract /.../chrX_info04.txt --keep myopia_list.txt --geno 0.05 --hwe 1e-6 --maf 0.01 --make-bed --threads 4 --out chrX
#merge chromosome
for CHR in {1..22}
do
  value1="chr"${CHR}
  value2="tmp_chr"${CHR}
  plink2 --bfile $value1 --rm-dup force-first --make-bed --threads 4 --out $value2
done
plink2 --bfile chrX --rm-dup force-first --make-bed --threads 4 --out tmp_chrX
./plink1.9/plink --merge-list batch.txt --make-bed --threads 4 --out merge
rm tmp*
plink2 --bfile merge --maf 0.01 --make-bed --out merge_maf001 --threads 4 --memory 50000

#estimate kinship with pihat > 0.2
plink2 --bfile merge_maf001 --indep-pairwise 50 5 0.2 --out indepSNP
./plink1.9/plink --bfile merge_maf001 --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2
plink2 --bfile merge_maf001 --missing --threads 2 --memory 50000
#remove the sample with lower call rate in a pair of kinship
plink2 --bfile merge_maf001 --remove pihat_low_call_rate.txt --geno 0.05 --hwe 1e-6 --maf 0.01 --make-bed --threads 2 --memory 50000 --out merge_maf001_no_relative

#PCA
plink2 --bfile merge_maf001_no_relative --pca 10 --out merge_maf001_no_relative_pca
#make covar.txt with age and 10 PCs
#GWAS logistic model
plink2 --bfile merge_maf001_no_relative --covar covar.txt --covar-variance-standardize --glm hide-covar --adjust --out ./result/AMD_10pc --threads 4



