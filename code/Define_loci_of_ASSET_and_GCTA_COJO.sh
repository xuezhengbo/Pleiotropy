/share/pub/xuezb/software/plink1.9/plink --clump Meta.txt \
--bfile /share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/merge \
--clump-r2 0.4 \
--clump-kb 500 \
--clump-p1 5e-08 \
--clump-p2 5e-02 \
--threads 1 \
--memory 10000 \
--out Meta

### make loci.bed for each clump, including chr, start, and end

#bedtools merge the overlap clumps
/share/apps/bedtools2/bin/bedtools merge -i loci.bed -c 1,4 -o count,collapse

#GCTA-COJO conditional analysis
#for example
/share/pub/xuezb/software/gcta_1.93.2beta/gcta64 \
--bfile /share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/merge \
--cojo-file /share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/Meta_for_cojo.txt \
--cojo-slct \
--cojo-p 5e-8 \
--extract-region-snp rs7678123 1000 \
--out /share/pub/xuezb/biobank/GWAS/ASSET/UKB_v1/cojo/rs7678123 \
--threads 8


