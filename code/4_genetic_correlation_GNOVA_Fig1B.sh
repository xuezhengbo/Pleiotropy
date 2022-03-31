python ./GNOVA/gnova.py ./AMD_info.sumstats.gz ./DR_info.sumstats.gz \
--N1 66387 \
--N2 62229 \
--bfile ./1000G_EUR_maf001/1000G_eur_chr@_maf001 \
--out ./GNOVA/AMD_DR_1000G_maf001.txt

python ./GNOVA/gnova.py ./AMD_info.sumstats.gz ./GLC_info.sumstats.gz \
--N1 66387 \
--N2 68390 \
--bfile ./1000G_EUR_maf001/1000G_eur_chr@_maf001 \
--out ./GNOVA/AMD_GLC_1000G_maf001.txt

python ./GNOVA/gnova.py ./AMD_info.sumstats.gz ./RD_info.sumstats.gz \
--N1 66387 \
--N2 64011 \
--bfile ./1000G_EUR_maf001/1000G_eur_chr@_maf001 \
--out ./GNOVA/AMD_RD_1000G_maf001.txt

python ./GNOVA/gnova.py ./AMD_info.sumstats.gz ./Myopia_info.sumstats.gz \
--N1 66387 \
--N2 64268 \
--bfile ./1000G_EUR_maf001/1000G_eur_chr@_maf001 \
--out /./GNOVA/AMD_Myopia_1000G_maf001.txt

python ./GNOVA/gnova.py ./DR_info.sumstats.gz ./GLC_info.sumstats.gz \
--N1 62229 \
--N2 68390 \
--bfile ./1000G_EUR_maf001/1000G_eur_chr@_maf001 \
--out ./GNOVA/DR_GLC_1000G_maf001.txt

python ./GNOVA/gnova.py ./DR_info.sumstats.gz ./RD_info.sumstats.gz \
--N1 62229 \
--N2 64011 \
--bfile ./1000G_EUR_maf001/1000G_eur_chr@_maf001 \
--out ./GNOVA/DR_RD_1000G_maf001.txt

python ./GNOVA/gnova.py ./DR_info.sumstats.gz ./Myopia_info.sumstats.gz \
--N1 62229 \
--N2 64268 \
--bfile ./1000G_EUR_maf001/1000G_eur_chr@_maf001 \
--out ./GNOVA/DR_Myopia_1000G_maf001.txt

python ./GNOVA/gnova.py ./GLC_info.sumstats.gz ./RD_info.sumstats.gz \
--N1 68390 \
--N2 64011 \
--bfile ./1000G_EUR_maf001/1000G_eur_chr@_maf001 \
--out ./GNOVA/GLC_RD_1000G_maf001.txt

python ./GNOVA/gnova.py ./GLC_info.sumstats.gz ./Myopia_info.sumstats.gz \
--N1 68390 \
--N2 64268 \
--bfile ./1000G_EUR_maf001/1000G_eur_chr@_maf001 \
--out ./GNOVA/GLC_Myopia_1000G_maf001.txt

python ./GNOVA/gnova.py ./RD_info.sumstats.gz ./Myopia_info.sumstats.gz \
--N1 64011 \
--N2 64268 \
--bfile ./1000G_EUR_maf001/1000G_eur_chr@_maf001 \
--out ./GNOVA/RD_Myopia_1000G_maf001.txt
