#!/bin/sh -l

#SBATCH -J 1000GenomesVCF
#SBATCH -o 1000g_chr"$chr"
#SBATCH --mail-user sai-gottam@0mrf.org
#SBATCH --mail-type=ALL
#SBATCH -p serial
#mem 10G
#SBATCH --cpus-per-task 4

#ml gcta
#ml plink
ml vcftools

#Merging of two files using plink merge command

#plink  --bfile GWAS1HRSprune --bmerge 1kg_chr.prune.bed 1kg_chr.prune.bim 1kg_chr.prune.fam --geno 0.95  --out GWASPCA1 --make-bed

#plink --bfile AAGWAS3ADDMerge --maf 0.05 --geno 0.95 --bmerge 1kg_chr_1.prune.bed 1kg_chr_1.prune.bim 1kg_chr_1.prune.fam --out AAGWAS5ADDMerge --make-bed

for chr  in {1..22}

do
{
#echo "$chr"

geno_file=/Volumes/hts_core/Shared/1000g/phase3_20130502/ALL.chr"$chr".phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz

#echo $geno_file

vcftools  --gzvcf $geno_file --positions positions.txt --plink-tped  --out 1000g_chr"$chr"

}

done


