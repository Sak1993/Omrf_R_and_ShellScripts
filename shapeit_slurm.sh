#!/bin/bash -l

#SBATCH -J MAC16MERGE
#SBATCH -o MAC16MERGE
#SBATCH --mail-user Sai-gottam@0mrf.org
#SBATCH --mail-type=ALL
#SBATCH -p serial
#SBATCH --mem 10G
#SBATCH --cpus-per-task 10




ml plink
ml shapeit
ml slurm

#file="GSA_5plates_Mexican_IBD25_uniq";
##plink --bfile  $file  --make-bed --hwe 0.0001 --maf 0.001  --geno 0.05 --allow-no-sex --out ${file}_SNPQC 

for chr in {1..22};

do 
{
map=/Volumes/hts_core/Shared/impute2_reference_panels/1000GP_Phase3_20141028/genetic_map_chr${chr}_combined_b37.txt;
hap=/Volumes/hts_core/Shared/impute2_reference_panels/1000GP_Phase3_20141028/1000GP_Phase3_chr${chr}.hap.gz; 
legend=/Volumes/hts_core/Shared/impute2_reference_panels/1000GP_Phase3_20141028/1000GP_Phase3_chr${chr}.legend.gz;
sample=/Volumes/hts_core/Shared/1000g/phase3_haplotypes_20141020/1000GP_Phase3.sample;


#if [[ ! -s  ${file}_chr${chr}.haps ]] ; then
#         echo ${chr}.haps  is empty
#       echo "$chr $pos1 $pos2 $chr "  >>empty.txt
#	continue
#fi


#plink --bfile GSA12plate_Mex --chr $chr  --make-bed --hwe 0.0001 --maf 0.001 --geno 0.1 --mind 0.2 --out GSA12plate_Mex_chr${chr}

#echo $Plink
#echo $Plink


shapeit -B GSA12plate_Mex_chr${chr} -M $map --input-ref $hap $legend $sample -T 8  --exclude-snp GSA12plate_Mex_chr${chr}.snp.strand.exclude -O GSA12plate_Mex_chr${chr}

}


done;

