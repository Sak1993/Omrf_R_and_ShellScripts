#!/bin/bash -l


#SBATCH -J MAC1111PLINK
#SBATCH -o MAC1111PLINK
#SBATCH --mail-user Sai-gottam@0mrf.org
#SBATCH --mail-type=ALL
#SBATCH -p serial
#SBATCH --mem 30G
#SBATCH --cpus-per-task 10


ml gtool
ml qctool
ml snptest


cat  hg19_autoChrsForImpute2.txt | while read chr pos1 pos2 p; do 
{

echo "$chr $pos1 $pos2 $p"; 

if [ -f ${p}.impute2.gz ];
then
{

awk '{if($6==0){ $6=0 }  if($7==1){ $7=0; print} else if($7==2){ $7=1;print} else print}'  GSA12plate_Mex_chr${chr}.sample  >${p}_recode.sample 
qctool -g ${p}.impute2.gz -s ${p}_recode.sample   -og  qc${p}.impute2  -maf 0.001 0.5  -info 0.3 1 -omit-chromosome 
awk -v c=$chr '{$1=c;print} ' qc${p}.impute2 >qc${p}.impute2.tmp
mv qc${p}.impute2.tmp qc${p}.impute2
gzip qc${p}.impute2
snptest -data qc${p}.impute2.gz  ${p}_recode.sample -o ${p}.snptest -frequentist 1  -method em -pheno plink_pheno -hwe 
sleep 0.1 ;
}

else
{
echo " file_${p}  missing" ;
echo "$chr $pos1 $pos2 $p">>snptest.missing;
}
fi

sleep 0.05;
}

done







