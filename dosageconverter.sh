# create the pheno.ped by "  awk '{print $1"_"$2,$1"_"$2,"0 0",$5,$6}' xxx.bam > pheno.ped "

ml dosageconverter
ml  slurm
ml mach2dat
password='tXD0WJSzF8kXwi';  ## provided by email.

Rsq_cutoff=0.3;

#echo "TRAIT           MARKER          ALLELES  FREQ1    RSQR   EFFECT1  OR      STDERR  WALDCHISQ PVALUE     LRCHISQ LRPVAL" >omrfIC_AA_mach2dat.Rsq30


for chr in {1..22};
do
{

jobID=$(sbatch --wrap " unzip -P '$password'  chr_${chr}.zip" | cut -f4 -d ' ' )

##sbatch --mem 30G --wrap " unzip -P '$password' chr_${chr}.zip"
#unzip -P $password chr_${chr}.zip


## convert vcf formach into mach dosage format
jobID=$(sbatch --mem 5G --dependency=afterok:$jobID  --wrap "DosageConvertor --vcfDose  chr${chr}.dose.vcf.gz --info  chr${chr}.info.gz --type mach --prefix chr${chr} --format  1 "  | cut -f4 -d ' ' )

##sbatch --mem 50G --dependency=afterok:$jobID  --wrap "mach2dat  -d pheno.dat -p pheno.ped -i  ${region}_tmp2.info -g  ${region}.dose.gz --rsqcutoff $Rsq_cutoff --frequency > AAgwas_chr${chr}_mach2dat.Rsq${Rsq_cutoff}
#sbatch --mail-type=ALL --mail-user=sunc@omrf.org --mem 20G --dependency=afterok:$jobID  --wrap "mach2dat  -d pheno.dat -p pheno_admixK5-4.ped -i  chr${chr}.mach.info -g  chr${chr}.mach.dose.gz --rsqcutoff $Rsq_cutoff --frequency > chr${chr}_mach2dat_admixK5_4.Rsq${Rsq_cutoff} "

#sbatch --mail-type=ALL --mail-user=sunc@omrf.org --mem 20G --dependency=afterok:$jobID  --wrap "mach2dat  -d pheno.dat -p pheno_admixK5-4_2ctrls.ped -i  chr${chr}.mach.info -g  chr${chr}.mach.dose.gz --rsqcutoff $Rsq_cutoff --frequency > chr${chr}_mach2dat_2ctrls_admixK5_4.Rsq${Rsq_cutoff} "

#sbatch --mail-type=ALL --mail-user=sunc@omrf.org --mem 20G --dependency=afterok:$jobID  --wrap "mach2dat  -d pheno.dat -p 386VS386_allCtrls.ped  -i  chr${chr}.mach.info -g  chr${chr}.mach.dose.gz --rsqcutoff $Rsq_cutoff --frequency > chr${chr}_mach2dat_386VS386_allCtrls_admixK5_4.Rsq${Rsq_cutoff} "


#sbatch --mail-type=ALL --mail-user=sunc@omrf.org --mem 20G --dependency=afterok:$jobID  --wrap "mach2dat  -d pheno.dat -p 212VS211_addctrls.ped -i  chr${chr}.mach.info -g  chr${chr}.mach.dose.gz --rsqcutoff $Rsq_cutoff --frequency > chr${chr}_mach2dat_212VS211_addctrls_admixK5_4.Rsq${Rsq_cutoff} "


sbatch --mail-type=ALL --mail-user=sunc@omrf.org --mem 20G --wrap "mach2dat  -d pheno.dat -p 175VS174_MEGActrls.ped  -i  chr${chr}.mach.info -g  chr${chr}.mach.dose.gz --rsqcutoff $Rsq_cutoff --frequency > chr${chr}_mach2dat_175VS174_MEGActrl_admixK5_4.Rsq${Rsq_cutoff} "



}
done


