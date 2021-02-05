
###impute the phrased GWAS haptypes from shapeit.
#input: the imputed regions from hg19_autoChrsForImpute2.txt and the bfile prefile name in the shapeit_slurm.sh
##output the imputed dosage data : chr*_ST_END.impute2.gz



ml slurm

file="MegaEgyptian_addGlobalCtrls_IBD25";   ## the bfile prefile name in the shapeit_slurm.sh 
rm *.impute2 *.o* *.e* empty.txt imputed_regions.txt
cat hg19_autoChrsForImpute2.txt | while read chr pos1 pos2 id;  ## read the each region for the imputation.

do 
{

map=/Volumes/hts_core/Shared/impute2_reference_panels/1000GP_Phase3_20141028/genetic_map_chr${chr}_combined_b37.txt;

hap=/Volumes/hts_core/Shared/impute2_reference_panels/1000GP_Phase3_20141028/1000GP_Phase3_chr${chr}.hap.gz; 
legend=/Volumes/hts_core/Shared/impute2_reference_panels/1000GP_Phase3_20141028/1000GP_Phase3_chr${chr}.legend.gz;

#if [[ ! -s  ${id}.haps ]] ; then
#         echo ${id}.haps  is empty
#       echo "$chr $pos1 $pos2 $id "  >>empty.txt
#	continue
#fi


echo "$chr $pos1 $pos2 $id "  >>imputed_regions.txt



##qsub -V  -N ${id}  -l h_vmem=12G,virtual_free=8000M -V -cwd  -b y  impute -use_prephased_g -known_haps_g  ${file}_chr${chr}.haps -m $map -h $hap -l  $legend -allow_large_regions  -int $pos1 $pos2    -filt_rules_l "'AFR<0.001'" -align_by_maf_g_ref   -o_gz -o ${id}.impute2

sbatch --mem=10G -p serial --wrap "impute -use_prephased_g -known_haps_g  ${file}_chr${chr}.haps -m $map -h $hap -l  $legend -allow_large_regions  -int $pos1 $pos2   -filt_rules_l 'ALL<0.002'  -align_by_maf_g_ref   -o_gz -o ${id}.impute2 "





sleep 0.1 

}


done;

