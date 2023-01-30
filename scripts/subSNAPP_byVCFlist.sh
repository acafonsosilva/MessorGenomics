#!/bin/bash
#SBATCH --cpus-per-task=24
#SBATCH --mem 48GB
#SBATCH -p fast
#SBATCH -o snapp.%N.%j.out
#SBATCH -e snapp.%N.%j.err

module load bcftools vcftools beast2
source /shared/ifbstor1/software/miniconda/etc/profile.d/conda.sh
conda activate ruby

### test to be sure hardFiltered vcflist is comparable with cram output

#bcftools merge -l qsl2_vcfList -m all --threads 24 -Oz  -o merged.vcf.gz
#tabix -p vcf merged.vcf.gz
#bcftools view --threads 24 -R filtered.bed -e 'AC==0 || AC==AN || F_MISSING > 0.1' --types snps -m2 -M2 -Oz -o qsl2_MaciRef_RNAseq_vcf_qual20dp10.vcf.gz merged.vcf.gz

#vcftools --gzvcf qsl2_MaciRef_RNAseq_vcf_qual20dp10.vcf.gz  --recode --thin 100 --out qsl2_MaciRef_RNAseq_vcf_filteredThin 
#mv qsl2_MaciRef_RNAseq_vcf_filteredThin.recode.vcf qsl2_MaciRef_RNAseq_vcf_filteredThin.vcf

#ruby ../../snapp_prep.rb -v qsl2_MaciRef_RNAseq_vcf_filteredThin.vcf -t ../spSNAPP/qsl2_popID.txt  -c ./spSNAPP/constraintsSP.txt -m 2000 -l 500000 -s ./spSNAPP/qsl2_busco_concatenatedIPGMA_sp.tree -o qsl2_sp_vcf

wait
beast -threads 24 snapp.xml
