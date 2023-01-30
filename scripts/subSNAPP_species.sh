#!/bin/bash
#SBATCH --cpus-per-task=24
#SBATCH --mem 48GB
#SBATCH -p fast
#SBATCH -o snapp.%N.%j.out
#SBATCH -e snapp.%N.%j.err

echo '########################################'
echo 'Date:' $(date --iso-8601=seconds)
echo 'User:' $USER
echo 'Host:' $HOSTNAME
echo 'Job Name:' $SLURM_JOB_NAME
echo 'Job Id:' $SLURM_JOB_ID
echo 'Directory:' $(pwd)
echo '########################################'

module load bcftools vcftools beast2
source /shared/ifbstor1/software/miniconda/etc/profile.d/conda.sh
conda activate ruby

##https://github.com/mmatschiner/tutorials/blob/master/species_tree_inference_with_snp_data/README.md
#bcftools mpileup -Ou -f /shared/projects/messor_wgs/messorBams/MaciRef/MessorNoMales_VCF/references/Maci.fasta --bam-list qsl2_cramFiles | bcftools call --threads 24 -mv -Oz -o qsl2_MaciRef_RNAseq_cram.vcf.gz 
#bcftools view -H qsl2_MaciRef_RNAseq_cram.vcf.gz | wc -l 

#bcftools filter --threads 24 -s LowQual -e '%QUAL<20 && INFO/DP<10' -Oz -o qsl2_MaciRef_RNAseq_cram_qual20dp10.vcf.gz qsl2_MaciRef_RNAseq_cram.vcf.gz

#bcftools view --threads 24 -e 'AC==0 || AC==AN || F_MISSING > 0.2' -m2 -M2 -O z -o qsl2_MaciRef_RNAseq_cram_filtered.vcf.gz qsl2_MaciRef_RNAseq_cram_qual20dp10.vcf.gz
#bcftools view -H qsl2_MaciRef_RNAseq_cram_filtered.vcf.gz | wc -l # 42819 SNPs

#vcftools --gzvcf qsl2_MaciRef_RNAseq_cram_filtered.vcf.gz  --recode --thin 100 --out qsl2_MaciRef_RNAseq_cram_filteredThin  # 35450 SNP 
#mv qsl2_MaciRef_RNAseq_cram_filteredThin.recode.vcf qsl2_MaciRef_RNAseq_cram_filteredThin.vcf

#ruby ../../snapp_prep.rb -v qsl2_MaciRef_RNAseq_cram_filteredThin.vcf -t qsl2_popID.txt -c constraintsSP.txt -m 1000 -l 100000 -s qsl2_busco_concatenatedIPGMA_sp.tree -o qsl2_sp

wait
beast -threads 24 -resume snapp.xml 
