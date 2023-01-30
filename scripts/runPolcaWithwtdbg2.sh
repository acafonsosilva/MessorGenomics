#!/bin/bash
#SBATCH -p long
#SBATCH --mem 34GB
#SBATCH --cpus-per-task 16
#SBATCH -o polca.%N.%j.out
#SBATCH -e polca.%N.%j.err  

#### polishing with polca from wtdbg2

module load masurca
module load bwa
module load samtools 
module load freebayes

cd /shared/projects/messor_wgs/assemblies/wtdbg2/polishing/polca/Mibe/
srun --mem 34GB --cpus-per-task 16 polca.sh -a  ../../../Mibe10kb/Mibe10kb.ctg.fa -r '/shared/projects/messor_wgs/illumina_fastq/clean_reads/Mibe_A3_cleaned_READ1.fastq.gz  /shared/projects/messor_wgs/illumina_fastq/clean_reads/Mibe_A3_cleaned_READ2.fastq.gz' -t 16 -m 2G

cd /shared/projects/messor_wgs/assemblies/wtdbg2/polishing/polca/Mtar/
srun --mem 34GB --cpus-per-task 16 polca.sh -a  ../../../Mtar10kb/Mtar10kb.ctg.fa -r '/shared/projects/messor_wgs/illumina_fastq/clean_reads/Mtar_A7_cleaned_READ1.fastq.gz  /shared/projects/messor_wgs/illumina_fastq/clean_reads/Mtar_A7_cleaned_READ2.fastq.gz' -t 16 -m 2G

cd /shared/projects/messor_wgs/assemblies/wtdbg2/polishing/polca/Maeg
srun --mem 34GB --cpus-per-task 16 polca.sh -a  ../../../Maeg10kb/Maeg10kb.ctg.fa -r '/shared/projects/messor_wgs/illumina_fastq/clean_reads/Maeg_A5_cleaned_READ1.fastq.gz  /shared/projects/messor_wgs/illumina_fastq/clean_reads/Maeg_A5_cleaned_READ2.fastq.gz' -t 16 -m 2G
