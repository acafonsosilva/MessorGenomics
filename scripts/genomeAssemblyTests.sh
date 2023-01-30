#### Analyses tested for genome assembly
#### this script is very messy, it's mostly a compilation of bits of code used until finish the assembly. Some bits have the code repeated for the 3 species and some other just for one example.

tar -xvf P201SC17122780-02-F001_1.tar
tar -xvf P201SC17122780-02-F001_2.tar


## convert subreads.bam to fasta and fastqz files
module load bam2fastx ## for some reason doesn’t call it without the full path  ###https://github.com/PacificBiosciences/bam2fastx
cd /shared/projects/messor_wgs/P201SC17122780-02-F001_2/raw_data/2_B03/
/shared/ifbstor1/software/miniconda/envs/bam2fastx-1.3.0/bin/bam2fasta m54217_201229_071552.subreads.bam -o Mtar_m54217_201229_071552
/shared/ifbstor1/software/miniconda/envs/bam2fastx-1.3.0/bin/bam2fastq m54217_201229_071552.subreads.bam -o Mtar_m54217_201229_071552
mv  *.gz  /shared/projects/messor_wgs/pacbio_fastx

cd /shared/projects/messor_wgs/P201SC17122780-02-F001_1/raw_data/1_B02/
/shared/ifbstor1/software/miniconda/envs/bam2fastx-1.3.0/bin/bam2fasta m54217_201228_110217.subreads.bam -o Maeg_m54217_201228_110217 
/shared/ifbstor1/software/miniconda/envs/bam2fastx-1.3.0/bin/bam2fastq m54217_201228_110217.subreads.bam -o Maeg_m54217_201228_110217
mv  *.gz  /shared/projects/messor_wgs/pacbio_fastx

cd /shared/projects/messor_wgs/P201SC17122780-02-F001_1/raw_data/1_D02/
/shared/ifbstor1/software/miniconda/envs/bam2fastx-1.3.0/bin/bam2fasta  m54217_201226_103221.subreads.bam -o Mibe_m54217_201226_103221
/shared/ifbstor1/software/miniconda/envs/bam2fastx-1.3.0/bin/bam2fastq m54217_201226_103221.subreads.bam -o Mibe_m54217_201226_103221
mv  *.gz  /shared/projects/messor_wgs/pacbio_fastx


### if need we can also convert fastq to fasta with seqkit
#module load seqkit
#seqkit fq2fa  Mtar_m54217_201229_071552.fastq.gz > Mtar_m54217_201229_071552.fasta

### check quality   # https://github.com/yfukasawa/LongQC
module load longqc/1.1.1 
srun --mem 80GB --cpus-per-task=12 longQC.py sampleqc -x pb-sequel Maeg_m54217_201228_110217.fastq.gz -o Maeg -s Maeg -p 12 --index 200M
srun --mem 80GB --cpus-per-task=12 longQC.py sampleqc -x pb-sequel Mibe_m54217_201226_103221.fastq.gz -o Mibe -s Mibe -p 12 --index 200M
srun --mem 80GB --cpus-per-task=12 longQC.py sampleqc -x pb-sequel Mtar_m54217_201229_071552.fastq.gz -o Mtar -s Mtar -p 12 --index 200M


##############################
########## Assembly ##########
##############################
##Wtdbg2
##Canu
##Masurca
##raven

##Wtdbg2
module load minimap2
module load samtools
module load bwa

/shared/home/asilva/wtdbg2/wtdbg2 -x sq -g 310m -p 0 -k 15 -AS 2 -s 0.05 -L 5000 -t 16 -i /shared/projects/messor_wgs/pacbio_fastx/Mibe_m54217_201226_103221.fastq.gz -fo Mibe5kb
/shared/home/asilva/wtdbg2/wtpoa-cns -t 16 -i Mibe5kb.ctg.lay.gz -fo Mibe5kb.ctg.fa

## initially assemblies were run with 10kb minimum reads given something some ant genome paper mentioning they got better results but eventually I compared using 5Kb which is the advised by the developers when using long reads. 
## used -g parameter - around 310 Mb given outcome from gaga project but I don't think the size changed much the outcome even thought these assemblies are smaller


##Canu

module load java-jdk/8.0.112 
srun --mem 200GB --cpus-per-task=32 --partition long /shared/home/asilva/canu-2.1.1/bin/canu -p Maeg genomeSize=254m -pacbio /shared/projects/messor_wgs/pacbio_fastx/Maeg_m54217_201228_110217.fastq.gz

srun --mem 200GB --cpus-per-task=32 --partition long /shared/home/asilva/canu-2.1.1/bin/canu -p Mibe genomeSize=256m -pacbio /shared/projects/messor_wgs/pacbio_fastx/Mibe_m54217_201226_103221.fastq.gz

srun --mem 200GB --cpus-per-task=32 --partition long /shared/home/asilva/canu-2.1.1/bin/canu  -p Mbar genomeSize=263m -pacbio /shared/projects/messor_wgs/pacbio_fastx/Mtar_m54217_201229_071552.fastq.gz

##Masurca

# to get library insert average length:
module load java-jdk/8.0.112 
/shared/home/asilva/bbmap/bbmerge-auto.sh in1=Maeg_A5_FDSW202587019-1r_HW2H2DSXY_L1_1.fq.gz in2=Maeg_A5_FDSW202587019-1r_HW2H2DSXY_L1_2.fq.gz ihist=Maeg_ihist.txt prefilter=2 rem extend2=100 k=62

##genome size to be used for Masurca:
#Maeg 254 -  287 insert size average
#Mibe 256 -  286 insert size average
#Mtar 263 -  313 insert size average

module load masurca
masurca config.txt  ## got deleted but the resulting assemble.sh file was kept
./assemble.sh  ## or sbatch subMasurca.sh

##raven
module load raven-assembler/1.4.0

srun --mem 50GB --cpus-per-task=32 raven -t 32 /shared/projects/messor_wgs/pacbio_fastx/Maeg_m54217_201228_110217.fastq.gz > Maeg.fasta
srun --mem 50GB --cpus-per-task=32 raven -t 32 /shared/projects/messor_wgs/pacbio_fastx/Mibe_m54217_201226_103221.fastq.gz > Mibe.fasta
srun --mem 50GB --cpus-per-task=32 raven -t 32 /shared/projects/messor_wgs/pacbio_fastx/Mtar_m54217_201229_071552.fastq.gz > Mtar.fasta

###############################
########## Polishing ##########
###############################
##Polca
##Ratatosk
##Wtdbg2Suggested

##NextPolish


##Polca
cd /shared/projects/messor_wgs/assemblies/wtdbg2/polishing/polca/Mibe/
module load masurca
module load  bwa
module load samtools 
module load freebayes
srun --mem 34GB --cpus-per-task 16 polca.sh -a  ../../../Mibe10kb/Mibe10kb.ctg.fa -r '/shared/projects/messor_wgs/illumina_fastq/clean_reads/Mibe_A3_cleaned_READ1.fastq.gz  /shared/projects/messor_wgs/illumina_fastq/clean_reads/Mibe_A3_cleaned_READ2.fastq.gz' -t 16 -m 2G

##Ratatosk

module load seqtk/ ### need in .fasta format

cd /shared/projects/messor_wgs/assemblies/wtdbg2/polishing/Ratatosk/Maeg/
srun --mem 20GB  --cpus-per-task=32 /shared/home/asilva/Ratatosk/build/src/Ratatosk -v -c 32 -i 287 -s /shared/projects/messor_wgs/illumina_fastq/clean_reads/Maeg_A5_cleaned_READ1.fastq.gz /shared/projects/messor_wgs/illumina_fastq/clean_reads/Maeg_A5_cleaned_READ2.fastq.gz -l /shared/projects/messor_wgs/assemblies/wtdbg2/Maeg10kb/Maeg10kb.ctg.fa -o Maeg

seqtk seq -a Maeg.fastq > Maeg.fasta


##Wtdbg2Suggested
sample=Mibe
srun --mem 100GB --cpus-per-task 24 minimap2 -t24 -ax map-pb -r2k  ../../../$sample'10kb/'$sample'10kb.ctg.fa' /shared/projects/messor_wgs/pacbio_fastx/Mibe_m54217_201226_103221.fastq.gz | samtools sort -@4 > $sample'10kb.bam'
srun --mem 100GB --cpus-per-task 24 samtools view -F0x900 $sample'10kb.bam' | $wtdbg2_path/wtpoa-cns -t 24 -d ../../../$sample'10kb/'$sample'10kb.ctg.fa' -i - -fo $sample'10kb.cns.fa'
srun --mem 100GB bwa index $sample'10kb.cns.fa'
srun --mem 100GB --cpus-per-task 24 bwa mem -t 24 $sample'10kb.cns.fa'  ../../../../../illumina_fastq/clean_reads/$sample'_A3_cleaned_READ1.fastq.gz' ../../../../../illumina_fastq/clean_reads/$sample'_A3_cleaned_READ2.fastq.gz' | samtools sort -O SAM | $wtdbg2_path/wtpoa-cns -t 24 -x sam-sr -d $sample'10kb.cns.fa' -i - -fo $sample'10kb.srp.fa'

##nextpolish ran in seed ### 

module load nextPolish

cd /home/ac.afonso/working_dir/polishing2x/nextpolish

ls /home/ac.afonso/working_dir/clean_reads/Mibe_A3_cleaned_READ1.fastq /home/ac.afonso/working_dir/clean_reads/Mibe_A3_cleaned_READ2.fastq > Mibe_A3_cleaned.fofn

ls /home/ac.afonso/working_dir/pacbio_fastx/Mibe_m54217_201226_103221.fastq.gz > Mibe.fastq.fofn


##make config file Mibe_wtdbg2_wtdbg2Suggested_run.cfg

#[General]
#job_type = local
#job_prefix = nextPolish
#task = best
#rewrite = yes
#rerun = 3
#parallel_jobs = 6
#multithread_jobs = 5
#genome = /home/ac.afonso/working_dir/assemblies/assembled_polished1x_fastas/Mibe_wtdbg2_wtdbg2Suggested.fasta
#genome_size = auto
#workdir = Mibe_wtdbg2_wtdbg2Suggested_nextpolish
#polish_options = -p {multithread_jobs}
#
#[sgs_option]
#sgs_fofn = Mibe_A3_cleaned.fofn
#sgs_options = -max_depth 100 -bwa
#
#[lgs_option]
#lgs_fofn = Mibe.fastq.fofn
#lgs_options = -min_read_len 1k -max_depth 100
#lgs_minimap2_options = -x map-pb

nextPolish Mibe_wtdbg2_wtdbg2Suggested_run.cfg

cat /home/ac.afonso/working_dir/polishing2x/nextpolish/Mibe_wtdbg2_wtdbg2Suggested_nextpolish/*.kmer_count/*.polish.ref.sh.work/polish_genome*/genome.nextpolish.*.fasta > /home/ac.afonso/working_dir/polishing2x/nextpolish/Mibe_wtdbg2_wtdbg2Suggested_nextpolish/Mibe_wtdbg2_wtdbg2Suggested_nextpolish.fasta

cp -r /home/ac.afonso/working_dir/polishing2x/nextpolish/Mibe_wtdbg2_wtdbg2Suggested_nextpolish/Mibe_wtdbg2_wtdbg2Suggested_nextpolish.fasta /home/ac.afonso/working_dir/assemblies/assembled_polished2x_fastas/


################################
###### merging assemblies ######
################################

### ran in Seed to test if merging Canu and wtdbg2 assemblies would give better assemblies but decided to just finish the draft genome with wtdbg2 
module load quickmerge

merge_wrapper.py ../../assemblies/assembled_polished2x_fastas/Mibe_wtdbg2_polca_nextpolish.fasta ../../assemblies/assembled_polished2x_fastas/Mibe_canu_polca_nextpolish.fasta  -p merged_Mibe_wtdbg2_canu_polcaNextpolishx2


#################################
########## Scaffolding ##########
#################################

##LRScaf - without reference
##RagTagScaf - with reference

################################
###### LRScaf scaffolding ######
################################

path=/shared/projects/messor_wgs
analysis=Mtar_wtdbg2_polca_nextpolish
draft=$path/assemblies/seed_draftAssemblies/$analysis.fasta
subreads=$path/pacbio_fastx/Mtar_m54217_201229_071552.fastq.gz
subreads_fasta=$path/pacbio_fastx/Mtar_m54217_201229_071552.fasta
illumina1=$path/illumina_fastq/clean_reads/Mtar_A7_cleaned_READ1.fastq.gz
illumina2=$path/illumina_fastq/clean_reads/Mtar_A7_cleaned_READ2.fastq.gz

source /shared/ifbstor1/software/miniconda/etc/profile.d/conda.sh
module load minimap2
module load java-jdk/8.0.112
module load masurca
module load bwa

### map draft to subreads ###
mkdir $path/scaffolding/LRScaf/$analysis
cd $path/scaffolding/LRScaf/$analysis


echo 'Start mapping draft assembly to subreads'
minimap2 -t 16 $draft $subreads > $analysis'_draft.mm'
wait 

echo 'Start running LRScaf'
java -Xms100g -Xmx100g -jar /shared/home/asilva/LRScaf-1.1.10.jar -c $draft -a $analysis'_draft.mm' -t mm -o $path/scaffolding/LRScaf/$analysis -i 0.1 -r 0.2 -misl 1 -micl 500 -mioll 400 -miolr 0.8 -mxel 1000 -mxer 0.1 -mxohl 1000 -mxohr 0.1 -mmcm 10 -p 16 ### for whatever reason -iqrt cannot be set as a command line
wait

mv scaffolds.fasta $analysis'_LRScaf.fasta'

## polca
echo 'Start running polca after LRScaf'
polca.sh -a  $analysis'_LRScaf.fasta' -r $illumina1\ $illumina2 -t 16 -m 2G
wait

mv $analysis'_LRScaf.fasta.PolcaCorrected.fa' $analysis'_LRScaf_polca.fasta'
cp -r $analysis'_LRScaf_polca.fasta' $path/assemblies/assembled_scaffolds_fastas/
rm -rf $analysis'_LRScaf.fasta.'*
rm *draft.mm


####################################
###### RagTagScaf scaffolding ######
####################################

source /shared/ifbstor1/software/miniconda/etc/profile.d/conda.sh
module load minimap2
module load masurca
module load bwa
module load ragtag

mkdir $path/scaffolding/RagTag/$analysis
cd $path/scaffolding/RagTag/$analysis

srun --cpus-per-task=16 --mem 50GB ragtag.py scaffold $path/assemblies/reference_assembly/GAGA-0413_Messor_capitatus.fasta $path/assemblies/seed_draftAssemblies/$analysis'.fasta' -o $analysis'_RagTagScaf' -u
wait

cp -r $analysis'_RagTagScaf/ragtag.scaffolds.fasta' $analysis'_RagTagScaf.fasta'

## polca
echo 'Start running polca after RagTagScaf'
polca.sh -a  $analysis'_RagTagScaf.fasta' -r $illumina1\ $illumina2 -t 16 -m 2G
wait

mv $analysis'_RagTagScaf.fasta.PolcaCorrected.fa' $analysis'_RagTagScaf_polca.fasta'
cp -r $analysis'_RagTagScaf_polca.fasta' $path/assemblies/assembled_scaffolds_fastas/
rm -rf *_RagTagScaf.fasta.*


#################################
########## Gap closing ##########
#################################

##TGS-GapCloser
##LR_Gapcloser

### Maeg
echo 'Start running Maeg'
analysis=Maeg_wtdbg2_polca_nextpolish
draft=$path/assemblies/seed_draftAssemblies/$analysis.fasta
subreads_fasta=$path/pacbio_fastx/Maeg_m54217_201228_110217.fasta
illumina1=$path/illumina_fastq/clean_reads/Maeg_A5_cleaned_READ1.fastq.gz
illumina2=$path/illumina_fastq/clean_reads/Maeg_A5_cleaned_READ2.fastq.gz

### Mibe

echo 'Start running Mibe'
analysis=Mibe_wtdbg2_polca_nextpolish
draft=$path/assemblies/seed_draftAssemblies/$analysis.fasta
subreads_fasta=$path/pacbio_fastx/Mibe_m54217_201226_103221.fasta
illumina1=$path/illumina_fastq/clean_reads/Mibe_A3_cleaned_READ1.fastq.gz
illumina2=$path/illumina_fastq/clean_reads/Mibe_A3_cleaned_READ2.fastq.gz

### Mtar
echo 'Start running Mtar'
analysis=Mtar_wtdbg2_polca_nextpolish
draft=$path/assemblies/seed_draftAssemblies/$analysis.fasta
subreads_fasta=$path/pacbio_fastx/Mtar_m54217_201229_071552.fasta
illumina1=$path/illumina_fastq/clean_reads/Mtar_A7_cleaned_READ1.fastq.gz
illumina2=$path/illumina_fastq/clean_reads/Mtar_A7_cleaned_READ2.fastq.gz

##########################
###### LR_Gapcloser ######
##########################

##### LR_Gapcloser with LRScaf scaffolding ######

mkdir $path/gapClosing/LR_Gapcloser/$analysis'_LRScaf'
cd $path/gapClosing/LR_Gapcloser/$analysis'_LRScaf'

echo 'Start running LR_Gapcloser'
srun --cpus-per-task=16 --mem 50GB /shared/home/asilva/LR_Gapcloser/src/LR_Gapcloser.sh -i $path/assemblies/assembled_scaffolds_fastas/$analysis'_LRScaf_polca.fasta' -l $subreads -s p -o $analysis'_LRScaf_polca_LR' -t 16
wait

cp -r $analysis'_LRScaf_polca_LR'/iteration-3/gapclosed.fasta $analysis'_LRScaf_polca_LR.fasta'


## polca
echo 'Start running polca after LR_Gapcloser'
srun --cpus-per-task=16 --mem 50GB polca.sh -a $analysis'_LRScaf_polca_LR.fasta' -r $illumina1\ $illumina2 -t 16 -m 2G
wait

mv $analysis'_LRScaf_polca_LR.fasta.PolcaCorrected.fa' $analysis'_LRScafpolca_LRpolca.fasta'
cp -r $analysis'_LRScafpolca_LRpolca.fasta' $path/assemblies/assembled_scaffolds_gapsClosed_fastas/
rm -rf *.fasta.*


## quast
module load quast
echo 'Start running quast of LR_Gapcloser + polca'
quast $analysis'_LRScafpolca_LRpolca.fasta' -o quast --large


## busco
module load busco
echo 'Start running busco of LR_Gapcloser + polca'
srun --cpus-per-task=16 --mem 50GB busco -i $analysis'_LRScafpolca_LRpolca.fasta' -o busco.$analysis'_LRScafpolca_LRpolca' -c 16 -m geno -l $path/buscoScoring/busco_downloads/lineages/hymenoptera_odb10
module unload busco


###### LR_Gapcloser with RagTagScaf scaffolding ######

mkdir $path/gapClosing/LR_Gapcloser/$analysis'_RagTagScaf'
cd $path/gapClosing/LR_Gapcloser/$analysis'_RagTagScaf'

echo 'Start running LR_Gapcloser'
/shared/home/asilva/LR_Gapcloser/src/LR_Gapcloser.sh -i $path/assemblies/assembled_scaffolds_fastas/$analysis'_RagTagScaf_polca.fasta' -l $subreads -s p -o $analysis'_RagTagScaf_polca_LR' -t 16
wait

cp -r $analysis'_RagTagScaf_polca_LR'/iteration-3/gapclosed.fasta $analysis'_RagTagScaf_polca_LR.fasta'

## polca
echo 'Start running polca after LR_Gapcloser'
polca.sh -a $analysis'_RagTagScaf_polca_LR.fasta' -r $illumina1\ $illumina2 -t 16 -m 2G
wait

mv $analysis'_RagTagScaf_polca_LR.fasta.PolcaCorrected.fa' $analysis'_RagTagScafpolca_LRpolca.fasta'
cp -r $analysis'_RagTagScafpolca_LRpolca.fasta' $path/assemblies/assembled_scaffolds_gapsClosed_fastas/
rm -rf *.fasta.*

## quast
module load quast
echo 'Start running quast of LR_Gapcloser + polca'
quast $analysis'_RagTagScafpolca_LRpolca.fasta' -o quast --large
wait

## busco
module load busco
echo 'Start running busco of LR_Gapcloser + polca'
busco -i $analysis'_RagTagScafpolca_LRpolca.fasta' -o busco.$analysis'_RagTagScafpolca_LRpolca' -c 16 -m geno -l $path/buscoScoring/busco_downloads/lineages/hymenoptera_odb10
module unload busco


#####################
### TGS-GapCloser ###
#####################
conda activate tgsgapcloser

##### TGS-GapCloser with LRScaf scaffolding ######

#module load seqkit
#seqkit fq2fa $subreads > $subreads_fasta

mkdir $path/gapClosing/tgsgapcloser/$analysis'_LRScaf'
cd $path/gapClosing/tgsgapcloser/$analysis'_LRScaf'

echo 'Start running TGS-GapCloser'
srun --mem 100GB --cpus-per-task 16 tgsgapcloser --scaff $path/assemblies/assembled_scaffolds_fastas/$analysis'_LRScaf_polca.fasta' --reads $subreads_fasta --output $analysis'_LRScaf_polca_TGS' --racon /shared/home/asilva/racon-v1.4.21/build/bin/racon --tgstype pb --thread 16 --r_round 3
rm -rf *paf

## polca
echo 'Start running polca after TGS-GapCloser'
polca.sh -a $analysis'_LRScaf_polca_TGS.scaff_seqs' -r $illumina1\ $illumina2 -t 16 -m 2G
wait

mv $analysis'_LRScaf_polca_TGS.scaff_seqs.PolcaCorrected.fa' $analysis'_LRScafpolca_TGSpolca.fasta'
cp -r $analysis'_LRScafpolca_TGSpolca.fasta' $path/assemblies/assembled_scaffolds_gapsClosed_fastas/
rm -rf *_LRScaf_polca_TGS.scaff_seqs.*   ### still to check if works well

## quast
module load quast
echo 'Start running quast of TGS-GapCloser + polca'
quast $analysis'_LRScafpolca_TGSpolca.fasta' -o quast --large


## busco
module load busco
echo 'Start running busco of TGS-GapCloser + polca'
srun --cpus-per-task=16 --mem 50GB busco -i $analysis'_LRScafpolca_TGSpolca.fasta' -o busco.$analysis'_LRScafpolca_TGSpolca' -c 16 -m geno -l $path/buscoScoring/busco_downloads/lineages/hymenoptera_odb10


##### TGS-GapCloser with RagTagScaf scaffolding ######

path=/shared/projects/messor_wgs

source /shared/ifbstor1/software/miniconda/etc/profile.d/conda.sh
module load minimap2
module load masurca
module load bwa
module load ragtag

mkdir $path/gapClosing/tgsgapcloser/$analysis'_RagTagScaf'
cd $path/gapClosing/tgsgapcloser/$analysis'_RagTagScaf'

echo 'Start running TGS-GapCloser'
srun --mem 100GB --cpus-per-task 16 tgsgapcloser --scaff $path/assemblies/assembled_scaffolds_fastas/$analysis'_RagTagScaf_polca.fasta' --reads $subreads_fasta --output $analysis'_RagTagScaf_polca_TGS' --racon /shared/home/asilva/racon-v1.4.21/build/bin/racon --tgstype pb --thread 16 --r_round 3
wait
rm -rf *paf
rm -rf *.ont*.fasta

## polca
echo 'Start running polca after TGS-GapCloser'
polca.sh -a $analysis'_RagTagScaf_polca_TGS.scaff_seqs' -r $illumina1\ $illumina2 -t 16 -m 2G
wait

mv $analysis'_RagTagScaf_polca_TGS.scaff_seqs.PolcaCorrected.fa' $analysis'_RagTagScafpolca_TGSpolca.fasta'
cp -r $analysis'_RagTagScafpolca_TGSpolca.fasta' $path/assemblies/assembled_scaffolds_gapsClosed_fastas/
rm -rf *.scaff_seqs.*


### The difference between LR_Gapcloser and TGS-GapCloser is very marginal with LR_Gapcloser much easier and faster to use. 

#################################
##### Assemblers comparison #####
#################################

module load quast
module load bwa
module load bedtools

path=/shared/projects/messor_wgs
drafts=$path/assemblies/assembled_scaffolds_gapsClosed_fastas
ref=$path/assemblies/reference_assembly/GAGA-0413_Messor_capitatus.fasta
subreads=$path/pacbio_fastx/Mtar_m54217_201229_071552.fastq.gz
illumina1=$path/illumina_fastq/clean_reads/Mtar_A7_cleaned_READ1.fastq.gz
illumina2=$path/illumina_fastq/clean_reads/Mtar_A7_cleaned_READ2.fastq.gz

srun --cpus-per-task=16 --mem 50GB quast -o Mtar --large -r $ref --pacbio $subreads --pe1 $illumina1 --pe2 $illumina2 -e -t 16 --fragmented --no-sv $drafts/Mtar_wtdbg2_polca_nextpolish_LRScafpolca_LRpolca.fasta $drafts/Mtar_canu_polca_nextpolish_LRScafpolca_LRpolca.fasta $drafts/merged_Mtar_wtdbg2_canu_polcaNextpolishx2_mPolcaNextpolish_LRScafpolca_LRpolca.fasta $drafts/Mtar_wtdbg2_polca_nextpolish_RagTagScafpolca_LRpolca.fasta $drafts/Mtar_wtdbg2_polca_nextpolish_LRScafpolca_TGSpolca.fasta
wait

path=/shared/projects/messor_wgs
drafts=$path/assemblies/assembled_scaffolds_gapsClosed_fastas
ref=$path/assemblies/reference_assembly/GAGA-0413_Messor_capitatus.fasta
subreads=$path/pacbio_fastx/Mibe_m54217_201226_103221.fastq.gz
illumina1=$path/illumina_fastq/clean_reads/Mibe_A3_cleaned_READ1.fastq.gz
illumina2=$path/illumina_fastq/clean_reads/Mibe_A3_cleaned_READ2.fastq.gz

srun --cpus-per-task=16 --mem 50GB quast -o Mibe --large -r $ref --pacbio $subreads --pe1 $illumina1 --pe2 $illumina2 -e -f -b -t 16 --fragmented --upper-bound-assembly --no-sv $drafts/Mibe_wtdbg2_polca_nextpolish_LRScafpolca_LRpolca.fasta $drafts/Mibe_canu_polca_nextpolish_LRScafpolca_LRpolca.fasta $drafts/merged_Mibe_wtdbg2_canu_polcaNextpolishx2_mPolcaNextpolish_LRScafpolca_LRpolca.fasta $drafts/Mibe_wtdbg2_polca_nextpolish_RagTagScafpolca_LRpolca.fasta $drafts/Mibe_wtdbg2_polca_nextpolish_LRScafpolca_TGSpolca.fasta
wait

path=/shared/projects/messor_wgs
drafts=$path/assemblies/assembled_scaffolds_gapsClosed_fastas
ref=$path/assemblies/reference_assembly/GAGA-0413_Messor_capitatus.fasta
subreads=$path/pacbio_fastx/Maeg_m54217_201228_110217.fastq.gz
illumina1=$path/illumina_fastq/clean_reads/Maeg_A5_cleaned_READ1.fastq.gz
illumina2=$path/illumina_fastq/clean_reads/Maeg_A5_cleaned_READ2.fastq.gz

srun --cpus-per-task=16 --mem 50GB quast -o Maeg --large -r $ref --pacbio $subreads --pe1 $illumina1 --pe2 $illumina2 -e -f -b -t 16 --fragmented --upper-bound-assembly --no-sv $drafts/Maeg_wtdbg2_polca_nextpolish_LRScafpolca_LRpolca.fasta $drafts/Maeg_canu_polca_nextpolish_LRScafpolca_LRpolca.fasta $drafts/merged_Maeg_wtdbg2_canu_polcaNextpolishx2_mPolcaNextpolish_LRScafpolca_LRpolca.fasta $drafts/Maeg_wtdbg2_polca_nextpolish_RagTagScafpolca_LRpolca.fasta $drafts/Maeg_wtdbg2_polca_nextpolish_LRScafpolca_TGSpolca.fasta
