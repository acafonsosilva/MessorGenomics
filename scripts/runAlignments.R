##xargs <- data.frame(input="/shared/projects/messor_phylo/RAxML/allNoMaleHybrid.txt", seqType="phasedH0",Nfil=0.5,bpMin=350,gapsPp=0.25)
##/shared/projects/messor_phylo/RAxML/allNoHybridPE_phasedH0_concatenatedRAxML

suppressWarnings(suppressMessages(library(seqinr)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(doFuture)))
suppressWarnings(suppressMessages(library(ape)))
suppressWarnings(suppressMessages(library(ips)))
library(argparse)
library(parallel)
registerDoFuture()
plan(multicore)

is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}


parser <- ArgumentParser()
parser$add_argument('--input', '-i')
parser$add_argument('--seqType', '-s')  ### busco | ambiguous | mitogenome | phasedH0 | phasedH01
parser$add_argument('--Nfil', '-N')
parser$add_argument('--bpMin', '-bp')
parser$add_argument('--gapsPp', '-g')
xargs <- parser$parse_args()

thr <- detectCores()
phyluce <- '/shared/home/asilva/.conda/envs/phyluce-1.7.0/bin/'
samples <- read.delim(xargs$input, header = F)
Nsamples=nrow(samples)

if(xargs$seqType == 'phasedH0' | xargs$seqType == 'phasedH01'){
  pha <- '/shared/projects/messor_dils/messor_data/messor_buscoMapped/snakePhasing/phased_fastas/'
  
  if(xargs$seqType == 'phasedH01'){
    fastas <- c(paste0(pha, samples$V1, '.0.fasta'),
                paste0(pha, samples$V1, '.1.fasta'))
  }

  if(xargs$seqType == 'phasedH0'){
    fastas <- paste0(pha, samples$V1, '.0.fasta')
  }

  bigFasta <- list()
  for (f in fastas){
    dna <- read.fasta(f, seqtype = c("DNA"),
                      seqonly = FALSE, whole.header = T,
                      forceDNAtolower = FALSE, 
                      set.attributes = TRUE)
    
    sample <- word(tools::file_path_sans_ext(f), sep ="/",-1)
    
  ids <- tibble(id = names(dna)) %>% 
    mutate(sample = sample, 
           locus = gsub("_pilon", "", word(id, sep = '@',2)), 
           newID = paste0(locus,"_", sample," |", locus))
  
  names(dna) <- ids %>% pull(newID)
  bigFasta <- c(bigFasta, dna)
  }
  write.fasta(bigFasta, names = names(bigFasta), file.out = 'toAlign.fasta')
}

if(xargs$seqType == 'mitogenome'){
  fastas <- paste0('/shared/projects/messor_dils/messor_data/mitogenomes/',samples$V1,'_cleaned.mitof.fasta')
  
  bigFasta <- list()
  for (f in fastas){
    dna <- read.fasta(f, seqtype = c("DNA"),
                      seqonly = FALSE, whole.header = T,
                      forceDNAtolower = FALSE, 
                      set.attributes = TRUE)
    
    sample <- word(tools::file_path_sans_ext(f), sep ="/",-1)
    
  ids <- tibble(id = names(dna)) %>% 
    mutate(sample = word(id, sep = '@',1), 
           locus = word(id, sep = '@',2), 
           newID = paste0(locus,"_",sample," |",locus))
  
  names(dna) <- ids %>% pull(newID)
  bigFasta <- c(bigFasta, dna)
  }
  write.fasta(bigFasta, names = names(bigFasta), file.out = 'toAlign.fasta')
}

if(xargs$seqType == 'busco'){
  fastas <- list.files('/shared/projects/messor_dils/messor_data/busco/', pattern = ".fasta", full.names = T)
  fastas <- fastas[str_detect(fastas, paste0(samples$V1, collapse = "|")) & !str_detect(fastas,"fasta.")]
  
  bigFasta <- list()
  for (f in fastas){
    dna <- read.fasta(f, seqtype = c("DNA"),
                      seqonly = FALSE, whole.header = T,
                      forceDNAtolower = FALSE, 
                      set.attributes = TRUE)
    
    sample <- word(tools::file_path_sans_ext(f), sep ="/",-1)
    
    ids <- tibble(id = names(dna)) %>% 
      mutate(sample = sample, 
             locus = gsub("_pilon", "", word(id, sep = '@',2)), 
             newID = paste0(locus,"_", sample," |", locus))
    
    names(dna) <- ids %>% pull(newID)
    bigFasta <- c(bigFasta, dna)
  }
  write.fasta(bigFasta, names = names(bigFasta), file.out = 'toAlign.fasta')
}

if(xargs$seqType == 'ambiguous'){
  fastas <- paste0('/shared/projects/messor_dils/messor_data/ambiguous/',samples$V1,'_amb.fasta')
  
  bigFasta <- list()
  for (f in fastas){
    
    dna <- read.fasta(f, seqtype = c("DNA"),
                      seqonly = FALSE, whole.header = T,
                      forceDNAtolower = FALSE, 
                      set.attributes = TRUE)
    
    ids <- tibble(id = names(dna)) %>% 
      mutate(sample = word(id, sep = '@',1), 
             locus = word(id, sep = '@',2), 
             newID = paste0(locus,"_",sample," |",locus))
    
    names(dna) <- ids %>% pull(newID)
    bigFasta <- c(bigFasta, dna)
  }
  write.fasta(bigFasta, names = names(bigFasta), file.out = 'toAlign.fasta')
}

allFasta <- list.files(pattern="toAlign.fasta", recursive = TRUE, full.names = TRUE)

argsFinal <- c(paste('--input', allFasta),
               '--output aligned',
               '--output-format fasta',
               paste("--taxa", Nsamples),
               "--aligner mafft",
               paste('--cores',thr),
               '--incomplete-matrix',
               '--ambiguous')
system2(paste0(phyluce,'phyluce_align_seqcap_align'), argsFinal)


#####################
### for some reason the aligning script produces some labels with _R_ which affected the outcomes later so this R script removes loci names from sequences and clean this issue = cleanLociFromFasta.r 

fastas <- list.files('aligned', pattern = ".fasta",full.names = T)

foreach(f = fastas) %dopar% {
  dna <- read.fasta(f, seqtype = c("DNA"),
                    seqonly = FALSE, whole.header = T,
                    forceDNAtolower = FALSE, 
                    set.attributes = TRUE)
  locus <- word(tools::file_path_sans_ext(f), sep ="/",-1)
  names(dna) <- gsub(paste0(locus,"_"),"",names(dna))
  names(dna) <- gsub("_R_","",names(dna))
  test <- ifelse(any(str_detect(names(dna),"_R_")), locus, NA)
  if(!is.na(test))
    print(paste("check this locus:",test))
  write.fasta(dna, names = names(dna), file.out = f)
}

#####################

argsFinal <- c('--alignments aligned',
               '--input-format fasta',
               '--output alignedGblocks',
               '--output-format fasta',
               paste('--cores',thr))
system2(paste0(phyluce,'phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed'), argsFinal)

#####################

fastas <- list.files('alignedGblocks', pattern = ".fasta",full.names = T)

dt <- foreach(f = fastas, .combine=rbind) %dopar% {
  dna <- read.dna(f, "fasta", as.character = FALSE, as.matrix = TRUE)
  locus <- word(tools::file_path_sans_ext(f), sep ="/",-1)
  data.frame(aln = locus, 
             informative_sites = pis(dna,"absolute"),
             differences = length(seg.sites(dna)),
             length = ncol(dna),
             nind = nrow(dna),
             t(base.freq(dna, freq = TRUE, all = TRUE)[c('n','-')]))
}

dtout <- dt %>% 
  mutate(diffPp = differences/length,
         sitesPp = informative_sites/length,
         gapsNPp = (X. + n)/(length*nind)) %>%
  filter(informative_sites > 0,
         length > as.numeric(350), 
         gapsNPp < as.numeric(0.15)) %>%  ### bpMin = 350, gapsNPp = 0.15 to exclude alignments too short and with too much missing data
  pivot_longer(cols = c("length","diffPp", "sitesPp")) %>%
  group_by(name) %>%
  mutate(outlier = is_outlier(value)) %>% 
  group_by(aln) %>%
  summarize(any(outlier == TRUE)) ### exclude any outlier loci given the distribution of loci size, proportion of parsimonious and segregating sites

dt <- left_join(dt, dtout) %>% 
  rename(gap = X., outlier = `any(outlier == TRUE)` ) %>% 
  replace_na(list(outlier = "TRUE"))

write.csv(dt, "stats.csv", quote = F, row.names = F)

oldpath <- paste0(getwd(),'/alignedGblocks')
newpath <- paste0(getwd(),"/alignedFiltered/")
dir.create(newpath)
fastas <- list.files(oldpath, pattern = ".fasta", full.names = T)

foreach(f = fastas) %dopar% {
  dna <- read.fasta(f, seqtype = c("DNA"),
                    seqonly = FALSE, whole.header = T,
                    forceDNAtolower = FALSE, 
                    set.attributes = TRUE)
  
  locus <- word(tools::file_path_sans_ext(f), sep ="/",-1)
  if(length(dna) > Nsamples*as.numeric(xargs$Nfil) & 
     dt[dt$aln %in% locus,'outlier'] == FALSE){  ### e.g. Nfil=50% of samples included
    
    if(!xargs$seqType %in% "phasedH01"){
      names(dna) <- ifelse(str_detect(names(dna),"_cleaned.nucl.megahit.busco.nc|_busco"), 
                           word(names(dna),sep = "_",1,1),
                           word(names(dna), sep = fixed("."),1,1))
      }

    write.fasta(dna, names = names(dna), file.out = paste0(newpath,locus,".fasta"))
  }
}

dir.create('phyluce')
system("mv aligned* phyluce/")
system("mv *log phyluce/")
system("mv toAlign.fasta phyluce/")

argsFinal <- c('--alignments phyluce/alignedFiltered',
               '--input-format fasta',
               '--output concatenated',
               '--nexus')
system2(paste0(phyluce,'phyluce_align_concatenate_alignments'), argsFinal)

argsFinal <- c(paste('--alignments concatenated'),
               '--input-format nexus',
               '--output-format fasta',
               '--output concatenated1')
system2(paste0(phyluce,'phyluce_align_convert_one_align_to_another'), argsFinal)

system("mv concatenated1/*fasta .")
system("rm -rf concatenated")
system("rm -rf concatenated1")

#system2('mv', c(xargs$input, getwd()))

