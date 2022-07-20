# DADA2 Workflow
# To run
samnames <- 1 #if the file names are separated by underscores, where are the samplenames?
primer_fwd <- c("AAACTYAAAKGAATTGRCGG")
primer_rev <- c("ACGGGCGGTGTGTRC")
filter_chl_metazoa <- TRUE #do you want to take out chl & metazoa?


#==========Load Libraries==========
cat("Loading Libraries... \n")
library(dada2)
library(Biostrings) 
library(tidyverse)
source("/data/WINONA/common_bioinformatics/extra_functions.R")

#==========Make Directories==========
cat("Making Directories... \n")
dir_blast <- "./blast/"
dir_dada2 <- "./dada2/"
dir_database <- "/data/WINONA/common_bioinformatics/database/"
dir_fastq <- "./fastq_cutadapt/"
dir_filtered <- "./fastq_filtered/"
dir_quality <- "./quality/"

dir.create(dir_filtered)
dir.create(dir_quality)
# dir.create(dir_dada2)
# dir.create(dir_blast)

#==========Get files==========
cat("Getting all files... \n")
fns <- list.files(dir_fastq, full.names = TRUE)

cat("Getting reverse & forward files... \n")
fns_R1 <- fns[str_detect(basename(fns), "_1.")]
fns_R2 <- fns[str_detect(basename(fns), "_2.")]

sample_names <- str_split(basename(fns_R1), pattern = "_", simplify = TRUE)
sample_names <- sample_names[ , samnames]
saveRDS(sample_names, "00_sample_names.RDS")

#==========Check Primers==========
cat("Checking Primers... \n")
primer_length_fwd <- str_length(primer_fwd)
primer_length_rev <- str_length(primer_rev)

FWD.orients <- allOrients(primer_fwd)
REV.orients <- allOrients(primer_rev)
FWD.orients
REV.orients

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fns_R1[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fns_R2[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fns_R1[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fns_R2[[1]]))

#==========compute number of paired reads==========
cat("Computing number of reads in each file \n")
#create empty data frame
df <- data.frame()
#loop through all the R1 files (no need to go through R2 which should be the same)
for (i in 1:length(fns_R1)){
  geom <- fastq.geometry(fns_R1[i]) #get  number of reads
  df_one_row <- data.frame(n_seq = geom[1], file_name = sample_names[i]) #make nice one row
  df_one_row$file_name <- as.character(df_one_row$file_name)
  df <- bind_rows(df,df_one_row) #make dataframe
  saveRDS(df, "01_paired_reads.RDS")
}
#display number of sequences and write data to a small file
cat("Number of reads in each file: \n")
df


#===========Plot Quality Profile==========
cat("Plotting quality profile... \n")
for (i in 1:4){
  qp_1 <- plotQualityProfile(fns[i])
  qp_1_filename <- paste0(dir_quality, basename(fns[i]), "_qual.pdf")
  ggsave(plot = qp_1, filename = qp_1_filename, device = "pdf", width = 15, height = 10, scale = 1, units = "cm")
}

#==========Filter & Trim Reads==========
cat("Filtering and Trimming reads\n")
#create filenames for the filtered reads
filt_R1 <- str_c(dir_filtered, sample_names, "_R1_filt.fastq")
filt_R2 <- str_c(dir_filtered, sample_names, "_R2_filt.fastq")

#filter and trim!
cat("Filtering & Trimming reads... \n")
out8 <- filterAndTrim(fwd = fns_R1, filt = filt_R1, rev = fns_R2, filt.rev = filt_R2,
                      compress = FALSE, truncQ = 8,
                      maxEE = c(50,50), multithread = TRUE)
head(out8)
saveRDS(out8, "02_filter_trim_Q8.RDS")
cat("\n Reads filtered & trimmed. \n")

#==========Learn Error Rate==========
cat("Learning forward reads error... \n")
error_R1 <- learnErrors(filt_R1, multithread = TRUE, randomize = TRUE)
saveRDS(error_R1, file = "03_error_R1_q8.RDS")

cat("Learning reverse reads error... \n")
error_R2 <- learnErrors(filt_R2, multithread = TRUE, randomize = TRUE)
saveRDS(error_R2, file = "03_error_R2_q8.RDS")

plotErrors(error_R1, nominalQ = TRUE)
ggsave("03_error_R1_q8.png", device = "png", dpi = "retina", width = 15, height = 15, units = "cm")

plotErrors(error_R2, nominalQ = TRUE)
ggsave("03_error_R2_q8.png", device = "png", dpi = "retina", width = 15, height = 15, units = "cm")

cat("Errors plotted & saved. \n")


#========== Dereplication, Denoising, Merging
if(length(sample_names) > 20){
  cat("Big Data: dereplicating, denoising, and merging... \n")
  # define a function to track
  getN <- function(x) {sum(getUniques(x))}
  
  derep_track <- rep(0, length = length(sample_names))
  dada_track <- rep(0, length = length(sample_names))
  
  # Sample inference and merger of paired-end reads
  mergers <- vector("list", length(sample_names))
  names(mergers) <- sample_names
  
  for(i in 1:(length(sample_names))){
    sam <- sample_names[i]
    cat("Processing:", sam, "\n")
    
    derep_R1 <- derepFastq(filt_R1[i])
    dada_R1 <- dada(derep_R1, err=error_R1, multithread=TRUE)
    
    derep_R2 <- derepFastq(filt_R2[i])
    dada_R2 <- dada(derep_R2, err=error_R2, multithread=TRUE)
    
    derep_track[i] <- getN(derep_R1)
    dada_track[i] <- getN(dada_R1)
    
    merger <- mergePairs(dada_R1, derep_R1, dada_R2, derep_R2, verbose = TRUE)
    mergers[[sam]] <- merger
  }
  rm(derep_R1, derep_R2)
  saveRDS(mergers, file = "06_mergers.RDS")
  cat("Merged files saved. \n")
  
  cat("\n Tracking number of reads: \n")
  track <- cbind(out8, derep_track, dada_track, sapply(mergers, getN))
  
  colnames(track) <- c("input", "filtered", "dereplicated", "denoised", "merged")
  rownames(track) <- sample_names
  
  print(track)
  write_tsv(data.frame(track), "09_track_numbers_part1.tsv")
  
  
}else{
  cat("Dereplicating... \n")
  derep_R1 <- derepFastq(filt_R1, verbose = TRUE)
  derep_R2 <- derepFastq(filt_R2, verbose = TRUE)
  names(derep_R1) <- sample_names
  names(derep_R2) <- sample_names
  
  saveRDS(derep_R1, file = "04_derep_R1.RDS")
  saveRDS(derep_R2, file = "04_derep_R2.RDS")
  
  cat("Denoising... \n")
  dada_R1 <- dada(derep_R1, err = error_R1, multithread = TRUE, pool = FALSE)
  dada_R2 <- dada(derep_R2, err = error_R2, multithread = TRUE, pool = FALSE)
  
  saveRDS(dada_R1, file = "05_dada_R1.RDS")
  saveRDS(dada_R2, file = "05_dada_R2.RDS")
  
  dada_R1[[1]]
  dada_R2[[1]]
  
  cat("Merging... \n")
  mergers <- mergePairs(dada_R1, derep_R1, dada_R2, derep_R2, verbose = TRUE)
  saveRDS(mergers, file = "06_mergers.RDS")
  head(mergers[[1]])
  cat("\n Merged files saved. \n")
  
  getN <- function(x) {sum(getUniques(x))}
  
  cat("\n Tracking number of reads: \n")
  track <- cbind(out8, sapply(derep_R1, getN), sapply(dada_R1, getN), sapply(mergers, getN))
  
  colnames(track) <- c("input", "filtered", "dereplicated", "denoised", "merged")
  rownames(track) <- sample_names
  
  print(track)
  write_tsv(data.frame(track), "09_track_numbers_part1.tsv")
  
}

sessionInfo()

