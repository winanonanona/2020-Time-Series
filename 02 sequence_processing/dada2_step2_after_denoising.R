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
# dir_fastq <- "./fastq_cutadapt/"
# dir_filtered <- "./fastq_filtered/"
# dir_quality <- "./quality/"

# dir.create(dir_filtered)
# dir.create(dir_quality)
dir.create(dir_dada2)
dir.create(dir_blast)

cat("\n Loading stuff you need... \n")
sample_namesv1 <- readRDS("./v1/00_sample_names.RDS")
sample_namesv2 <- readRDS("./v2/00_sample_names.RDS")
sample_names <- c(sample_namesv1, sample_namesv2)

cat("\n Combining sequence table... \n")
mergedv1 <- readRDS("./v1/06_mergers.RDS")
mergedv2 <- readRDS("./v2/06_mergers.RDS")
mergers <- c(mergedv1, mergedv2)

#==========Make Sequence Table!==========
cat("\n Making sequence table... \n")
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "07_seqtab.RDS")

#inspect distribution of sequence lengths
dim(seqtab)
table(nchar(getSequences(seqtab)))
plot(table(nchar(getSequences(seqtab))))
ggsave("07_sequence_length.png", device = "png", dpi = "retina", width = 15, height = 12, units = "cm")

#==========Removing Chimeras==========
cat("\n Removing Chimeras... \n")
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "pooled", multithread = TRUE, verbose = TRUE)
saveRDS(seqtab_nochim, "08_seqtab_nochim.RDS")

#compute percentage of non-chimeras
cat("\n percentage of non-chimeras: ", sum(seqtab_nochim)/sum(seqtab)*100, "%", sep = "")
cat("\n total number of sequences: ", sum(seqtab_nochim), sep = "")

#==========Track Number of Reads==========
getN <- function(x) {sum(getUniques(x))}

cat("\n Tracking number of reads: \n")
track <- cbind(sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab_nochim))

colnames(track) <- c("merged", "tabled", "no chimeras")
rownames(track) <- sample_names

print(track)
write_tsv(data.frame(track), "09_track_numbers_part2.tsv")


#==========Transforming & Saving ASVs==========
cat("\n Transforming & Saving ASVs... \n")
seqtab_nochim_transpose <- as.data.frame(t(seqtab_nochim), stringsAsFactors = FALSE) %>%
  rownames_to_column(var = "sequence") %>%
  rowid_to_column(var = "ASVNumber") %>%
  mutate(ASVNumber = sprintf("asv%04d", ASVNumber)) %>%
  mutate(sequence = str_replace_all(sequence, "(-|\\.)", ""))
saveRDS(seqtab_nochim_transpose, "10_seqtab_nochim_transpose.RDS")

df <- seqtab_nochim_transpose
seq_out <- Biostrings::DNAStringSet(df$sequence)

names(seq_out) <- df$ASVNumber

Biostrings::writeXStringSet(seq_out, str_c(dir_dada2, "10_ASV_no_taxa.fasta"), compress = FALSE, width = 20000)

#==========Removing stuff I don't need==========
rm(dada_R1, dada_R2, derep_R1, derep_R2, df, df_one_row, error_R1, error_R2, merger, mergers, out8, qp_1, seq_out, track, z, dada_track, derep_track, geom, i, qp_1_filename, sam, seqtab)

#==========Assigning Taxonomy!==========
##=========Silva=========
cat("\n Assigning Taxonomy with Silva... \n")
silva_dat <- paste0(dir_database,"silva_nr_v132_train_set.fa.gz")
silva_taxa <- assignTaxonomy(seqtab_nochim, refFasta = silva_dat, outputBootstraps = TRUE, verbose = TRUE, multithread = TRUE)
saveRDS(silva_taxa, file = str_c("11_silva_taxa.RDS"))

write_tsv(as_tibble(silva_taxa$tax), file = str_c(dir_dada2, "silva_taxa.txt"))
write_tsv(as_tibble(silva_taxa$boot), file = str_c(dir_dada2, "silva_taxa_boot.txt"))

silva_taxa_tax <- as.data.frame(silva_taxa$tax, stringsAsFactors = FALSE)
silva_taxa_boot <- as.data.frame(silva_taxa$boot, stringsAsFactors = FALSE)
columns <- colnames(silva_taxa_boot)
colnames(silva_taxa_boot) <- str_c(columns,"_boot")

seqtab_nochim_transpose_silva <- silva_taxa_tax %>%
  bind_cols(silva_taxa_boot) %>%
  bind_cols(seqtab_nochim_transpose)

write_csv(seqtab_nochim_transpose_silva, paste0(dir_dada2, "12_silva_seqtab_taxa.csv"))
saveRDS(seqtab_nochim_transpose_silva, paste0("12_silva_seqtab_taxa.RDS"))

rownames(seqtab_nochim_transpose_silva) <- NULL

##=========Removing copepods & chloroplasts=========
if(filter_chl_metazoa == TRUE){
  ###find all chloroplast & metazoa
  index_arthropods <- which(seqtab_nochim_transpose_silva$Phylum == "Arthropoda")
  index_chloroplast <- which(seqtab_nochim_transpose_silva$Order == "Chloroplast")
  cat("Number of arthropods:", length(index_arthropods))
  cat("Number of chloroplast:", length(index_chloroplast))
  
  ###take them out
  seqtab_nochim_transpose_silva_filtered <- seqtab_nochim_transpose_silva[-c(index_chloroplast, index_arthropods),]
  cat("\n Arthropods & Chloroplast taken out \n")
}else{
  seqtab_nochim_transpose_silva_filtered <- seqtab_nochim_transpose_silva
}

###find all remaining euks & unassigned
index_euk <- which(seqtab_nochim_transpose_silva_filtered$Kingdom == "Eukaryota")
index_unassigned <- which(is.na(seqtab_nochim_transpose_silva_filtered$Kingdom))
cat("\n Number of remaining eukaryotes:", length(index_euk))
cat("\n Number of unassigned sequences:", length(index_unassigned))

###take them out
seqtab_nochim_transpose_silva_filtered_noeuk <- seqtab_nochim_transpose_silva_filtered[-c(index_euk, index_unassigned), ]
cat("\n Remaining euks & unassigned taken out \n")

###prepare for reassignment
#get remaining euks & unassigned
seqtab_nochim_transpose_silva_filtered_euks <- seqtab_nochim_transpose_silva_filtered[c(index_euk, index_unassigned), ]
to_reassign <- seqtab_nochim_transpose_silva_filtered_euks$sequence

# seqtab for assigntaxonomy
seqtab_nochim_reassign <- seqtab_nochim[ , (which(colnames(seqtab_nochim) %in% to_reassign))]

#seqtab for bind_cols
seqtab_nochim_reassign_transpose <- seqtab_nochim_transpose %>%
  filter(sequence %in% to_reassign)

##=========PR2=========
cat("Assigning remaining Eukaryotes & unassigned sequences with PR2... \n")
PR2_tax_levels <- c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")
PR2_dat <- paste0(dir_database,"pr2_version_4.13.0_18S_dada2.fasta.gz")
PR2_taxa <- assignTaxonomy(seqtab_nochim_reassign, refFasta = PR2_dat, taxLevels = PR2_tax_levels, outputBootstraps = TRUE, verbose = TRUE, multithread = TRUE)
saveRDS(PR2_taxa, file = str_c("13_PR2_18s_taxa.RDS"))

write_tsv(as_tibble(PR2_taxa$tax), path = str_c(dir_dada2, "PR2_taxa.txt"))
write_tsv(as_tibble(PR2_taxa$boot), path = str_c(dir_dada2, "PR2_taxa_boot.txt"))

PR2_taxa_tax <- as.data.frame(PR2_taxa$tax, stringsAsFactors = FALSE)
PR2_taxa_boot <- as.data.frame(PR2_taxa$boot, stringsAsFactors = FALSE) %>%
  rename_all(funs(str_c(.,"_boot")))

seqtab_nochim_transpose_PR2 <- PR2_taxa_tax %>%
  bind_cols(PR2_taxa_boot) %>%
  bind_cols(seqtab_nochim_reassign_transpose)

colnames(seqtab_nochim_transpose_PR2)

write_csv(seqtab_nochim_transpose_PR2, paste0(dir_dada2, "14_PR2_seqtab.csv"))
saveRDS(seqtab_nochim_transpose_PR2, paste0("14_PR2_seqtab_taxa.RDS"))

rownames(seqtab_nochim_transpose_PR2) <- NULL

##==========Stick Both Together==========
#add supergroup/division to silva-ed
cat("\n Sticking both together... \n")
seqtab_nochim_transpose_silva_edited <- seqtab_nochim_transpose_silva_filtered_noeuk %>%
  mutate(Supergroup = NA, .after = Kingdom) %>%
  mutate(Species = NA, .after = Genus) %>%
  mutate(Supergroup_boot = Kingdom_boot, .after = Kingdom_boot) %>%
  mutate(Species_boot = Genus_boot, .after = Genus_boot)

seqtab_nochim_transpose_PR2_edited <- seqtab_nochim_transpose_PR2 %>%
  rename(c("Phylum" = "Division", "Phylum_boot" = "Division_boot"))

# STICK THEM TOGETHER WOOHOO
seqtab_nochim_transpose_combined <- seqtab_nochim_transpose_PR2_edited %>%
  bind_rows(seqtab_nochim_transpose_silva_edited) %>%
  arrange(ASVNumber)
saveRDS(seqtab_nochim_transpose_combined, "15_seqtab_combined.RDS")

cat("\n Total Number of ASVs:", nrow(seqtab_nochim_transpose_combined))

#==========Removing unimportant ASVs============
cat("\n Removing ASVs with less than 10 reads across all samples... \n")
# less than 10 features
seqtab_final <- seqtab_nochim_transpose_combined %>%
  mutate(features.sum = rowSums(.[19:(length(sample_names)+18)])) %>%
  filter(features.sum > 9) %>%
  mutate(features.sum = NULL)

saveRDS(seqtab_final, "16_seqtab_final.RDS")
write_csv(seqtab_final, paste0(dir_dada2, "16_seqtab_final.csv"))
cat("\n Number of ASVs retained:", nrow(seqtab_final))


index_shit <- which(nchar(seqtab_final[,18]) < 400)
for(i in index_shit){
  print(sum(seqtab_final[i,19:(length(sample_names)+18)]))
}

rownames(seqtab_final) <- NULL

#==========Remove NAs to make pretty phyloseq objects==========
cat("\n Removing NAs to make pretty phyloseq objects... \n")
seqtab_final_noNAs <- seqtab_noNAs(seqtab_final,taxa_levels = "pr2")
saveRDS(seqtab_final_noNAs, "17_seqtab_final_noNAs.RDS")
write_csv(seqtab_final_noNAs, paste0(dir_dada2, "17_seqtab_final_noNAs.csv"))

#==========Remove Phylum bootstrap <80 to make nicer bar plots==========
cat("\n Removing Phylum bootstrap <80 to make nicer bar plots... \n")
seqtab_final_phy80 <- seqtab_final_noNAs %>%
  filter(Phylum_boot >= 80)
saveRDS(seqtab_final_phy80, "17_seqtab_final_phy80.RDS")
write_csv(seqtab_final_phy80, paste0(dir_dada2, "17_seqtab_final_phy80.csv"))

sessionInfo()
