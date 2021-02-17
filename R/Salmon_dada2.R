############################################
### V4 processing for Salmon2
### Lou LaMartina, started Feb 15, 2020
############################################


library(dada2)
library(phyloseq)
library(decontam)
library(ggplot2)


###################################
### prepare data for processing ###
###################################

# set working directory
setwd("~/Desktop/Lab/Projects/Aeromonas/Seafood/Salmon2")


# set file paths
path <- "./cutadapt"
pathF <- "./cutadapt/fastqF"
pathR <- "./cutadapt/fastqR"


# set paths for filtered reads
filtered_pathF <- "./cutadapt/fastqF/Filtered"
filtered_pathR <- "./cutadapt/fastqR/Filtered"


# sort forward and reverse reads into file paths
fastqFs <- sort(list.files(pathF, pattern = "_pass_1.fastq.gz", full.names = TRUE))
fastqRs <- sort(list.files(pathR, pattern = "_pass_2.fastq.gz", full.names = TRUE))


# extract file names
Sample_names <- sapply(strsplit(basename(fastqFs), "_"), '[', 1)




################################
### Inspect sequence quality ###
################################

# visualize quality of reads
qualityF <- plotQualityProfile(fastqFs[1:4]); qualityF
qualityR <- plotQualityProfile(fastqRs[1:4]); qualityR


# plot quality profiles of filtered reads
filtF.plot <- plotQualityProfile(filtered_Fs[1:4]); filtF.plot
filtR.plot <- plotQualityProfile(filtered_Rs[1:4]); filtR.plot


# save quality profiles
ggsave("./Plots/filt_qualityF.pdf", plot = filtF.plot, device = "pdf", width = 12, height = 8, units = "in")
ggsave("./Plots/filt_qualityR.pdf", plot = filtR.plot, device = "pdf", width = 12, height = 8, units = "in")


# set sample names to the ID only
names(filteredFs) <- Sample_names
names(filteredRs) <- Sample_names




#######################
### Quality control ###
#######################

# give filtered files new names and paths
filteredFs <- file.path(filtered_pathF, paste0(Sample_names, "_F_filt.fastq.gz"))
filteredRs <- file.path(filtered_pathR, paste0(Sample_names, "_R_filt.fastq.gz"))


# filter based on quality and read length
filtered_out <- filterAndTrim(fastqFs, filteredFs, fastqRs, filteredRs, 
                              maxEE = 2, maxN = 0, truncQ = 10,
                              rm.phix = TRUE, compress = TRUE, verbose = TRUE, multithread = TRUE)


# inspect how many reads were filtered out of each sample
filtered_out
(1 - (filtered_out[,2] / filtered_out[,1])) * 100
mean((1 - (filtered_out[,2] / filtered_out[,1])) * 100)
# [1] 1.600213 % removed


# plot quality profiles of filtered reads
filtF.plot <- plotQualityProfile(filtered_Fs[1:4]); filtF.plot
filtR.plot <- plotQualityProfile(filtered_Rs[1:4]); filtR.plot


# save quality profiles
ggsave("./Plots/filt_qualityF.pdf", plot = filtF.plot, device = "pdf", width = 12, height = 8, units = "in")
ggsave("./Plots/filt_qualityR.pdf", plot = filtR.plot, device = "pdf", width = 12, height = 8, units = "in")


# set sample names to the ID only
names(filteredFs) <- Sample_names
names(filteredRs) <- Sample_names




############################
### Learning error rates ###
############################

# learn and visualize error rates of F reads
errorF <- learnErrors(filtered_s, multithread = TRUE)
errF.plot <- plotErrors(errorF, nominalQ = TRUE, err_in = TRUE, err_out = TRUE); errorF


# learn and visualize error rates of R reads
errorR <- learnErrors(filteredRs, multithread = TRUE)
errR.plot <- plotErrors(errorR, nominalQ = TRUE, err_in = TRUE, err_out = TRUE); errorR


# save error plots
ggsave("./Plots/errorF.pdf", plot = errF.plot, device = "pdf", width = 12, height = 8, units = "in")
ggsave("./Plots/errorR.pdf", plot = errR.plot, device = "pdf", width = 12, height = 8, units = "in")




################################
### Merging paired-end reads ###
################################

# create list of merged reads
Mergers <- vector("list", length(Sample_names))
names(Mergers) <- Sample_names


# sample inference and merging paired-end reads
for(i in Sample_names) {
  cat("\nProcessing", i, "(", match(i, Sample_names),"/", 
      length(Sample_names), ") :", format(Sys.time(), "%H:%M %p"), "\n")
  derepF <- derepFastq(filteredFs[[i]])
  dadaF <- dada(derepF, err = errorF, multithread = TRUE)
  derepR <- derepFastq(filteredRs[[i]])
  dadaR <- dada(derepR, err = errorR, multithread = TRUE)
  Merger <- mergePairs(dadaF, derepF, dadaR, derepR)
  Mergers[[i]] <- Merger
}


# removing dereps to save memory
rm(derepF, derepR)


# contruct a sequence table
Sequence_table <- makeSequenceTable(Mergers)


# dimensions: num rows x num columns
dim(Sequence_table)
# [1]   23 2666




##################################
### Quality control: processed ###
##################################


########
### trim

# inspect sequence distribution
seq_distribution <- data.frame(table(nchar(getSequences(Sequence_table)))) 
colnames(seq_distribution) <- c("SeqLength", "Frequency")
seq_distribution$SeqLength <- as.numeric(as.character(seq_distribution$SeqLength))


# visualize and save
dist <- 
  ggplot(seq_distribution, aes(x = SeqLength, y = Frequency)) +
  geom_bar(stat = "identity", width = 0.5, fill = "black") +
  geom_point(aes(x = median(nchar(getSequences(Sequence_table))), 
                 y = max(seq_distribution$Frequency)), color = "red") +
  geom_text(data = data.frame(x = median(nchar(getSequences(Sequence_table))), 
                              y = max(seq_distribution$Frequency)), 
            aes(x, y, label = "median"), vjust = 2, color = "red") +
  scale_x_continuous(breaks = ceiling(quantile(seq_distribution$SeqLength)),
                     labels = ceiling(quantile(seq_distribution$SeqLength)))  +
  labs(x = "Sequence length (bp)", y = "Frequency")
dist

ggsave("seq_distribution.pdf", plot = dist, device = "pdf", width = 16, height = 4, units = "in")


# remove reads of non target length, 5% above and below the median 
median(nchar(getSequences(Sequence_table)))
# [1] 284

min_len <- floor(median(nchar(getSequences(Sequence_table))) * 0.95)
max_len <- ceiling(median(nchar(getSequences(Sequence_table))) * 1.05)


# modify sequence table with new guidelines
Seq_table_trimmed <- Sequence_table[ ,nchar(colnames(Sequence_table)) 
                                     %in% seq(min_len, max_len)]


dim(Seq_table_trimmed)
#[1]   23 2281



###################
### Remove chimeras

# removing chimeras with denovo screening
Sequence_table_nochim <- removeBimeraDenovo(Seq_table_trimmed, method = "consensus",
                                            multithread = TRUE, verbose = TRUE)


# how many unique sequences were moved?
ncol(Sequence_table)               # [1] 2666
ncol(Sequence_table_nochim)        # [1] 1758
ncol(Sequence_table) - 
       ncol(Sequence_table_nochim) # [1] 908


# what percentage of reads were identified as chimeras?
(1 - sum(Sequence_table_nochim) / sum(Sequence_table)) * 100
# [1] 89.99269 % 




#######################
### Assign taxonomy ###
#######################

Sequence_taxa_table <- assignTaxonomy(Sequence_table_nochim, 
                                      "~/Desktop/Lab/Projects/Misc/silva_nr_v132_train_set.fa.gz", 
                                      multithread = TRUE)


# remove nonspecific ASVs
nonspecifics <- rbind(subset(data.frame(Sequence_taxa_table), Kingdom == "Eukaryota"),
                      subset(data.frame(Sequence_taxa_table), Order == "Chloroplast"),
                      subset(data.frame(Sequence_taxa_table), Family == "Mitochondria"))

Sequence_taxa_table <- Sequence_taxa_table[! rownames(Sequence_taxa_table) %in% rownames(nonspecifics),]
Sequence_table_nochim <- Sequence_table_nochim[, ! colnames(Sequence_table_nochim) %in% rownames(nonspecifics)]
Sequence_table_relabun <- Sequence_table_nochim / rowSums(Sequence_table_nochim)




############
### save ###
############

write.csv(data.frame(Sample_name = rownames(Sequence_table_nochim), Sequence_table_nochim), 
          "Salmon_counts.csv", row.names = FALSE)
write.csv(data.frame(FASTA = rownames(Sequence_taxa_table), Sequence_taxa_table), 
          "Salmon_taxa.csv", row.names = FALSE, na = "")
write.csv(data.frame(Sample_name = rownames(Sequence_table_relabun), Sequence_table_relabun), 
          "Salmon_relabun.csv", row.names = FALSE)




# # # # # # # # # # # # # # # # # # # # # # 
save.image("Salmon_dada2_environment.RData")
# # # # # # # # # # # # # # # # # # # # # # 


