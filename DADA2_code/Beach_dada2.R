#######################################
### V4V5 processing of beach samples
### samples from Beattie et al 2020
### Lou LaMartina, started Jan 27, 2020
#######################################



######################
### load libraries ###
######################

library(dada2)
library(phyloseq)
library(decontam)
library(ggplot2)



###################################
### prepare data for processing ###
###################################

# set working directory
setwd("~/Desktop/Aeromonas/Beach")


# set file paths
path <- "./cutadapt"
pathF <- "./cutadapt/fastqF"
pathR <- "./cutadapt/fastqR"


# set paths for filtered reads
filtered_pathF <- "./cutadapt/fastqF/Filtered"
filtered_pathR <- "./cutadapt/fastqR/Filtered"


# sort forward and reverse reads into file paths
fastqFs <- sort(list.files(pathF, pattern = "R1.fastq.gz", full.names = TRUE))
fastqRs <- sort(list.files(pathR, pattern = "R2.fastq.gz", full.names = TRUE))


# extract file names
Sample_names <- sub("_R1.fastq.gz", "", sub("1St_", "", basename(fastqFs)))




################################
### Inspect sequence quality ###
################################

# visualize quality of reads
qualityF <- plotQualityProfile(fastqFs[1:4]); qualityF
qualityR <- plotQualityProfile(fastqRs[1:4]); qualityR


# save quality profiles
ggsave("./Plots/qualityF.pdf", plot = qualityF, device = "pdf", width = 12, height = 8, units = "in")
ggsave("./Plots/qualityR.pdf", plot = qualityR, device = "pdf", width = 12, height = 8, units = "in")


# check: if there are not the same number of F and R files, stop.
if(length(fastqFs) != length(fastqRs)) 
  stop("Forward and reverse files do not match.")




#######################
### Quality control ###
#######################

# give filtered files new names and paths
filteredFs <- file.path(filtered_pathF, paste0(Sample_names, "_F_filt.fastq"))
filteredRs <- file.path(filtered_pathR, paste0(Sample_names, "_R_filt.fastq"))


# # # # # # # # # # # # # # # # # # # # # # # # # #
### FOR LARGE DATASETS (>100 files)

# what is the maximum number of files you want to 
# filterAndTrim at the same time?
# for example, if you have 110 files, 
# and you want to process only 20 at a time,
# we will process five sets of 20 and
# one set of 10.
# put your answer here:

maxN = 20

# then run the code below.
# # # # # # # # # # # # # # # # # # # # # # # # # # 


# define variables
intN = length(fastqFs) %/% maxN
remN = length(fastqFs) %% maxN



##############
### unfiltered

# create sets of files
fastqFsets.ls <- list()

for (i in c(seq(length(fastqFs) %/% maxN), intN + 1)) {
  fastqFsets.ls[[i]] <- fastqFs[(maxN * (i - 1) + 1) : (i * maxN)]
  fastqFsets.ls[[intN + 1]] <- fastqFs[(maxN * intN + 1) : length(fastqFs)]
}


# do same for Rs by changing file names
fastqRsets.ls <- list()
for (f in 1:length(fastqFsets.ls)) {
  fastqRsets.ls[[f]] <- sub("fastqF", "fastqR", sub("_R1", "_R2", fastqFsets.ls[[f]]))
}



############
### filtered

# create sets of files
filtFsets.ls <- list()
for(i in c(seq(length(filtered_Fs) %/% maxN), intN + 1)) {
  filtFsets.ls[[i]] <- filtered_Fs[(maxN * (i - 1) + 1) : (i * maxN)]
  filtFsets.ls[[intN + 1]] <- filtered_Fs[(maxN * intN + 1) : length(filtered_Fs)]
}


# do same for Rs by changing file names
filtRsets.ls <- list()
for(f in 1:length(filtFsets.ls)){
  filtRsets.ls[[f]] <- sub("fastqF", "fastqR", sub("F_filt", "R_filt", filtFsets.ls[[f]]))
}



###################
### filter and trim

filtered_out.ls <- list()
for(j in 1:length(fastqFsets.ls)) {
  cat("\n", format(Sys.time(), "%H:%M %p"), " - Filtering set ", j, " . . . \n")
  filtered_out.ls[[j]] <- filterAndTrim(fastqFsets.ls[[j]], filtFsets.ls[[j]],
                                        fastqRsets.ls[[j]], filtRsets.ls[[j]],
                                        maxEE = 2, maxN = 0, truncQ = 10,
                                        rm.phix = TRUE, compress = TRUE, 
                                        verbose = TRUE, multithread = TRUE)
}


# inspect how many reads were filtered out of each sample
filtered_out <- do.call(rbind, filtered_out.ls)
(1 - (filtered_out[,2] / filtered_out[,1])) * 100
mean((1 - (filtered_out[,2] / filtered_out[,1])) * 100)
# [1] 23.16235 % , kind of a lot
# 29.16279 mins with this method :)


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
# [1]   132 60655




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
# [1] 374 perfect!

min_len <- floor(median(nchar(getSequences(Sequence_table))) * 0.95)
max_len <- ceiling(median(nchar(getSequences(Sequence_table))) * 1.05)


# modify sequence table with new guidelines
Seq_table_trimmed <- Sequence_table[ ,nchar(colnames(Sequence_table)) 
                                     %in% seq(min_len, max_len)]


dim(Seq_table_trimmed)
# [1]   132 59889



###################
### Remove chimeras

# removing chimeras with denovo screening
Sequence_table_nochim <- removeBimeraDenovo(Seq_table_trimmed, method = "consensus",
                                            multithread = TRUE, verbose = TRUE)


# how many unique sequences were moved?
ncol(Sequence_table) - ncol(Sequence_table_nochim)
# [1] 16336


# what percentage of reads were identified as chimeras?
(1 - sum(Sequence_table_nochim) / sum(Sequence_table)) * 100
# [1] 2.293696 %




#######################
### Assign taxonomy ###
#######################

Sequence_taxa_table <- assignTaxonomy(Sequence_table_nochim, 
                                      "~/Desktop/Lab/Projects/Misc/silva_nr_v138_train_set.fa.gz", 
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
          "Beach_counts.csv", row.names = FALSE)
write.csv(data.frame(FASTA = rownames(Sequence_taxa_table), Sequence_taxa_table), 
          "Beach_taxa.csv", row.names = FALSE, na = "")
write.csv(data.frame(Sample_name = rownames(Sequence_table_relabun), Sequence_table_relabun), 
          "Beach_relabun.csv", row.names = FALSE)




# # # # # # # # # # # # # # # # # # # # # # 
save.image("Beach_dada2_environment.RData")
# # # # # # # # # # # # # # # # # # # # # # 



