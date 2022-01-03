#########################################
### Aeromonas 16S community meta-analysis
### Lou LaMartina, finalized Jan 3, 2022
#########################################
  

library(RColorBrewer)
library(reshape2)
setwd("~/Desktop/aero")



#######################
### getting started ###
#######################


# ~ or ~ load all files list - 3.971166 secs (RECOMMENDED)
all_files.ls <- readRDS("./RData/All_files_list.RData")


# subset ASV counts, relative abundances, taxonomy, and metadata files
counts_files.ls <- all_files.ls[grep("counts", names(all_files.ls))]
taxa_files.ls <- all_files.ls[grep("taxa", names(all_files.ls))]


# counts -
# change row names to sample names,
# remove sample name column
for (i in names(counts_files.ls)) {
  rownames(counts_files.ls[[i]]) <- counts_files.ls[[i]]$Sample_name
  counts_files.ls[[i]] <- counts_files.ls[[i]][-1]
}


# taxonomy -
# change row names to FASTA sequences
for (i in names(taxa_files.ls)) {
  rownames(taxa_files.ls[[i]]) <- taxa_files.ls[[i]]$FASTA
}


# remove "all files" list to save memory
rm(all_files.ls)


# load all info
info_all <- read.csv("./RData/Aeromonas_info_ALL.csv")


# get file names
fileNames_all <- list()
for(i in 1:length(counts_files.ls)){
  fileNames_all[[i]] <- rownames(counts_files.ls[[i]])
}
fileNames_all <- unlist(fileNames_all)


# get file names from sample info
fileNames_info <- as.character(info_all$File)


# subset sample info to files i have
info_all <- subset(info_all, File %in% fileNames_all)
fileNames_info <- as.character(info_all$File)


# keep files that i have info for
filt_files.ls <- list()
for(i in 1:length(counts_files.ls)){
  filt_files.ls[[i]] <- counts_files.ls[[i]][rownames(counts_files.ls[[i]]) %in% fileNames_info,]
}


fileNames_filt <- list()
for(i in 1:length(filt_files.ls)){
  fileNames_filt[[i]] <- rownames(filt_files.ls[[i]])
}
fileNames_filt <- unlist(fileNames_filt)


# still a difference of 10
fileNames_good <- intersect(fileNames_filt, fileNames_info)


# keep only those files
for(i in 1:length(filt_files.ls)){
  filt_files.ls[[i]] <- filt_files.ls[[i]][rownames(filt_files.ls[[i]]) %in% fileNames_good,]
}
fileNames_filt <- list()
for(i in 1:length(filt_files.ls)){
  fileNames_filt[[i]] <- rownames(filt_files.ls[[i]])
}
fileNames_filt <- unlist(fileNames_filt)


counts_files.ls <- filt_files.ls


# convert to relative abundance
relabun_files.ls <- list()
for(i in 1:length(counts_files.ls)){
  relabun_files.ls[[i]] <- counts_files.ls[[i]] / rowSums(counts_files.ls[[i]])
}




########################
### subset aeromonas ###
########################

aero_FASTAs.ls <- list()
aero_relabun.ls <- list()

for (i in 1:length(relabun_files.ls)) {
  aero_FASTAs.ls[[i]] <- subset(taxa_files.ls[[i]], Genus == "Aeromonas")$FASTA
  
  if (length(aero_FASTAs.ls[[i]]) > 1) {
    aero_relabun.ls[[i]] <- relabun_files.ls[[i]][, colnames(relabun_files.ls[[i]]) %in% aero_FASTAs.ls[[i]]]
  
  } else {
    aero_relabun.ls[[i]] <- data.frame(relabun_files.ls[[i]][, colnames(relabun_files.ls[[i]]) %in% aero_FASTAs.ls[[i]]])
    colnames(aero_relabun.ls[[i]]) <- aero_FASTAs.ls[[i]]
    rownames(aero_relabun.ls[[i]]) <- rownames(relabun_files.ls[[i]])
    
  }
}




##################
### merge data ###
##################

# transpose so FASTA are row names, add source
aero_relabun.t.ls <- list()

for (i in 1:length(aero_relabun.ls)) {
  aero_relabun.t.ls[[i]] <- data.frame(t(aero_relabun.ls[[i]]))
  aero_relabun.t.ls[[i]]$FASTA <- rownames(aero_relabun.t.ls[[i]])
}



###################
### merge by FASTAs

merged_relabun.ls <- list()

for (i in 1:length(aero_relabun.t.ls)) {
  if (i == 1) {
    merged_relabun.ls[[i]] <- merge(aero_relabun.t.ls[[i]], 
                                  aero_relabun.t.ls[[i + 1]], by = "FASTA", all = TRUE)
  } else if (i > 1 & i < length(aero_relabun.t.ls)) {
    merged_relabun.ls[[i]] <- merge(merged_relabun.ls[[i - 1]], 
                                  aero_relabun.t.ls[[i + 1]], by = "FASTA", all = TRUE)
  }
}



# extract the last one - merge() makes you do it one by one
merged_relabun <- data.frame(merged_relabun.ls[[length(aero_relabun.t.ls) - 1]])


# change "NA" to zero
merged_relabun[is.na(merged_relabun)] <- 0


# make FASTA row names, transpose
rownames(merged_relabun) <- merged_relabun$FASTA
merged_relabun <- data.frame(t(merged_relabun[-1]))


# remove empty ASVs
dim(merged_relabun) # [1] 251 152
relabun_filt <- merged_relabun[, colSums(merged_relabun) > 0]
dim(relabun_filt) # [1] 251 148




###################
### proportions ###
###################

# calculate proportion aeromonas in each sample
propAero <- relabun_filt
propAero$totalAero <- rowSums(propAero)
propAero <- data.frame(File = rownames(propAero), totalAero = propAero$totalAero)


# add sample info
propAero <- merge(propAero, info_all, by = "File")
propAero$Biome <- "Seafood"
propAero$Biome[propAero$SampleType %in% c("Sand", "Water")] <- "Environmental"



# labels
labs <- as.character(unique(propAero$SampleSource))
names(labs) <- c("Yellow\nPerch", "Atlantic\nSalmon", "Nase", "Backshore",
              "Berm", "Submerged", "Lake", "River", 
              "Post-\nchlorinated\nwastewater", "Pre-\nchlorinated\nwastewater",
              "Brown\nTrout", "Tench")

propAero$labs <- NA
for(i in 1:length(labs)){
  propAero$labs[propAero$SampleSource == labs[[i]]] <- names(labs)[[i]]
}


# colors
colors <- brewer.pal(11, "Spectral")[c(1,2,11,3,10)]


# plot
dot.plot <-
  ggplot(propAero, aes(x = labs, y = totalAero, 
                       color = SampleType, fill = SampleType)) +
  geom_boxplot(width = 0.4, size = 0.5, outlier.color = NA, alpha = 0.25, show.legend = F) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.5, show.legend = F) +
  facet_grid(. ~ Biome + SampleType, scales = "free", space = "free") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  theme_classic() +
  scale_y_continuous(trans = scales::pseudo_log_trans(0.001, 10),
                     breaks = c(0, 0.01, 0.1, 1), 
                     labels = c("0%", "1%", "10%", "100%")) +
  theme(axis.text.x = element_text(size = 6, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 8, color = "black", face = "bold"),
        axis.title.x = element_text(size = 8, color = "black", face = "bold"),
        strip.text = element_text(size = 8, color = "black", face = "bold"),
        strip.background = element_rect(colour = NA, fill = NA),
        panel.border = element_rect(size = 0.25, color = "grey80", fill = NA),
        axis.ticks = element_line(size = 0.25),
        axis.line = element_line(size = 0.25)) +
  labs(y = "Percent Aeromonas in samples", x = "Sample source") 
dot.plot

# ggsave("dots.pdf", plot = dot.plot, device = "pdf", width = 8, height = 3, units = "in")

