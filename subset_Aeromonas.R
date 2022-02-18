#################################################
### Subset Aeromonas genus from total communities
### for Jones et al., 2022
### Lou LaMartina, finalized Jan 3, 2022
#################################################


setwd("~/Desktop/Aeromonas")


#################
### load data ###
#################

# ASV counts
counts_files.ls <- lapply(list.files("./RData/Counts/Total", 
                                     pattern = ".csv", full.names = T), function(i) read.csv(i))

# ASV taxonomy
taxa_files.ls <- lapply(list.files("./RData/Taxonomy", 
                                     pattern = ".csv", full.names = T), function(i) read.csv(i))

# ASV taxonomy
info_files.ls <- lapply(list.files("./RData/Metadata", 
                                   pattern = ".csv", full.names = T), function(i) read.csv(i))


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


# simplify dataset names
Datasets <- sapply(strsplit(names(counts_files.ls), "_"), '[', 1)
names(counts_files.ls) <- Datasets
names(relabun_files.ls) <- Datasets
names(taxa_files.ls) <- Datasets
names(info_files.ls) <- Datasets
Datasets
# [1] "Aquaculture" "Beach"       "Effluent"    "Nase"        "River"       "Salmon"      "Tench"       "Trout"     


# remove "all files" list to save memory
rm(all_files.ls)


# load compiled sample info
info_all <- read.csv("./RData/Aeromonas_16S_metadata.csv")




########################
### subset aeromonas ###
########################

aero_FASTAs.ls <- list()
aero_counts.ls <- list()

for (i in Datasets) {
  aero_FASTAs.ls[[i]] <- subset(taxa_files.ls[[i]], Genus == "Aeromonas")$FASTA
  
  if (length(aero_FASTAs.ls[[i]]) > 1) {
    aero_counts.ls[[i]] <- counts_files.ls[[i]][, colnames(counts_files.ls[[i]]) %in% aero_FASTAs.ls[[i]]]
    
  } else {
    aero_counts.ls[[i]] <- data.frame(counts_files.ls[[i]][, colnames(counts_files.ls[[i]]) %in% aero_FASTAs.ls[[i]]])
    colnames(aero_counts.ls[[i]]) <- aero_FASTAs.ls[[i]]
    rownames(aero_counts.ls[[i]]) <- rownames(counts_files.ls[[i]])

  }
}




##################
### merge data ###
##################

# transpose so FASTA are row names, add source
aero_counts.t.ls <- list()

for (i in Datasets) {
  aero_counts.t.ls[[i]] <- data.frame(t(aero_counts.ls[[i]]))
  aero_counts.t.ls[[i]]$FASTA <- rownames(aero_counts.t.ls[[i]])
}


# merge by FASTAs
merged_counts.ls <- list()

for (i in 1:length(aero_counts.t.ls)) {
  if (i == 1) {
    merged_counts.ls[[i]] <- merge(aero_counts.t.ls[[i]], 
                                   aero_counts.t.ls[[i + 1]], by = "FASTA", all = TRUE)
  } else if (i > 1 & i < length(aero_counts.t.ls)) {
    merged_counts.ls[[i]] <- merge(merged_counts.ls[[i - 1]], 
                                   aero_counts.t.ls[[i + 1]], by = "FASTA", all = TRUE)
  }
}


# extract the last one - merge() makes you do it one by one
merged_counts <- data.frame(merged_counts.ls[[length(aero_counts.t.ls) - 1]])


# change "NA" to zero
merged_counts[is.na(merged_counts)] <- 0


# transpose -> samples as rows
rownames(merged_counts) <- merged_counts$FASTA
merged_counts <- t(merged_counts[-1])


# add sample name variable
merged_counts <- data.frame(Sample_name = rownames(merged_counts), merged_counts)


# extract sample names & store in new list
names.ls <- list()

for(i in 1:length(aero_counts.t.ls)) {
  names.ls[[i]] <- names(aero_counts.t.ls[[i]][-1])
}


# add dataset name
names(names.ls) <- Datasets


# turn into data frame
names <- data.frame(Sample_name = unlist(names.ls))
names$Source <- gsub("[[:digit:]]+", "", rownames(names))


# add the counts set name to final counts by merging
merged_counts <- merge(names, merged_counts, by = "Sample_name")
rownames(merged_counts) <- merged_counts$Sample_name


# remove empty samples / ASVs
dim(merged_counts) # [1] 349 154
counts_filt <- merged_counts[-c(1:2)]
counts_filt <- counts_filt[rowSums(counts_filt) > 0, colSums(counts_filt) > 0]
dim(counts_filt) # [1] 227 150


# add sample info back
counts_filt$Sample_name <- rownames(counts_filt)
counts_filt <- merge(merged_counts[,1:2], counts_filt, by = "Sample_name")
dim(counts_filt) # [1] 227 152
rownames(counts_filt) <- counts_filt$Sample_name


# save!
write.csv(counts_filt, "./RData/Aeromonas_16S_ASVcounts.csv", row.names = FALSE, na = "")

