###########################################
### combine & analyze Aeromonas communities
### from 7 datasets
### Lou LaMartina, started Jan 29 2021
###########################################



library(dada2)
library(reshape2)

setwd("~/Desktop/Lab/Projects/Aeromonas/FINAL/Relabun")



########################
### load & prep data ###
########################

# read in files
relabun_files.ls <- lapply(list.files(), function(i) read.csv(i))


# name each dataset
names(relabun_files.ls) <- sapply(strsplit(basename(list.files()), "_"), '[', 1)


# make FASTAs their own column:
# make sample names row names, 
# remove sample name column, 
# transpose -> FASTAs as rows, 
# make new row names (FASTAs) their own column
for(i in 1:length(relabun_files.ls)){
  rownames(relabun_files.ls[[i]]) <- relabun_files.ls[[i]]$Sample_name
  relabun_files.ls[[i]] <- relabun_files.ls[[i]][-1]
  relabun_files.ls[[i]] <- data.frame(t(relabun_files.ls[[i]]))
  relabun_files.ls[[i]] <- data.frame(FASTA = rownames(relabun_files.ls[[i]]), relabun_files.ls[[i]])
}




##################
### merge data ###
##################

# merge all datasets by their FASTAs
merged_files.ls <- list()

for (i in 1:length(relabun_files.ls)) {
  if (i == 1) {
    merged_files.ls[[i]] <- merge(relabun_files.ls[[i]], 
                                  relabun_files.ls[[i + 1]], by = "FASTA", all = TRUE)
  } else if (i > 1 & i < length(relabun_files.ls)) {
    merged_files.ls[[i]] <- merge(merged_files.ls[[i - 1]], 
                                  relabun_files.ls[[i + 1]], by = "FASTA", all = TRUE)
  }
}


# extract the last one - merge() makes you do it one by one
merged_data <- data.frame(merged_files.ls[[length(relabun_files.ls) - 1]])
merged_data[is.na(merged_data)] <- 0




############################
### add dataset variable ###
############################

# transpose -> samples as rows
rownames(merged_data) <- merged_data$FASTA
merged_data <- t(merged_data[-1])


# add sample name variable
merged_data <- data.frame(Sample_name = rownames(merged_data), merged_data)


# extract sample names & store in new list
names.ls <- list()

for(i in 1:length(relabun_files.ls)) {
  names.ls[[i]] <- names(relabun_files.ls[[i]][-1])
}


# add dataset name (line 33)
names(names.ls) <- names(relabun_files.ls)


# turn into data frame
names <- data.frame(Sample_name = unlist(names.ls))
names$Source <- gsub("[[:digit:]]+", "", rownames(names))


# add the data set name to final data by merging
merged_data <- merge(names, merged_data, by = "Sample_name")



























# # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # 
# 
# ### test merging
# 
# # same # of samples?
# ncol(relabun_files.ls[[1]]) +
#   ncol(relabun_files.ls[[2]]) - 1 == ncol(merged_files.ls[[1]])
# # [1] TRUE
# 
# 
# # same # of ASVs?
# length(unique(c(as.character(relabun_files.ls[[1]]$FASTA),
#                 as.character(relabun_files.ls[[2]]$FASTA)))) == 
#   nrow(merged_files.ls[[1]])
# # [1] TRUE
# 
# 
# # does merging work?
# identical(merged_files.ls[[1]],
#           merge(relabun_files.ls[[1]], relabun_files.ls[[2]], by = "FASTA", all = TRUE))
# # [1] TRUE
# 
# 
# # visualize
# gplots::venn(list(file1 = relabun_files.ls[[1]]$FASTA, file2 = relabun_files.ls[[2]]$FASTA))


# ################
# ### checking ###
# ################
# 
# ### merge 1
# 
# ncol(relabun_files.ls[[1]]) + 
#   ncol(relabun_files.ls[[2]]) - 1 # 171
# 
# length(unique(c(as.character(relabun_files.ls[[1]]$FASTA),
#                 as.character(relabun_files.ls[[2]]$FASTA)))) 
# # 11
# 
# test1 <- merge(relabun_files.ls[[1]], relabun_files.ls[[2]], by = "FASTA", all = TRUE)
# dim(test1) # 11 171
# 
# gplots::venn(list(f1 = relabun_files.ls[[1]]$FASTA, f2 = relabun_files.ls[[2]]$FASTA))
# 
# 
# ### merge 2
# 
# ncol(relabun_files.ls[[1]]) + 
#   ncol(relabun_files.ls[[2]]) + 
#   ncol(relabun_files.ls[[3]]) - 2
# # 265
# 
# length(unique(c(as.character(relabun_files.ls[[1]]$FASTA),
#                 as.character(relabun_files.ls[[2]]$FASTA),
#                 as.character(relabun_files.ls[[3]]$FASTA)))) 
# # 34
# 
# test2 <- merge(test1, relabun_files.ls[[3]], by = "FASTA", all = TRUE)
# dim(test2) 
# # 34 265
# 
# gplots::venn(list(f2 = test1$FASTA, f3 = relabun_files.ls[[3]]$FASTA))
# 
# 
# 
# ### merge 3
# 
# ncol(relabun_files.ls[[1]]) + 
#   ncol(relabun_files.ls[[2]]) + 
#   ncol(relabun_files.ls[[3]]) + 
#   ncol(relabun_files.ls[[4]]) - 3
# # 301
# 
# length(unique(c(as.character(relabun_files.ls[[1]]$FASTA),
#                 as.character(relabun_files.ls[[2]]$FASTA),
#                 as.character(relabun_files.ls[[3]]$FASTA),
#                 as.character(relabun_files.ls[[4]]$FASTA)))) 
# # 48
# 
# test3 <- merge(test2, relabun_files.ls[[4]], by = "FASTA", all = TRUE)
# dim(test3) 
# # 48 301
# 
# gplots::venn(list(f2 = test2$FASTA, f3 = relabun_files.ls[[4]]$FASTA))
# 
# 
# 
# 
# ### merge 4
# 
# ncol(relabun_files.ls[[1]]) + 
#   ncol(relabun_files.ls[[2]]) + 
#   ncol(relabun_files.ls[[3]]) + 
#   ncol(relabun_files.ls[[4]]) + 
#   ncol(relabun_files.ls[[5]]) - 4
# # 598
# 
# length(unique(c(as.character(relabun_files.ls[[1]]$FASTA),
#                 as.character(relabun_files.ls[[2]]$FASTA),
#                 as.character(relabun_files.ls[[3]]$FASTA),
#                 as.character(relabun_files.ls[[4]]$FASTA),
#                 as.character(relabun_files.ls[[5]]$FASTA)))) 
# # 95
# 
# test4 <- merge(test3, relabun_files.ls[[5]], by = "FASTA", all = TRUE)
# dim(test4) 
# # [1]  95 598
# 
# gplots::venn(list(f2 = test2$FASTA, f3 = relabun_files.ls[[4]]$FASTA))
# 
# 
# 
# ### final merge
# 
# ncol(relabun_files.ls[[1]]) + 
#   ncol(relabun_files.ls[[2]]) + 
#   ncol(relabun_files.ls[[3]]) + 
#   ncol(relabun_files.ls[[4]]) + 
#   ncol(relabun_files.ls[[5]]) + 
#   ncol(relabun_files.ls[[6]]) + 
#   ncol(relabun_files.ls[[7]]) - 6
# # 803
# 
# length(unique(c(as.character(relabun_files.ls[[1]]$FASTA),
#                 as.character(relabun_files.ls[[2]]$FASTA),
#                 as.character(relabun_files.ls[[3]]$FASTA),
#                 as.character(relabun_files.ls[[4]]$FASTA),
#                 as.character(relabun_files.ls[[5]]$FASTA),
#                 as.character(relabun_files.ls[[6]]$FASTA),
#                 as.character(relabun_files.ls[[7]]$FASTA)))) 
# # 95
# 
# 
# 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 
# 
# aqua <- data.frame(relabun_files.ls[[1]])
# beach <- data.frame(relabun_files.ls[[2]])
# unique(as.character(beach$FASTA), as.character(aqua$FASTA))
# 
# 
# 
# # merge each by FASTA
# length(unique(c(as.character(relabun_files.ls[[1]]$FASTA),
#                 as.character(relabun_files.ls[[2]]$FASTA))))
# 
# gplots::venn(list(beach = relabun_files.ls[[1]]$FASTA, aqua = relabun_files.ls[[2]]$FASTA))
# 
# merge1 <- merge(relabun_files.ls[[1]], relabun_files.ls[[2]], by = "FASTA", all = TRUE)
# 
# 
# 
# length(unique(c(as.character(merge1$FASTA),
#                 as.character(relabun_files.ls[[3]]$FASTA))))
# 
# gplots::venn(list(merge1 = merge1$FASTA, aqua = relabun_files.ls[[3]]$FASTA))
# 
# merge2 <- merge(merge1, relabun_files.ls[[3]], by = "FASTA", all = TRUE)
# 
# 
# 
# 
# length(unique(c(as.character(merge2$FASTA),
#                 as.character(relabun_files.ls[[4]]$FASTA))))
# 
# gplots::venn(list(merge2 = merge2$FASTA, aqua = relabun_files.ls[[4]]$FASTA))
# 
# merge3 <- merge(merge2, relabun_files.ls[[4]], by = "FASTA", all = TRUE)

