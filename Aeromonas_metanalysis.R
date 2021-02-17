##################################
### Aeromonas meta-analysis
### Lou LaMartina, started 2/15/21
##################################

# all data trimmed with cutadapt, processed with dada2, assigned with silva v138

# data not shown (for future studies): 
# perch guts, WWTP influent time series, WWTP influent across US, lake/harbor

# aeromonas not found:
# clams, HMP, cow

# studies:
# bartelme 2017 = aquaculture
# cloutier 2015 = beach
# beattie 2020 = wwtpeff
# mcclary IP = river
# guivier 2020 = nase
# webster 2018 = salmon
# dulski 2020 = tench
# vasemägi 2017 = trout


setwd("~/Desktop/Lab/Projects/Aeromonas/FINAL")

library(readr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)


#######################
### getting started ###
#######################

# read in files - 7.033162 mins
# all_files.ls <- lapply(list.files(pattern = ".csv"), function(i) read_csv(i))


# add file names
# names(all_files.ls) <- sapply(strsplit(list.files(pattern = ".csv"), "\\."), '[', 1)


# ~ or ~ load all files list - 3.971166 secs (RECOMMENDED)
all_files.ls <- readRDS("./RData/All_files_list.RData")


# subset ASV counts, relative abundances, taxonomy, and metadata files
counts_files.ls <- all_files.ls[grep("counts", names(all_files.ls))]
relabun_files.ls <- all_files.ls[grep("relabun", names(all_files.ls))]
taxa_files.ls <- all_files.ls[grep("taxa", names(all_files.ls))]
info_files.ls <- all_files.ls[grep("info", names(all_files.ls))]


# counts -
# change row names to sample names,
# remove sample name column
for (i in names(counts_files.ls)) {
  rownames(counts_files.ls[[i]]) <- counts_files.ls[[i]]$Sample_name
  counts_files.ls[[i]] <- counts_files.ls[[i]][-1]
}


# relatives -
# change row names to sample names,
# remove sample name column
for (i in names(relabun_files.ls)) {
  rownames(relabun_files.ls[[i]]) <- relabun_files.ls[[i]]$Sample_name
  relabun_files.ls[[i]] <- relabun_files.ls[[i]][-1]
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


# remove "all files" list to save memory
rm(all_files.ls)




########################
### subset aeromonas ###
########################

aero_FASTAs.ls <- list()
aero_counts.ls <- list()
aero_relabun.ls <- list()

for (i in Datasets) {
  aero_FASTAs.ls[[i]] <- subset(taxa_files.ls[[i]], Genus == "Aeromonas")$FASTA
  
  if (length(aero_FASTAs.ls[[i]]) > 1) {
    aero_counts.ls[[i]] <- counts_files.ls[[i]][, colnames(counts_files.ls[[i]]) %in% aero_FASTAs.ls[[i]]]
    aero_relabun.ls[[i]] <- relabun_files.ls[[i]][, colnames(relabun_files.ls[[i]]) %in% aero_FASTAs.ls[[i]]]
  
  } else {
    aero_counts.ls[[i]] <- data.frame(counts_files.ls[[i]][, colnames(counts_files.ls[[i]]) %in% aero_FASTAs.ls[[i]]])
    colnames(aero_counts.ls[[i]]) <- aero_FASTAs.ls[[i]]
    rownames(aero_counts.ls[[i]]) <- rownames(counts_files.ls[[i]])
    
    aero_relabun.ls[[i]] <- data.frame(relabun_files.ls[[i]][, colnames(relabun_files.ls[[i]]) %in% aero_FASTAs.ls[[i]]])
    colnames(aero_relabun.ls[[i]]) <- aero_FASTAs.ls[[i]]
    rownames(aero_relabun.ls[[i]]) <- rownames(relabun_files.ls[[i]])
  }
}




##################
### merge data ###
##################

# transpose so FASTA are row names, add source
aero_counts.t.ls <- list()
aero_relabun.t.ls <- list()

for (i in Datasets) {
  aero_counts.t.ls[[i]] <- data.frame(t(aero_counts.ls[[i]]))
  aero_counts.t.ls[[i]]$FASTA <- rownames(aero_counts.t.ls[[i]])

  aero_relabun.t.ls[[i]] <- data.frame(t(aero_relabun.ls[[i]]))
  aero_relabun.t.ls[[i]]$FASTA <- rownames(aero_relabun.t.ls[[i]])
}



###################
### merge by FASTAs

### counts ###
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


### relatives ###
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
merged_counts <- data.frame(merged_counts.ls[[length(aero_counts.t.ls) - 1]])
merged_relabun <- data.frame(merged_relabun.ls[[length(aero_relabun.t.ls) - 1]])


# change "NA" to zero
merged_counts[is.na(merged_counts)] <- 0
merged_relabun[is.na(merged_relabun)] <- 0



########################
### add dataset variable


### counts ###

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



### relatives ###
 
# transpose -> samples as rows
rownames(merged_relabun) <- merged_relabun$FASTA
merged_relabun <- t(merged_relabun[-1])


# add sample name variable
merged_relabun <- data.frame(Sample_name = rownames(merged_relabun), merged_relabun)


# extract sample names & store in new list
names.ls <- list()

for(i in 1:length(aero_relabun.t.ls)) {
  names.ls[[i]] <- names(aero_relabun.t.ls[[i]][-1])
}


# add dataset name
names(names.ls) <- Datasets


# turn into data frame
names <- data.frame(Sample_name = unlist(names.ls))
names$Source <- gsub("[[:digit:]]+", "", rownames(names))


# add the relabun set name to final relabun by merging
merged_relabun <- merge(names, merged_relabun, by = "Sample_name")
rownames(merged_relabun) <- merged_relabun$Sample_name


# remove empty samples / ASVs
dim(merged_relabun) # [1] 349 154
relabun_filt <- merged_relabun[-c(1:2)]
relabun_filt <- relabun_filt[rowSums(relabun_filt) > 0, colSums(relabun_filt) > 0]
dim(relabun_filt) # [1] 227 150


# add sample info back
relabun_filt$Sample_name <- rownames(relabun_filt)
relabun_filt <- merge(merged_relabun[,1:2], relabun_filt, by = "Sample_name")
dim(relabun_filt) # [1] 227 152
rownames(relabun_filt) <- relabun_filt$Sample_name




###################
### proportions ###
###################

# calculate proportion aeromonas in each sample
propAero <- relabun_filt
propAero$totalAero <- rowSums(relabun_filt[-c(1:2)])
propAero <- propAero[c(1, 2, ncol(propAero))]
propAero[ncol(propAero)]


# add biome variable
propAero$Biome <- "Environment"
propAero$Biome[propAero$Source == "Tench" |
                 propAero$Source == "Trout" |
                 propAero$Source == "Salmon" |
                 propAero$Source == "Nase"] <- "Seafood"
propAero$Source[propAero$Source == "Effluent"]  <- "Treated\nwastewater"


# plot
dot.plot <-
  ggplot(propAero, aes(x = Source, y = log10(totalAero))) +
  geom_boxplot(width = 0.3, size = 0.3, outlier.color = NA, alpha = 0.5) +
  geom_jitter(width = 0.2, shape = 1, size = 1, alpha = 0.5) +
  facet_grid(. ~ Biome, scales = "free_x") +
  theme_classic() +
  scale_y_continuous(breaks = seq(0, -5), labels = c("100%", "10%", "1%", "0.1%", "0.01%", "0.001%")) +
  theme(axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 8, color = "black", face = "bold"),
        axis.title.x = element_text(size = 8, color = "black", face = "bold"),
        panel.border = element_rect(color = "grey80", fill = NA, size = 0.5), 
        strip.text = element_text(size = 8, color = "black", face = "bold"),
        strip.background = element_rect(colour = "grey80", fill = "grey80"),
        axis.ticks = element_line(size = 0.25),
        axis.line = element_line(size = 0.25)) +
  labs(y = "Percent Aeromonas in samples")
dot.plot

# ggsave("./Plots/dots.pdf", plot = dot.plot, device = "pdf", width = 6, height = 4, units = "in")




############
### save ###
############

# counts
for(i in Datasets){
  write.csv(data.frame(Sample_name = rownames(aero_counts.ls[[i]]), aero_counts.ls[[i]]), 
            paste0("./Data/Counts/Aeromonas/", i,"_Aero_counts.csv"),
            row.names = FALSE)
}


# relatives
for(i in Datasets){
  write.csv(data.frame(Sample_name = rownames(aero_relabun.ls[[i]]), aero_relabun.ls[[i]]), 
            paste0("./Data/Relatives/Aeromonas/", i,"_Aero_relabun.csv"),
            row.names = FALSE)
}








# # # # # # # # # # # # # # # # # # # # # #
# save.image("Aeromonas_environment.RData")
# # # # # # # # # # # # # # # # # # # # # #




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# colnames(taxa_files.ls[["Trout_taxa"]])[1] <- "FASTA"
# colnames(taxa_files.ls[["Tench_taxa"]])[1] <- "FASTA"
# colnames(taxa_files.ls[["Salmon_taxa"]])[1] <- "FASTA"
# 
# colnames(counts_files.ls[["Trout_counts"]])[1] <- "Sample_name"
# colnames(counts_files.ls[["Tench_counts"]])[1] <- "Sample_name"
# colnames(counts_files.ls[["Salmon_counts"]])[1] <- "Sample_name"
# 
# colnames(relabun_files.ls[["Trout_relabun"]])[1] <- "Sample_name"
# colnames(relabun_files.ls[["Tench_relabun"]])[1] <- "Sample_name"
# colnames(relabun_files.ls[["Salmon_relabun"]])[1] <- "Sample_name"
# 
# 
# write.csv(counts_files.ls[["Tench_counts"]], "Tench_counts.csv", row.names = FALSE)
# write.csv(taxa_files.ls[["Tench_taxa"]], "Tench_taxa.csv", row.names = FALSE, na = "")
# write.csv(relabun_files.ls[["Tench_relabun"]], "Tench_relabun.csv", row.names = FALSE)
# 
# write.csv(counts_files.ls[["Trout_counts"]], "Trout_counts.csv", row.names = FALSE)
# write.csv(taxa_files.ls[["Trout_taxa"]], "Trout_taxa.csv", row.names = FALSE, na = "")
# write.csv(relabun_files.ls[["Trout_relabun"]], "Trout_relabun.csv", row.names = FALSE)
# 
# write.csv(counts_files.ls[["Salmon_counts"]], "Salmon_counts.csv", row.names = FALSE)
# write.csv(taxa_files.ls[["Salmon_taxa"]], "Salmon_taxa.csv", row.names = FALSE, na = "")
# write.csv(relabun_files.ls[["Salmon_relabun"]], "Salmon_relabun.csv", row.names = FALSE)




# # transpose so FASTA are row names
# relabun_files.t.ls <- list()
# for (i in names(relabun_files.ls)) {
#   relabun_files.t.ls[[i]] <- data.frame(t(relabun_files.ls[[i]]))
# }


# counts_files.t.ls <- list()
# for (i in names(counts_files.ls)) {
#   counts_files.t.ls[[i]] <- data.frame(t(counts_files.ls[[i]]))
# }



# Hospital_info <- info_files.ls[["Hospital_info"]]
# Hospital_info <- Hospital_info[grep("chlor", Hospital_info$`Library Name`, ignore.case = TRUE),]
# 
# Hospital_counts <- counts_files.ls[["Hospital_counts"]]
# Hospital_counts <- subset(Hospital_counts, rownames(Hospital_counts) %in% Hospital_info$Run)
# Hospital_counts <- Hospital_counts[, colSums(Hospital_counts) > 0]
# 
# Hospital_relabun <- relabun_files.ls[["Hospital_relabun"]]
# Hospital_relabun <- subset(Hospital_relabun, rownames(Hospital_relabun) %in% Hospital_info$Run)
# Hospital_relabun <- Hospital_relabun[, colSums(Hospital_relabun) > 0]

# Hospital_taxa <- taxa_files.ls[["Hospital_taxa"]]
# Hospital_taxa <- subset(Hospital_taxa, rownames(Hospital_taxa) %in% colnames(Hospital_counts))

# info_files.ls[["Hospital_info"]] <- Hospital_info
# counts_files.ls[["Hospital_counts"]] <- Hospital_counts
# relabun_files.ls[["Hospital_relabun"]] <- Hospital_relabun
# taxa_files.ls[["Hospital_taxa"]] <- Hospital_taxa

# Hospital_info <- info_files.ls[["Hospital_info"]]
# Hospital_counts <- counts_files.ls[["Hospital_counts"]]
# Hospital_relabun <- relabun_files.ls[["Hospital_relabun"]]
# 
# 
# write.csv(counts_files.ls[["Hospital_counts"]], "Hospital_counts.csv", row.names = TRUE)
# write.csv(taxa_files.ls[["Hospital_taxa"]], "Hospital_taxa.csv", row.names = TRUE, na = "")
# write.csv(relabun_files.ls[["Hospital_relabun"]], "Hospital_relabun.csv", row.names = TRUE)
# write.csv(info_files.ls[["Hospital_info"]], "Hospital_info.csv", row.names = TRUE)


# for (i in names(all_files.ls)) {
#   print(colnames(all_files.ls[[i]])[1])
# }
# 
# temp1 <- all_files.ls[["Trout_relabun"]]
# colnames(all_files.ls[["Trout_relabun"]])[1] <- "Sample_name"
# 
# temp2 <- all_files.ls[["Trout_taxa"]]
# colnames(all_files.ls[["Trout_taxa"]])[1] <- "FASTA"

# # convert from tibbles to data frames
# for (i in names(all_files.ls)) {
#   all_files.ls[[i]] <- data.frame(all_files.ls[[i]])
# }

# names(all_files.ls)[9] <- "Effluent_counts"
# names(all_files.ls)[10] <- "Effluent_info"
# names(all_files.ls)[11] <- "Effluent_relabun"
# names(all_files.ls)[12] <- "Effluent_taxa"



# ncol(aero_counts.t.ls[[1]]) - 1 +
#   ncol(aero_counts.t.ls[[2]]) - 1 +
#   ncol(aero_counts.t.ls[[3]]) - 1 +
#   ncol(aero_counts.t.ls[[4]]) - 1 +
#   ncol(aero_counts.t.ls[[5]]) - 1 +
#   ncol(aero_counts.t.ls[[6]]) - 1 +
#   ncol(aero_counts.t.ls[[7]]) - 1 +
#   ncol(aero_counts.t.ls[[8]]) - 1
# # [1] 357
# 
# 
# # FASTA
# length(unique(c(rownames(aero_counts.t.ls[[1]]),
#                 rownames(aero_counts.t.ls[[2]]),
#                 rownames(aero_counts.t.ls[[3]]),
#                 rownames(aero_counts.t.ls[[4]]),
#                 rownames(aero_counts.t.ls[[5]]),
#                 rownames(aero_counts.t.ls[[6]]),
#                 rownames(aero_counts.t.ls[[7]]),
#                 rownames(aero_counts.t.ls[[8]]))))
# # [1] 152
# 
# 
# 
# 
# 
# ### test merging
# 
# # same # of samples?
# ncol(aero_counts.t.ls[[1]]) +
#   ncol(aero_counts.t.ls[[2]]) - 1 == ncol(test)
# # [1] TRUE
# 
# 
# # same # of ASVs?
# length(unique(c(as.character(aero_counts.t.ls[[1]]$FASTA),
#                 as.character(aero_counts.t.ls[[2]]$FASTA)))) == 
#   nrow(test)
# # [1] TRUE
# 
# 
# # does merging work?
# identical(test,
#           merge(aero_counts.t.ls[[1]], aero_counts.t.ls[[2]], by = "FASTA", all = TRUE))
# # [1] TRUE
# 
# 
# # visualize
# gplots::venn(list(file1 = aero_counts.t.ls[[1]]$FASTA, file2 = aero_counts.t.ls[[2]]$FASTA))



####################
### stacked bars ###
####################

# 
# 
# 
# topFASTAs <- data.frame(totalRel = colSums(relabun_filt[-c(1:2)]),
#                         FASTA = colnames(relabun_filt[-c(1:2)]))
# topFASTAs <- topFASTAs[order(topFASTAs$totalRel, decreasing = TRUE),]
# topFASTAs$ASV <- paste("Aero", sprintf("%0.3d", 1:nrow(topFASTAs)))
# topFASTAs$topASV <- topFASTAs$ASV
# topFASTAs$topASV[21:nrow(topFASTAs)] <- "Other"
# 
# 
# relabun_filt.m <- merge(relabun_filt.m, topFASTAs, by = "FASTA")
# 
# colors <- c(colorRampPalette(brewer.pal(9, "Spectral"))(length(unique(topFASTAs$topASV)) - 1), "grey80")
# 
# 
# bars <- 
#   ggplot(relabun_filt.m, aes(x = Source, y = Relabun, color = topASV, fill = topASV)) +
#   geom_bar(stat = "identity", position = "fill", show.legend = TRUE) +
#   scale_fill_manual(values = colors) +
#   scale_color_manual(values = colors) +
#   theme_classic() +
#   theme(legend.text = element_text(size = 8)) +
#   guides(color = guide_legend(keyheight = 1, keywidth = 1, units = "in", ncol = 1),
#          fill = guide_legend(keyheight = 1, keywidth = 1, units = "in", ncol = 1))
# 
# ggsave("bars.pdf", plot = bars, device = "pdf", width = 10, height = 8, units = "in")



