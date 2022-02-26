#################################################
### Subset Aeromonas genus from total communities
### for Jones et al., 2022
### Lou LaMartina, finalized Feb 26, 2022
#################################################


setwd("~/Desktop/Lab/Projects/Others/Aeromonas/Final")


library(ggplot2)
library(RColorBrewer)


#################
### load data ###
#################

# ASV counts
counts_files.ls <- lapply(list.files("./Counts/", pattern = ".csv"), function(i) read.csv(i))

# ASV taxonomy
taxa_files.ls <- lapply(list.files("./Taxa/", pattern = ".csv"), function(i) read.csv(i))


# simplify dataset names
Datasets <- names(counts_files.ls)
Datasets
# [1] "Aquaculture" "Beach"       "Effluent"    "Nase"        "River"       "Salmon"      "Tench"       "Trout" 


# change row names to sample names
for (i in Datasets) {
  rownames(counts_files.ls[[i]]) <- counts_files.ls[[i]]$Sample_name
  counts_files.ls[[i]] <- counts_files.ls[[i]][-1]
}


# convert to relative abundance
relabun_files.ls <- list()
for(i in 1:length(counts_files.ls)){
  relabun_files.ls[[i]] <- counts_files.ls[[i]] / rowSums(counts_files.ls[[i]])
}


# change row names to FASTA sequences
for (i in Datasets) {
  rownames(taxa_files.ls[[i]]) <- taxa_files.ls[[i]]$FASTA
}


# load all info
info <- read.csv("AeromonasMetanalysis_16S_metadata.csv")




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


# transpose so FASTA are row names, add source
aero_relabun.t.ls <- list()

for (i in 1:length(aero_relabun.ls)) {
  aero_relabun.t.ls[[i]] <- data.frame(t(aero_relabun.ls[[i]]))
  aero_relabun.t.ls[[i]]$FASTA <- rownames(aero_relabun.t.ls[[i]])
}


# merge by FASTAs
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
dim(merged_relabun) # [1] 333 151




###################
### proportions ###
###################

# calculate proportion aeromonas in each sample
propAero <- merged_relabun
propAero$totalAero <- rowSums(propAero)
propAero <- data.frame(File = rownames(propAero), totalAero = propAero$totalAero)


# add sample info
propAero <- merge(propAero, info, by = "File")


# calculate mean & std deviation
stats <- propAero[c(2,7,8)]
stats <- data.frame(mean = aggregate(totalAero ~  SampleSource + SampleType, mean, data = stats),
                    sd = aggregate(totalAero ~ SampleSource + SampleType, sd, data = stats))[-c(4,5)]
#   mean.SampleSource  mean.SampleType mean.totalAero sd.totalAero
# 1                Gut  Atlantic Salmon   0.0515949753 0.1643431706
# 2               Skin  Atlantic Salmon   0.0067877637 0.0146675815
# 3              Beach       Beach sand   0.0001207746 0.0003177870
# 4                Gut      Brown trout   0.0026353595 0.0062711770
# 5        Aquaculture           Filter   0.0002607541 0.0004491563
# 6              Beach       Lake water   0.0008704838 0.0013614855
# 7                Gut             Nase   0.0491976537 0.0449155558
# 8               Skin             Nase   0.0248607623 0.0241214670
# 9         Wastewater Post-chlorinated   0.0011812171 0.0007419777
# 10        Wastewater  Pre-chlorinated   0.0030126622 0.0036916766
# 11             River      River water   0.0028700847 0.0031514172
# 12       Aquaculture             Tank   0.0034283053 0.0089282294
# 13               Gut            Tench   0.2040527504 0.2695263323



# most in each dataset
for(i in Datasets){
  top <- names(sort(rowSums(merged_relabun[rownames(merged_relabun) %in% rownames(counts_files.ls[[i]]),]), decreasing = T)[1])
  cat("\n", i)
  print(info[info$File == top,])
  print(rowSums(merged_relabun[rownames(merged_relabun) %in% top,]))
}
# Aquaculture      File  BioProject    BioSample Country  Location SampleSource SampleType
# 27 RA0110A PRJNA491504 SAMN10097279     USA Milwaukee  Aquaculture       Tank
# RA0110A 
# 0.02711231 
# 
# Beach                        File  BioProject    BioSample Country   Location SampleSource SampleType
# 146 MI_NB1_09_23_13_20140219 PRJNA272391 SAMN03282217     USA Ferrysburg        Beach Lake water
# MI_NB1_09_23_13_20140219 
# 0.005676409 
# 
# Effluent           File  BioProject    BioSample Country  Location SampleSource      SampleType
# 272 SRR11487887 PRJNA622864 SAMN14531099     USA Milwaukee   Wastewater Pre-chlorinated
# SRR11487887 
# 0.009420787 
# 
# Nase         File BioProject    BioSample Country Location SampleSource SampleType
# 71 ERR3458748 PRJEB33555 SAMEA5816803  France  St Just          Gut       Nase
# ERR3458748 
# 0.1372662 
# 
# River       File  BioProject    BioSample Country      Location SampleSource  SampleType
# 242 NL00296 PRJNA665728 SAMN16265636     USA Nichols Creek        River River water
# NL00296 
# 0.01253992 
# 
# Salmon         File BioProject      BioSample        Country Location SampleSource      SampleType
# 32 ERR2129783 PRJEB22688 SAMEA104307779 United Kingdom    Frome          Gut Atlantic Salmon
# ERR2129783 
# 0.5467512 
# 
# Tench          File  BioProject    BioSample Country Location SampleSource SampleType
# 322 SRR9039948 PRJNA542255 SAMN11618244  Poland   Lezyny          Gut      Tench
# SRR9039948 
# 0.8916264 
# 
# Trout          File  BioProject    BioSample Country Location SampleSource  SampleType
# 299 SRR6354471 PRJNA388139 SAMN08136043 Estonia    Voesu          Gut Brown trout
# SRR6354471 
# 0.02722898 


# ordering
propAero$combo <- paste(propAero$SampleSource, propAero$SampleType)
labs  <- sort(unique(propAero$combo))
names(labs) <- c(18, 17, 21, 20, 12, 14, 13, 11, 19, 15, 16, 22, 23)
labs <- data.frame(labs = paste(names(labs), labs), combo = labs)
propAero <- merge(propAero, labs, by = "combo")

facs <- sort(unique(propAero$SampleSource))
names(facs) <- c(3, 5, 1, 4, 2, 6)
facs <- data.frame(facs = paste(names(facs), facs), SampleSource = facs)
propAero <- merge(propAero, facs, by = "SampleSource")

colors <- c("#08519C", "#4EB3D3", "#A50F15", "#41AB5D", "#6A51A3", "#F16913")


# plot
dot.plot <-
  ggplot(propAero, aes(x = labs, y = totalAero)) +
  geom_hline(yintercept = c(0.01, 0.1), size = 0.25, color = "grey80", linetype = "dotted") +
  geom_violin(size = 0.5, alpha = 0.5, show.legend = F, 
              aes(fill = SampleSource), color = NA) +
  geom_jitter(width = 0.1, size = 1, show.legend = F,
              aes(fill = SampleSource, color = SampleSource)) +
  facet_grid(. ~ facs, scales = "free", space = "free") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  theme_classic() +
  scale_y_continuous(trans = scales::pseudo_log_trans(0.001, 10),
                     breaks = c(0, 0.01, 0.1, 1), 
                     labels = c("0%", "1%", "10%", "100%")) +
  theme(axis.text.x = element_text(size = 6, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.y = element_text(size = 10, color = "black", face = "bold"),
        axis.title.x = element_blank(),
        strip.text = element_blank(),
        strip.background = element_rect(colour = NA, fill = NA),
        panel.border = element_rect(size = .75, color = "grey80", fill = NA),
        axis.ticks = element_line(size = 0.25),
        axis.line = element_line(size = 0.25)) +
  labs(y = "Percent Aeromonas in samples", x = "Sample source") 
dot.plot

#ggsave("dots.pdf", plot = dot.plot, device = "pdf", width = 10, height = 4, units = "in")


