# Systematic Review and Meta-Analysis of Antimicrobial Susceptibility among  <i>Aeromonas spp.</i> from a One Health perspective

Daniel Christopher Jones, Jenna Lewis,  Emily Lou LaMartina, Andrew James Dahl, Nischala Nagesh Holavanahalli, Ryan Newton, Troy A Skwor



## Microbial community 16S rRNA gene sequence datasets

![image](https://github.com/loulanomics/Aeromonas_metanalysis/blob/main/SRA_code/Aeromonas_dataset_overview.png)


## Estimate proportions of <i>Aeromonas</i> in each dataset

### 1.  Obtain raw FASTQs
- <b>NCBI Short Read Archive.  </b> Find published studies with 16S rRNA gene amplicon sequencing done on seafood-associated samples (see ./SRA_code).
- <b>UWM archives.  </b> Compile aquatic-environmental datasets from previously published works.
- see Aeromonas_16S_metadata.csv


### 2.  Process raw 16S rRNA gene reads
- <b>DADA2.  </b> Quality filter, merge paired-end reads, remove chimeric sequences, assign [taxonomy](https://www.arb-silva.de/documentation/release-138/) with [DADA2](https://www.bioconductor.org/packages/release/bioc/html/dada2.html) (see ./DADA2_code).


### 3.  Parse sequences classified as genus <i>Aeromonas</i>
- see Aeromonas_metanalyis.R


![image](https://github.com/NewtonLabUWM/Aeromonas_metanalysis/blob/main/DADA2_code/Fig1.png)



