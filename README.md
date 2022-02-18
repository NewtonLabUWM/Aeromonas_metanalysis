# Systematic Review and Meta-Analysis of Antimicrobial Susceptibility among  <i>Aeromonas spp.</i> from a One Health perspective

Daniel Christopher Jones, Jenna Lewis,  Emily Lou LaMartina, Andrew James Dahl, Nischala Nagesh Holavanahalli, Ryan Newton, Troy A Skwor



## Microbial community 16S rRNA gene sequence datasets

### Treated wastewater (effluent)

Beattie, R. E., Skwor, T., & Hristova, K. R. (2020). Survivor microbial populations in post-chlorinated wastewater are strongly associated with untreated hospital sewage and include ceftazidime and meropenem resistant populations. [Science of The Total Environment](https://www.sciencedirect.com/science/article/pii/S0048969720337074?casa_token=iv1F7xNgfiAAAAAA:gU5u_5YeFvghDMprAboJJppcGjHLi0bVumTizm2T97Y8S42JHilexx9VlJ6_P27r4RPf_nbPoQ), 740, 140186.


### Aquaculture

Bartelme, R. P., McLellan, S. L., & Newton, R. J. (2017). Freshwater recirculating aquaculture system operations drive biofilter bacterial community shifts around a stable nitrifying consortium of ammonia-oxidizing archaea and comammox Nitrospira. [Frontiers in microbiology](https://www.frontiersin.org/articles/10.3389/fmicb.2017.00101/full), 8, 101.


### Rivers

McClary-Gutierrez, J. S., Driscoll, Z., Nenn, C., & Newton, R. J. (2021). Human Fecal Contamination Corresponds to Changes in the Freshwater Bacterial Communities of a Large River Basin. [Microbiology spectrum](https://journals.asm.org/doi/full/10.1128/Spectrum.01200-21), 9(2), e01200-21.



### Seafood

1.  <b>Nase.  </b> Guivier, E., Pech, N., Chappaz, R., & Gilles, A. (2020). Microbiota associated with the skin, gills, and gut of the fish Parachondrostoma toxostoma from the Rhône basin. [Freshwater Biology](https://onlinelibrary.wiley.com/doi/abs/10.1111/fwb.13437), 65(3), 446-459.

2.  <b>Salmon.  </b> Webster, T. M. U., Consuegra, S., Hitchings, M., & de Leaniz, C. G. (2018). Interpopulation variation in the Atlantic salmon microbiome reflects environmental and genetic diversity. [Applied and environmental microbiology](https://aem.asm.org/content/84/16/e00691-18), 84(16).

3.  <b>Tench.  </b> Dulski, T., Kozłowski, K., & Ciesielski, S. (2020). Habitat and seasonality shape the structure of tench (Tinca tinca L.) gut microbiome. [Scientific reports](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7064478/), 10(1), 1-11.

4.  <b>Trout.  </b> Vasemägi, A., Visse, M., & Kisand, V. (2017). Effect of environmental factors and an emerging parasitic disease on gut microbiome of wild salmonid fish. [MSphere](https://msphere.asm.org/content/2/6/e00418-17), 2(6).



### Beaches


Cloutier, D. D., Alm, E. W., & McLellan, S. L. (2015). Influence of land use, nutrients, and geography on microbial communities and fecal indicator abundance at Lake Michigan beaches. [Applied and environmental microbiology](https://aem.asm.org/content/81/15/4904.short), 81(15), 4904-4913.


## Estimate proportions of <i>Aeromonas</i> in each dataset

### 1.  Obtain raw FASTQs
- <b>NCBI Short Read Archive.  </b> Find published studies with 16S rRNA gene amplicon sequencing done on seafood-associated samples (see ./SRA_code).
- <b>UWM archives.  </b> Compile aquatic-environmental datasets from previously published works.
- see Aeromonas_16S_metadata.csv


### 2.  Process raw 16S rRNA gene reads
- <b>DADA2.  </b> Quality filter, merge paired-end reads, remove chimeric sequences, assign taxonomy with DADA2 (see ./DADA2_code).


### 3.  Parse sequences classified as genus <i>Aeromonas</i>
- see Aeromonas_metanalyis.R


![image](https://github.com/NewtonLabUWM/Aeromonas_metanalysis/blob/main/DADA2_code/Fig1.png)
