This GitHub Repository contains code and data included in:

## [Sea ice presence is linked to higher carbon export and vertical microbial connectivity in the Eurasian Arctic Ocean]()
Eduard Fadeev, Andreas Rogge, Simon Ramondenc, Eva-Maria Nöthig, Claudia Wekerle, Christina Bienhold, Ian Salter, Anya M. Waite, Laura Hehemann, Antje Boetius, Morten H. Iversen

### Abstract:
Arctic Ocean sea-ice cover is shrinking due to warming. Long-term sediment trap data shows higher export efficiency of particulate organic carbon in regions with seasonal sea ice compared to regions without sea ice. To investigate this sea-ice enhanced export, we compared how different early summer phytoplankton communities in seasonally ice-free and ice-covered regions of the Fram Strait affect carbon export and vertical dispersal of microbes. In situ collected aggregates revealed two-fold higher carbon export of diatom-rich aggregates in ice-covered regions, compared to Phaeocystis aggregates in the ice-free region. Using microbial source tracking, we found that ice-covered regions were also associated with more surface-born microbial clades exported to the deep-sea. Taken together, our results showed that ice-covered regions are responsible for high export efficiency and provide strong vertical microbial connectivity. Therefore, continuous ice-loss may decrease the vertical export efficiency, and thus the pelagic-benthic coupling, with potential repercussions for Arctic deep-sea ecosystems._


### Content:
**dada2** - output files from the 16S dada2 workflow and metadata. \
**Tables** - metadata. \
**scripts** - supporting scripts for the statistical analysis. \
```Long_term_sed_flux.R``` -  Sediment trap fluxes (Figure 1B).\
```Trajectories_summary_calculations.R``` -  Calculations particle sinking characteristics included in Table 1. \
```Vert_conn_16S_data_preprocess.R``` - Processing of dada2 output into phyloseq object.  \
```Vert_conn_16S_alpha_div.R``` - Mic. comm. alpha diversity calculations (included in supplementary material). \
```Vert_conn_16S_beta_div.R``` -  Mic. comm. beta diversity (Figure 3). \
```Vert_conn_16S_MST.R``` - Microbial source tracking (Figure 4). \
```Vert_conn_16S_enrichment.R``` - Taxonomic enrochment in microbial communities along the water column (Figure 5). \
```Vert_conn_16S_sediment.R``` - Taxonomic composition of shared ASVs between the water column and the deepsea sediment (Figure 6). \

### Author:
Eduard Fadeev([dr.eduard.fadeev@gmail.com](mailto:dr.eduard.fadeev@gmail.com)) 

### Software Versions:
R version 4.1.1 (2021-08-10)\
RStudio version: 1.4.1717
