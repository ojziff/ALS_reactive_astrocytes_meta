# ALS_A1_astrocytes
This repository contains scripts to analyse the data and reproduce the figures from the ALS and A1 astrocyte comparison paper

ALS astrocytes share transcriptomic signatures and functional deficits with A1 reactive astrocytes Benjamin E. Clarke*, Oliver J. Ziff1*, Doaa M. Taha*, Hamish Crerar, Hannah Franklin, Nicholas M. Luscombe, Rickie Patani

The scripts are written in Rmarkdown documents for readability and are organised in order of the Figures in the paper.

All RNA sequencing data generated for this study is deposited at NCBI GEO under accession number GSE160133. RAW Mass Spectrometry data have been deposited to the ProteomeXchange Consortium via the PRIDE partner repository with the dataset identifier PXD022604.

Previously published iPSC derived astrocytes carrying ALS mutations are available at GSE142730 (C9orf72), GSE102902 and GSE99843 (SOD1 mutants and control respectively). Cytokine-stimulated iPSC derived astrocytes are available at syn21861181. Mouse SOD1 astrocyte TRAP-seq is available at GSE74724. Astrocyte-specific TARDBP and membralin knockout mouse spinal cord RNA-seq is available at GSE156542 and GSE130763 respectively.

For each dataset we process samples with nf-core/rnaseq v3.0 utilising alignment with STAR and read quantification with salmon. Differential gene expression was performed use DESeq2 as per rmarkdown script. Mass spectrometry data were analysed with DEP. Schematics were created using Biorender.com and merged into figures with Adobe Illustrator.


