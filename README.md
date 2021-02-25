# ALS_A1_astrocytes
This repository contains scripts to analyse the data and reproduce the figures from the ALS astrocyte meta-analysis and correlation with A1 astrocyte project.

ALS astrocytes share multiomic signatures and functional deficits with A1 reactive astrocytes. Benjamin E. Clarke*, Oliver J. Ziff*, Doaa M. Taha*, Hamish Crerar, Nicholas M. Luscombe, Rickie Patani

The scripts are written in Rmarkdown documents for readability and are organised in order of the Figures and Tables in the paper.

All raw and processed mRNA sequencing data generated in this study have been deposited in the NCBI Sequence Read Archive (BioProject Gene Expression Omnibus) under accession number GSE160133. Raw Mass Spectrometry data have been deposited to the ProteomeXchange Consortium (http://proteomecentral.proteomexchange.org) via the PRIDE partner repository with the dataset identifier PXD022604. 

For analysis of human primary central nervous system cell types we downloaded RNA sequencing from GSE73721. For comparison with regionally defined astrocytes we downloaded RNAseq from GSE133489. Cytokine-stimulated iPSC derived astrocyte RNA sequencing is available at syn21861181. iPSC derived astrocytes carrying ALS mutations are available at GSE142730 (C9orf72), GSE102902 and GSE99843 (SOD1 mutants and control respectively). Mouse SOD1 astrocyte TRAP-seq is available at GSE74724. Astrocyte-specific TARDBP and membralin knockout mouse spinal cord RNA-seq is available at GSE156542 and GSE130763 respectively.

For each dataset we process samples with nf-core/rnaseq v3.0 utilising alignment with STAR and read quantification with salmon. Differential gene expression was performed use DESeq2 and mass spectrometry data were analysed with DEP, as per rmarkdown script. Schematics were created using Biorender.com and merged into figures with Adobe Illustrator.


