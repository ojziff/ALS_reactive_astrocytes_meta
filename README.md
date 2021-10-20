# Meta-analysis of human and mouse ALS astrocytes reveals multiomic signatures of inflammatory reactive states
This repository contains scripts to analyse the data and reproduce the figures from the ALS astrocyte meta-analysis and correlation with A1 astrocyte project.

Meta-analysis of human and mouse ALS astrocytes reveals multiomic signatures of inflammatory reactive states.
Oliver J. Ziff*, Benjamin E. Clarke*, Doaa M. Taha*, Hamish Crerar, Nicholas M. Luscombe, Rickie Patani

![GitHub Logo](/working_model.fig5.png)

The scripts are written in Rmarkdown documents for readability and are organised in order of the Figures and Tables in the paper.

All raw and processed mRNA sequencing data generated in this study have been deposited in the NCBI Sequence Read Archive (BioProject Gene Expression Omnibus) under accession number GSE160133. Raw Mass Spectrometry data have been deposited to the ProteomeXchange Consortium via the PRIDE partner repository with the dataset identifier PXD022604. 

Meta-analysis results can be browsed in the interactive web application at https://shiny.crick.ac.uk/ALS_reactive_astrocytes_meta/

In total 16 public sequencing datasets were used in this study:
| Comparison                                     | Sample Type              | Region                                      | Accession #          | Ref                     |
|------------------------------------------------|--------------------------|---------------------------------------------|----------------------|-------------------------|
| Purified CNS cell types from different regions | human biopsy             | Multiple CNS regions                        | GSE73721             | Zhang et al., 2016      |
| Region specific                                | hiPSC                    | Ventral & dorsal; forebrain & spinal cord   | GSE133489            | Bradley et al., 2019    |
| VCP vs Control                                 | hiPSC                    | Ventral spinal cord                         | GSE160133            | Ziff et al., 2021       |
| C9orf72 vs Control                             | hiPSC                    | Ventral spinal cord                         | GSE142730            | Birger et al., 2019     |
| SOD1 vs Control                                | hiPSC                    | Ventral spinal cord                         | GSE102902; GSE99843  | Tyzack et al., 2017     |
| FUS vs Control                                 | hiPSC                    | Ventral spinal cord                         | GSE162892            | Neyrinck et al., 2021   |
| Sporadic ALS vs Control                        | human fibroblast derived | Ventral spinal cord                         | GSE87385             | Ferraiuolo et al., 2016 |
| A1 (treated) vs A0 (untreated)                 | hiPSC                    | Ventral spinal cord                         | synapse ID: 21861229 | Barbar et al., 2020     |
| A1 (treated) vs A0 (untreated)                 | hiPSC                    | Ventral spinal cord (4 different protocols) | GSE182307            | Leng et al., 2020       |
| A1 (treated) vs A0 (untreated)                 | Mouse                    | Cortex                                      | GSE143598            | Guttenplan et al., 2020 |
| MCAO vs A0 (sham) 24 hour samples              | Mouse                    | Cortex                                      | GSE35338             | Zamanian et al., 2012   |
| Spinal cord injury vs uninjured                | Mouse                    | Spinal cord                                 | GSE76097             | Anderson et al., 2017   |
| Astrocyte TDP-43 knockout vs Wildtype          | Mouse                    | Spinal cord                                 | GSE156542            | Peng et al., 2020       |
| Astrocyte Membralin knockout vs Wildtype       | Mouse                    | Motor cortex                                | GSE130763            | Jiang et al., 2019      |
| SOD1 G93A vs Wildtype                          | Mouse                    | Forebrain                                   | GSE143598            | Guttenplan et al., 2020 |
| SOD1 G37R vs Wildtype (TRAP-seq)               | Mouse                    | Spinal cord                                 | GSE74724             | Sun et al., 2015        |

For each RNAseq dataset we process samples with nf-core/rnaseq v3.0 utilising alignment with STAR and read quantification with salmon. Differential gene expression was performed use DESeq2 and mass spectrometry data were analysed with DEP, as per the rmarkdown script. Schematics were created using Biorender.com and merged into figures with Adobe Illustrator.


