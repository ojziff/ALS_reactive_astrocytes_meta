.libPaths("/camp/lab/luscomben/home/users/ziffo/.conda/envs/r4.0.3/lib/R/library")
library(tidyverse) # loads ggplot2, tibble, readr, tidyr, purr, dplyr, stringr, forcats
library(usethis)
library(devtools)
library(DESeq2)
library(BiocParallel)
library("grid")
library("ggplotify")
library("RColorBrewer")
library(scales)
library(ggpubr)
library(ggsci)
library(readxl)
library(vroom)
library(rhdf5)
library("tximeta")
library("tximport")
library(rtracklayer)
library("rmarkdown")
library("data.table")
library("gsubfn")
library("pheatmap")
library(ComplexHeatmap)
library(circlize)
library('ggrepel')
library("clusterProfiler")
library("org.Hs.eg.db")
library("limma") 
library("AnnotationDbi")
library("Glimma")
library("colorRamps")
library("GO.db")
library("fgsea")
library("geneplotter")
library("genefilter") 
library("sva")
library("vsn")
library("RUVSeq")
library("gprofiler2")
library("topGO") 
library(Biobase)
library(SummarizedExperiment)
library(VennDiagram)
library(patchwork)
library(rstatix)
library(GenomicFeatures)
library(Gviz)
library(trackViewer)
library(GeneOverlap)
library(ggdendro)
library(dendextend)
library(biovizBase)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Rnorvegicus.UCSC.rn6)
library(GenomicScores) # GC content
library(phastCons100way.UCSC.hg38)
library(DEP)
library("EnhancedVolcano")
# library("systemPipeR")
library(ggforce)
library(pcr)
library(gplots)
library(plyranges)
library(seqinr)
library(GenomicRanges)
library(GeneStructureTools) # devtools::install_github("betsig/GeneStructureTools")
library(notNMD) # devtools::install_github('betsig/notNMD')
library(gdata)
library(fgsea)
library(piano)
library(parallel)
library(GSEABase)
library(snowfall)
library("GencoDymo") # remotes::install_github("monahton/GencoDymo")
library(VarCon) # remotes::install_github("caggtaagtat/VarCon")  https://bioconductor.org/packages/devel/bioc/manuals/VarCon/man/VarCon.pdf
library(crayon)
library(ggvenn)
# library("ggVennDiagram")
library(progeny)
library(dorothea)
library(biomaRt)
library(irlib) # install_github("jsha129/irlib")
library(tidygraph)
library(ggraph)
library(conflicted) # for package conflicts
select <- dplyr::select # let dplyr::select win package conflict
filter <- dplyr::filter
n <- dplyr::n
list <- base::list
rename <- dplyr::rename
unique <- base::unique
options(stringsAsFactors = F)

# GO terms ----------------------

gprofiler_database <- gmtPathways("/camp/home/ziffo/home/genomes/gene-ontology/gprofiler_full_hsapiens.name.edit.gmt") # IDs removed to see term names
# https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp GO pathways
curated.gene.sets.msigdb <- gmtPathways("/camp/home/ziffo/home/genomes/gene-ontology/c2.all.v7.1.symbols.gmt") # includes KEGG, REACTOME, PID, BIOCARTA, CGP
kegg.pathways.msigdb <- gmtPathways("/camp/home/ziffo/home/genomes/gene-ontology/c2.cp.kegg.v7.1.symbols.gmt")
reactome.pathways.msigdb <- gmtPathways("/camp/home/ziffo/home/genomes/gene-ontology/c2.cp.reactome.v7.1.symbols.gmt")
go.pathways.msigdb <- gmtPathways("/camp/home/ziffo/home/genomes/gene-ontology/c5.all.v7.1.symbols.gmt")
go.bp.pathways.msigdb<- gmtPathways("/camp/home/ziffo/home/genomes/gene-ontology/c5.bp.v7.1.symbols.gmt")
go.mf.pathways.msigdb<- gmtPathways("/camp/home/ziffo/home/genomes/gene-ontology/c5.mf.v7.1.symbols.gmt")
go.cc.pathways.msigdb <- gmtPathways("/camp/home/ziffo/home/genomes/gene-ontology/c5.cc.v7.1.symbols.gmt")
immune.pathways.msigdb <- gmtPathways("/camp/home/ziffo/home/genomes/gene-ontology/c7.all.v7.1.symbols.gmt")
RNA_BINDING_PROTEINS <- read_excel("/camp/home/ziffo/home/genomes/gene-ontology/rbp_consensus_gerstberger_2014.xls", sheet = "RBP table") %>% dplyr::pull(`gene name`)
nucleocytoplasmic_genes_of_interest <- c(RNA_BINDING_PROTEINS, 
                                         go.pathways.msigdb$GO_NUCLEAR_PORE, go.pathways.msigdb$GO_TRANSCRIPTION_EXPORT_COMPLEX, go.pathways.msigdb$GO_NUCLEAR_TRANSCRIBED_MRNA_CATABOLIC_PROCESS_NONSENSE_MEDIATED_DECAY, go.pathways.msigdb$GO_RNA_EXPORT_FROM_NUCLEUS, go.pathways.msigdb$GO_NUCLEOCYTOPLASMIC_CARRIER_ACTIVITY, go.pathways.msigdb$GO_NEGATIVE_REGULATION_OF_NUCLEOCYTOPLASMIC_TRANSPORT)
ALS_RBPs <- c("TARDBP", "SFPQ", "FUS", "ELAVL3", "HNRNPA1", "HNRNPA2B1", "MATR3", "ATXN2", "TAF15", "TIA1", "EWSR1")
glia.genes <- c("S100B", "SOX9")
astrocyte.genes <- c("GFAP","AQP4","PLA2G7","SLC39A12","MLC1","DIO2","SLC14A1","ALDH1L1","ALDOC","TTPA","ACSBG1","CHRDL1","SLC4A4","SLC1A2","SLC25A18","SLC1A3","F3","PPP1R3G","FZD2","MERTK","EZR","EVA1A","GJB6","HAPLN1","RFX4","PAPSS2","SLC15A2","PPP1R3C","TLR3","ACOT11","ATP1A2","BMPR1B","PRODH","GLI3","TMEM47","SLC9A3R1","CTH","NTSR2","SLC7A10","VCAM1","FGFR3","CCDC80","ENTPD2","CYBRD1","KCNE5","FAM20A","TNC","TLCD1","S1PR1","CBS","PBXIP1","GRIN2C","ADHFE1","AGT","GLDC","SLC7A2","GJA1","PDK4","EGFR","SOX9","CLDN10","PLCD4","ID4","FMO1","EMP2","LONRF3","HTRA1","MGST1","THRSP")
panreactive.astrocyte.genes <- c("LCN2", "STEAP4", "S1PR3", "TIMP1", "HSPB1", "CXCL10", "CD44", "OSMR", "CP", "VIM", "GFAP",
                                 "SERPINA3", "ASPG") # "SERPINA5",
A1.astrocyte.genes <- c("C3", "HLA-A", "HLA-B", "HLA-C","HLA-E","HLA-F", "MICA", "H2-T23","H2-D1", #, "HCP5", "HLA-H", "HLA-G", "HLA-K", "HLA-L", "AL645929.2",
                        "GBP2", "AMIGO2", "SERPING1","GGTA1P","GGTA1", #"GLT6D1", "A3GALT2", "IRGC",
                        "FBLN5", "UGT1A1", "FKBP5", "PSMB8", "SRGN","IIGP1")#, "MX1")
A2.astrocyte.genes <- c("S100A10", "EMP1", "CLCF1", "TGM1", "PTX3", "SPHK1", "CD109", "PTGS2", "SLC10A6", "TM4SF1", "B3GNT5", "CD14")
astrocyte.subtype.genes <- list(panreactive.astrocyte.genes,A1.astrocyte.genes,A2.astrocyte.genes)
astrocyte.reactive.genes <- c(panreactive.astrocyte.genes, A1.astrocyte.genes, A2.astrocyte.genes) %>% unique
astrocyte.reactivity.markers = data.frame("gene_name" = astrocyte.reactive.genes) %>% 
  mutate(group = case_when(gene_name %in% A1.astrocyte.genes ~ "A1",
                           gene_name %in% A2.astrocyte.genes ~ "A2",
                           gene_name %in% panreactive.astrocyte.genes ~ "Pan-reactive"))
# astrocyte.reactivity.markers.mouse = data.frame("gene_name.mouse" = astrocyte.reactive.genes) %>% 
#   mutate(group = case_when(gene_name %in% A1.astrocyte.genes ~ "A1",
#                            gene_name %in% A2.astrocyte.genes ~ "A2",
#                            gene_name %in% panreactive.astrocyte.genes ~ "Pan-reactive"))

cytoskeleton = c("GFAP", "NES", "SYNM", "VIM")
metabolism = c("ALDOC", "FABP7", "MAOB", "TSPO")
chaperone = c("CRYAB", "HSPB1")
secreted = c("C3", "CHI3L1", "LCN2", "SERPINA3", "THBS1", "THBS2", "MT2A", "MT1E") #, "MT1F",  "MT1X", "MT1G", "MT1M", "MT3", "MT1H", "MT1A")
signalling = c("NFATC4","NFAT5", "NTRK2", "IL17R", "S100B", "SOX9", "STAT3", "IRAS1", "IRAS2") # "NFATC1", "NFATC2","NFATC3",
channels = c("SLC1A3", "SLC1A2", "KCNJ10")
escartin.markers = tibble(gene_name = c(cytoskeleton, metabolism, chaperone, secreted, signalling, channels), group = case_when(gene_name %in% cytoskeleton ~ "cytoskeleton", gene_name %in% metabolism ~ "metabolism", gene_name %in% chaperone ~ "chaperone", gene_name %in% secreted ~ "secreted", gene_name %in% signalling ~ "signalling", gene_name %in% channels ~ "channels"))

all.astrocyte.markers <- c("GRN", "PSEN1", "LRP1", "APP", "VIM", "IFNGR1", "GRN", "LAMC3", "PSEN1", "MT3", "TREM2", "CDK6", "IFNG", "LRP1", "PLP1", "IL1B", "NR1D1", "ADORA2A", "SMO", "LDLR", "GFAP", "KRAS", "TSPAN2", "IL6", "TLR4", "APP", "S100A8", "EIF2B5", "TTBK1", "EGFR", "ROR2", "FPR2", "LAMB2", "C1QA", "BACE2", "POU3F2", "DRD1", "ROR1", "MAPT", "NF1", "C5AR1", "DLL1", "AGER", "MIR181C", "MIR181B2", "MIR181B1", "TNF", "CNTF", "MIR142", "KRAS", "NF1", "CCR2", "GPR183", "WNT1", "FZD1", "CTNNB1", "HEXB", "CCL2", "CCR2", "APCDD1", "MMP14", "GPR183", "SCRIB", "CCL3", "IFNGR1", "GRN", "PSEN1", "TREM2", "IFNG", "LRP1", "IL1B", "NR1D1", "ADORA2A", "SMO", "LDLR", "IL6", "APP", "TTBK1", "EGFR", "FPR2", "C1QA", "BACE2", "MAPT", "C5AR1", "AGER", "MIR181C", "MIR181B2", "MIR181B1", "TNF", "CNTF", "MIR142", "SOX8", "PAX6", "VIM", "IFNGR1", "GRN", "LAMC3", "PSEN1", "MT3", "TREM2", "ABL1", "MAPK1", "MAPK3", "MAG", "CDK6", "PRPF19", "SOX6", "IFNG", "NR2E1", "HES1", "ID2", "EPHA4", "LRP1", "PLP1", "TTC21B", "SOX9", "IL1B", "NKX2-2", "BMP2", "NR1D1", "ADORA2A", "LIF", "SMO", "LDLR", "GFAP", "KRAS", "TSPAN2", "IL6ST", "SERPINE2", "IL6", "BIN1", "TLR4", "GCM1", "NTRK3", "MBD1", "APP", "S100A8", "EIF2B5", "TTBK1", "EGFR", "NOTCH1", "S100B", "TAL1", "SHH", "STAT3", "MAP2K1", "ROR2", "GPR37L1", "FPR2", "LAMB2", "ID4", "C1QA", "DAB1", "CLCF1", "PTPN11", "F2", "BACE2", "NOG", "CNTN2", "POU3F2", "DRD1", "ROR1", "MAPT", "NF1", "C5AR1", "HES5", "DLL1", "AGER", "MIR181C", "MIR181B2", "MIR181B1", "TNF", "CNTF", "MIR142", "MAG", "PRPF19", "NR2E1", "HES1", "ID2", "EPHA4", "BMP2", "NR1D1", "LIF", "LDLR", "IL6ST", "SERPINE2", "IL6", "BIN1", "NTRK3", "MBD1", "TTBK1", "NOTCH1", "GPR37L1", "ID4", "DAB1", "CLCF1", "F2", "NOG", "CNTN2", "NF1", "HES5", "MIR181C", "MIR181B2", "MIR181B1", "MIR142", "MAG", "PRPF19", "HES1", "ID2", "BMP2", "LIF", "IL6ST", "SERPINE2", "BIN1", "TTBK1", "NOTCH1", "CLCF1", "MIR142", "NR2E1", "NR1D1", "LDLR", "NTRK3", "MBD1", "GPR37L1", "ID4", "DAB1", "F2", "NOG", "NF1", "HES5", "MIR181C", "MIR181B2", "MIR181B1", "SOX8", "SOX9", "GCM1", "TAL1", "NR1D1", "LDLR", "IL6", "TTBK1", "MIR181C", "MIR181B2", "MIR181B1", "MIR142", "NR1D1", "LDLR", "MIR181C", "MIR181B2", "MIR181B1", "TTBK1", "MIR142", "CCR2", "GPR183", "KCNK2", "MT3", "EZR", "MLC1", "SLC1A2", "ATP1B2", "GFAP", "SYT4", "EIF2S1", "APP", "SLC7A11", "PINK1", "GRM2", "GJB2", "AQP4", "KCNJ10", "SLC17A8", "GRM3", "DMD", "ADGRG1", "MT3", "MLC1", "ATP1B2", "GFAP", "EIF2S1", "AQP4", "SLC17A8", "ADGRG1", "MLH1", "MSH2", "TSC2", "IFNG", "MSH3", "MSH6", "PMS2", "POT1", "APC", "IDH1", "BRCA2", "TP53", "ERBB2", "CDKN2A", "AIFM1", "TSC1", "IDH2", "NF2", "NF1")
all.astrocyte.markers <- c(all.astrocyte.markers,astrocyte.reactive.genes)
microglia.genes <- c("ITGAM", "CX3CR1", "CCL3", "CSFIR", "CCL4", "P2RY12", "C1QB", "PLEK", "GPR183")
myeloid.markers = c("TMEM119", "AIF1", "C1QA", "CSF1R", "ITGAM", "PTPRC", "CEMIP2", "RUNX1", "CD83", "TNF", "CCL3", "CCL2", "IL1A", "TLR2") #  "CX3CR1",
all.microglia.go <- c(go.pathways.msigdb$GO_MICROGLIA_DIFFERENTIATION, go.pathways.msigdb$GO_MICROGLIAL_CELL_PROLIFERATION, go.pathways.msigdb$GO_MICROGLIAL_CELL_MIGRATION, microglia.genes, myeloid.markers)
# GO:0002265	astrocyte activation involved in immune response	GRN	PSEN1	LRP1	APP
# GO:0014002	astrocyte development	VIM	IFNGR1	GRN	LAMC3	PSEN1	MT3	TREM2	CDK6	IFNG	LRP1	PLP1	IL1B	NR1D1	ADORA2A	SMO	LDLR	GFAP	KRAS	TSPAN2	IL6	TLR4	APP	S100A8	EIF2B5	TTBK1	EGFR	ROR2	FPR2	LAMB2	C1QA	BACE2	POU3F2	DRD1	ROR1	MAPT	NF1	C5AR1	DLL1	AGER	MIR181C	MIR181B2	MIR181B1	TNF	CNTF	MIR142
# GO:0021896	forebrain astrocyte differentiation	KRAS	NF1
# GO:0021897	forebrain astrocyte development	KRAS	NF1
# GO:0035700	astrocyte chemotaxis	CCR2	GPR183
# GO:0036520	astrocyte-dopaminergic neuron signaling	WNT1	FZD1	CTNNB1
# GO:0043615	astrocyte cell migration	HEXB	CCL2	CCR2	APCDD1	MMP14	GPR183	SCRIB	CCL3
# GO:0048143	astrocyte activation	IFNGR1	GRN	PSEN1	TREM2	IFNG	LRP1	IL1B	NR1D1	ADORA2A	SMO	LDLR	IL6	APP	TTBK1	EGFR	FPR2	C1QA	BACE2	MAPT	C5AR1	AGER	MIR181C	MIR181B2	MIR181B1	TNF	CNTF	MIR142
# GO:0048708	astrocyte differentiation	SOX8	PAX6	VIM	IFNGR1	GRN	LAMC3	PSEN1	MT3	TREM2	ABL1	MAPK1	MAPK3	MAG	CDK6	PRPF19	SOX6	IFNG	NR2E1	HES1	ID2	EPHA4	LRP1	PLP1	TTC21B	SOX9	IL1B	NKX2-2	BMP2	NR1D1	ADORA2A	LIF	SMO	LDLR	GFAP	KRAS	TSPAN2	IL6ST	SERPINE2	IL6	BIN1	TLR4	GCM1	NTRK3	MBD1	APP	S100A8	EIF2B5	TTBK1	EGFR	NOTCH1	S100B	TAL1	SHH	STAT3	MAP2K1	ROR2	GPR37L1	FPR2	LAMB2	ID4	C1QA	DAB1	CLCF1	PTPN11	F2	BACE2	NOG	CNTN2	POU3F2	DRD1	ROR1	MAPT	NF1	C5AR1	HES5	DLL1	AGER	MIR181C	MIR181B2	MIR181B1	TNF	CNTF	MIR142
# GO:0048710	regulation of astrocyte differentiation	MAG	PRPF19	NR2E1	HES1	ID2	EPHA4	BMP2	NR1D1	LIF	LDLR	IL6ST	SERPINE2	IL6	BIN1	NTRK3	MBD1	TTBK1	NOTCH1	GPR37L1	ID4	DAB1	CLCF1	F2	NOG	CNTN2	NF1	HES5	MIR181C	MIR181B2	MIR181B1	MIR142
# GO:0048711	positive regulation of astrocyte differentiation	MAG	PRPF19	HES1	ID2	BMP2	LIF	IL6ST	SERPINE2	BIN1	TTBK1	NOTCH1	CLCF1	MIR142
# GO:0048712	negative regulation of astrocyte differentiation	NR2E1	NR1D1	LDLR	NTRK3	MBD1	GPR37L1	ID4	DAB1	F2	NOG	NF1	HES5	MIR181C	MIR181B2	MIR181B1
# GO:0060018	astrocyte fate commitment	SOX8	SOX9	GCM1	TAL1
# GO:0061888	regulation of astrocyte activation	NR1D1	LDLR	IL6	TTBK1	MIR181C	MIR181B2	MIR181B1	MIR142
# GO:0061889	negative regulation of astrocyte activation	NR1D1	LDLR	MIR181C	MIR181B2	MIR181B1
# GO:0061890	positive regulation of astrocyte activation	TTBK1	MIR142
# GO:2000458	regulation of astrocyte chemotaxis	CCR2	GPR183
# GO:2000464	positive regulation of astrocyte chemotaxis	CCR2
# GO:0097449	astrocyte projection	KCNK2	MT3	EZR	MLC1	SLC1A2	ATP1B2	GFAP	SYT4	EIF2S1	APP	SLC7A11	PINK1	GRM2	GJB2	AQP4	KCNJ10	SLC17A8	GRM3	DMD	ADGRG1
# GO:0097450	astrocyte end-foot	MT3	MLC1	ATP1B2	GFAP	EIF2S1	AQP4	SLC17A8	ADGRG1
# HP:0100707	Abnormality of the astrocytes	MLH1	MSH2	TSC2	IFNG	MSH3	MSH6	PMS2	POT1	APC	IDH1	BRCA2	TP53	ERBB2	CDKN2A	AIFM1	TSC1	IDH2	NF2	NF1

focal.adhesion.go <- c("LAP3", "CD99", "LASP1", "ABCB4", "ITGA3", "ITGA2B", "RALA", "CD9", "MRC2", "TSPAN9", "PLAUR", "EHD3", "CAPN1", "FHL1", "VIM", "CD44", "ARHGAP31", "VCL", "TNC", "CTNNA1", "HSPA5", "LIMA1", "BCAR1", "CYBA", "SYNE2", "GDI2", "PPP1R12A", "NCKAP1", "RPL18", "CNN2", "SLC9A3R2", "TLE2", "RHOA", "FGFR3", "PABPC1", "CDC42", "MAP4K4", "RPL31", "ACTN1", "LIMS2", "PVR", "FERMT2", "CLASP1", "HACD3", "ACTB", "REXO2", "MCAM", "USP33", "APBB1IP", "ACTN2", "ITGA8", "FAP", "TNS1", "SENP1", "DNM2", "EPB41L2", "RAB21", "PTPRC", "ITGB5", "RPS5", "FAT1", "RAB10", "CD59", "CPNE3", "CTTN", "NOX4", "TRIP6", "ADD1", "CASS4", "ASAP3", "RPL6", "RPLP0", "PXN", "SLC9A1", "ICAM1", "ITGA6", "SNAP23", "EZR", "SORBS1", "JAK2", "NRP1", "MISP", "MAPK1", "TRIOBP", "PACSIN2", "RPL3", "MYH9", "ZFYVE21", "PROCR", "FERMT1", "HCK", "MAPRE1", "CD99L2", "ARHGEF7", "FLT1", "MAPK3", "CORO2B", "RPS16", "RPS19", "ITGB8", "CAV2", "CAV1", "HSPB1", "LIMK1", "ENG", "PDLIM1", "GIT1", "RPL19", "PFN1", "SLC6A4", "YWHAE", "LAMTOR3", "HSPA8", "LPXN", "CBL", "CD81", "RPS13", "PPFIBP1", "CORO1C", "TNS2", "TRPV4", "ARPC3", "NEDD9", "OPRM1", "WASF1", "PTK7", "HSPA9", "PDGFRB", "PRKAR2A", "ACTR3", "EPB41L5", "ITGB6", "ITGA4", "RPS15", "IL1RL1", "FHL2", "RND3", "RPL22", "RHOU", "ARHGEF2", "DOCK7", "CD46", "CNN3", "GNA13", "TEK", "SORBS3", "PTK2B", "CAT", "RPL5", "PLAU", "ADGRE5", "LRP1", "SDC4", "STX16", "RPS10", "AHNAK", "EFNB2", "RPL23", "VASP", "FLRT3", "RRAS", "AP006333.1", "AIF1L", "MAP2K2", "PTPN12", "CDC42EP1", "RAC2", "FLNC", "ARHGAP22", "PALLD", "AJUBA", "STARD8", "CNN1", "NECTIN2", "ACTN4", "ARPC1B", "PAK4", "AKAP12", "CAP1", "RPL27", "PPFIA1", "TNS4", "ITGB4", "FLOT2", "PTPRA", "DCTN4", "MPRIP", "KRAS", "RRAS2", "NUMB", "YWHAQ", "ANXA1", "TES", "AVIL", "FLNB", "LMO7", "LCP1", "TNS3", "RAC1", "ARPC5L", "TLN1", "HMGA1", "FLOT1", "MDC1", "ATAT1", "SDCBP", "RDX", "KIF23", "ITGA11", "RPLP1", "ADAM10", "BCAR3", "ACTR2", "ITGAV", "ARHGAP24", "SCARB2", "PARVG", "GIT2", "ITGB7", "IQGAP1", "TGFB1I1", "ARMC5", "CDH13", "RPS2", "PDPK1", "CLTC", "GRB7", "RPS11", "RPL13A", "EPHA2", "HSPG2", "RPS8", "MTF2", "DCAF6", "PRUNE1", "PIP5K1A", "S100A7", "ARF1", "RHOB", "FBLN7", "LIMD1", "PHLDB2", "LPP", "RPS3A", "ARHGAP26", "G3BP1", "GNA12", "EGFR", "SH3KBP1", "CASK", "MSN", "ZNF185", "RPL7", "GSN", "RPL7A", "RSU1", "CAPN5", "PAK1", "RPS3", "HYOU1", "ITGB1", "DIXDC1", "TWF1", "ADAM17", "DST", "ARL14EP", "TADA1", "GJA1", "DAB2", "THY1", "PGM5", "ENAH", "SORBS2", "RPL30", "MMP14", "FZD1", "SHROOM4", "CSRP1", "ACTC1", "ZYX", "ITGB2", "RPL8", "NPHS1", "ITGA5", "JAK1", "FBLIM1", "NEXN", "ARPC5", "DDR2", "NCSTN", "CAPN2", "XIRP2", "ARPC2", "NFASC", "CLASP2", "RPL9", "ANXA5", "ITGA2", "RPS14", "DLC1", "SLC4A2", "YWHAZ", "HNRNPK", "ARF6", "ILK", "HSP90B1", "B2M", "PPIB", "YWHAB", "MAPRE2", "PDIA3", "TPM4", "SRP68", "CTNNB1", "FAM107A", "IRF2", "XIRP1", "ADAM9", "SNTB2", "TM4SF20", "MAP2K1", "PTK2", "LIMS1", "ALCAM", "YWHAG", "PDCD6IP", "CDH2", "RPS9", "RPS7", "TLN2", "KLF11", "SNTB1", "BSG", "GNB2", "SYNPO2", "CORO1B", "CFL1", "RPL38", "DAG1", "PEAK1", "CSPG4", "JUP", "RPL4", "CSRP2", "YES1", "RHOG", "RPLP2", "CD151", "FLII", "PLEC", "GAK", "CALR", "FZD2", "NPM1", "ADGRB1", "FES", "CAV3", "RPS17", "FHL3", "ACTG1", "UBOX5", "FLRT2", "LMLN", "P4HB", "ATP6V0C", "PIP5K1C", "PPP1CC", "CHP1", "SPRY4", "NHS", "PEAK3", "FOCAD", "PARVB", "PPIA", "EVL", "AFAP1", "MME", "PDLIM7", "FLNA", "ANXA6", "IGF2R", "PCBP2", "SRC", "SVIL", "DPP4", "PARVA", "RPL37A", "RPL12", "MPZL1", "RPS4X", "ITGBL1", "RPL10A", "L1CAM", "TGM2", "LAYN", "HSPA1B", "HSPA1A", "ARL2", "PPP1CB", "RPS29", "ITGA1", "TSPAN4", "RPS18", "ALKBH6", "PI4KA", "SCARF2", "ACTN3", "LIMS4", "LIMS3", "ITGB3", "AC068234.1", "CYFIP1", "PRAG1", "MARCKS")
collagen.formation.go <- c("FOXC1","COL11A1","COL5A3","ADAMTS2","TGFB2","MMP11","CHADL","COMP","PLOD3","AEBP1","TGFBR1","COL1A1","COL12A1","LOX","LOXL3","FMOD","P4HA1","COLGALT1","PXDN","COL5A1","LOXL2","CYP1B1","EMILIN1","LOXL4","ADAMTS14","COL2A1","LUM","RB1","P3H4","DPT","SFRP2","SERPINH1","VIPAS39","ADAMTS3","ACAN","DDR2","COL1A2","ATP7A","GREM1","SERPINF2","TNXB","COL3A1","FOXC2","ANXA2","COL14A1","OPTC","NF1","COL11A2","COL5A2","SCX","MIR29B1","MIR29B2 MMP25","MRC2","FAP","COL19A1","ADAMTS2","MMP2","MMP11","MMP9","CST3","MMP15","VSIR","CTSD","MMP8","MMP19","PEPD","MMP24","CTSL","MMP7","MMP20","MMP27","MMP13","ADAMTS14","FURIN","CTSK","ADAM15","MMP3","ITGB1","MMP21","MMP16","ADAMTS3","MMP14","CTSS","CTSB","MMP10","MMP26","KLK6","PHYKPL","TMPRSS6","MMP23B","PRTN3","MMP1","COL13A1","MMP17","COL15A1","MMP12","MMP28","PRSS2","MMP25","MRC2","VIM","TRAM2","FAP","COL19A1","ADAMTS2","MMP2","P2RX7","P3H2","SUCO","MMP11","HIF1A","MMP9","CST3","RGCC","MMP15","SMPD3","TGFB1","PLOD3","ENG","VSIR","COL1A1","P3H3","TNS2","PPARD","PDGFRB","ERRFI1","RAP1A","P3H1","CTSD","MMP8","MYB","ARG1","CCN2","TGFB3","GOT1","MMP19","PEPD","AMELX","BMP4","MMP24","ID1","COL5A1","PPARG","CTSL","IL6","MMP7","MMP20","MMP27","MMP13","EMILIN1","ADAMTS14","FURIN","ARRB2","CBX8","P3H4","RCN3","CTSK","ADAM15","SERPINH1","MMP3","ITGB1","VIPAS39","MMP21","MMP16","ADAMTS3","MMP14","CREB3L1","CYGB","WNT4","CTSS","NPPC","IHH","UCN","ITGA2","COL1A2","CTSB","LARP6","SERPINB7","MFAP4","MMP10","MMP26","SERPINF2","KLK6","TNXB","PHYKPL","CIITA","F2","F2R","VPS33B","TMPRSS6","MMP23B","PRTN3","HDAC2","MMP1","COL13A1","MMP17","COL15A1","MIR149","MIR218-1","MIR218-2","SCX","MMP12","MMP28","PRSS2","MIR145","MIR92A1","MIR29B1","MIR29A","MIR21","MIR29B2","MIR92A2","COL4A4","COL1A1","COL4A2","ITGA11","UBASH3B","DDR2","ITGA2","SYK","COL4A3","OSCAR","COL4A1","COL4A5","COL4A6","DDR1","TLL1","COL9A2","COL23A1","LAMA3","LAMC2","COL11A1","COL17A1","P4HA2","COL5A3","COL4A4","COL19A1","PLOD1","COL16A1","ADAMTS2","P3H2","ITGA6","COL9A3","TLL2","MMP9","COL20A1","PCOLCE","PLOD3","COL1A1","P3H3","COL12A1","COL9A1","LOX","COL7A1","LOXL3","P3H1","P4HA1","COL10A1","COL21A1","LOXL1","COLGALT1","PXDN","COL5A1","ITGB4","LOXL2","COL4A2","CTSL","CTSV","MMP7","MMP20","MMP13","LOXL4","ADAMTS14","COL2A1","COL6A1","COL6A2","COL8A1","SERPINH1","P4HA3","MMP3","PLOD2","ADAMTS3","COL26A1","CTSS","COL6A3","PCOLCE2","COL1A2","CTSB","PPIB","BMP1","COL3A1","COL4A3","COL22A1","CRTAP","COL24A1","COL8A2","COL6A5","CD151","PLEC","COL18A1","P4HB","COL4A1","COL14A1","COL4A5","COL25A1","COL27A1","LAMB3","COL13A1","COL4A6","COLGALT2","COL11A2","COL5A2","COL15A1","COL6A6","COL28A1")
extracellular.matrix.go <- c("CFLAR","ST7","ITGAL","ITGA3","ITGA2B","NOX1","MMP25","DCN","CAPN1","PHLDB1","CD44","TNFRSF1B","IBSP","TIMP2","TLL1","VCAN","CDH1","TNC","GPM6B","COL9A2","ELN","LAMC3","COL23A1","LAMA3","FOXC1","LAMC2","COL11A1","COL17A1","TNFRSF1A","BCL3","CLASP1","NTN4","FSCN1","ICAM3","NFKB2","FBLN1","ITGA8","FAP","COL5A3","CPB2","COL4A4","COL19A1","ITGB5","ITGAE","COL16A1","B4GALT1","ERO1B","ADAMTS2","MMP2","NID2","KIF9","ICAM1","LAMB4","LAMB1","ITGA6","CCDC80","CMA1","COL9A3","TGFB2","TMEM38B","TLL2","ABL1","MADCAM1","MMP11","PDGFB","CHADL","CTSG","MMP9","FERMT1","CST3","LAMA1","TIMP1","RGCC","MMP15","HAS3","SMPD3","CRISPLD2","FOXF1","ERCC2","TGFB1","ICAM4","ICAM5","HAS1","COMP","HPN","ITGB8","CAV2","CAV1","DNAJB6","SERPINE1","PLOD3","AEBP1","MEGF9","TGFBR1","ECM2","ENG","SPOCK2","KAZALD1","SH3PXD2A","ICAM2","COL1A1","VTN","VWF","MYF5","COL12A1","ADTRP","COL9A1","NR2E1","SMOC2","LAMA4","LOX","SPARC","COL7A1","MPV17","ITGB6","ITGA4","LOXL3","FN1","TNR","QSOX1","EXOC8","NID1","MFAP2","MMP8","TTR","CCN2","SPP1","TGFBI","FMOD","PLG","SERAC1","RECK","P4HA1","PRDX4","MMP19","LRP1","COL10A1","MATN4","SOX9","TCF15","MMP24","CAPNS1","KDR","LOXL1","NCAN","COLGALT1","PXDN","COL5A1","LAMA5","POMT1","GFAP","RAMP2","MATN3","ITGB4","BCAN","HAPLN2","POSTN","MYH11","SPINK5","LOXL2","PDGFRA","COL4A2","ETS1","CTSL","ADAM19","ITGA7","AGT","LAMC1","LCP1","IL6","CTSV","FOXF2","FLOT1","SULF1","MMP7","MMP20","MMP27","MMP13","THBS1","ITGA11","ADAM10","CYP1B1","EMILIN1","LOXL4","ADAMTS14","ITGAV","FGF2","SLC39A8","FBN2","COL2A1","LUM","ITGB7","RB1","FBLN5","FURIN","ITGAX","GFOD2","P3H4","COL6A1","COL6A2","APP","HSPG2","CCN1","ITGA10","DPT","ADAMTSL4","CTSK","ADAM15","ITGA9","COL8A1","PHLDB2","SFRP2","HAPLN1","SCUBE3","CSGALNACT1","NOTCH1","ADAM12","GAS2","HSD17B12","SERPINH1","MMP3","ITGB1","VIPAS39","ADAM8","DSPP","DMP1","ABI3BP","WNT3A","MMP21","JAM2","ADAMTS5","MMP16","ADAMTS3","ITGAD","MMP14","MYO1E","CREB3L1","ACAN","F11R","ADAMTS4","SCUBE1","CARMIL2","ITGB2","MPZL3","FGFR4","NPHS1","ITGA5","PDPN","MATN1","VCAM1","DDR2","OLFML2B","CAPN2","CTSS","COL6A3","ELF3","IHH","FBLN2","CLASP2","ADAMTS9","PTX3","APBB2","MELTF","ITGA2","EGFLAM","KLKB1","COL1A2","TNFRSF11B","ATP7A","HTRA1","JAM3","SPINT1","FBN1","WDR72","MFAP4","MMP10","GREM1","SMAD3","MMP26","SPINT2","SERPINF2","KLK4","KLK2","KLK5","TNXB","BMP1","COL3A1","NPNT","CTRB1","CTRB2","COL4A3","KLK7","PTK2","COL22A1","ANTXR1","ITGAM","FSHR","HAS2","COL24A1","RXFP1","FGG","FGA","FGB","COL8A2","ANGPTL7","LAMB2","TPSAB1","BSG","EFEMP2","HPSE2","ADAMTS20","NDNF","DAG1","SH3PXD2B","A2M","FOXC2","RIC8A","VWA1","WASHC1","OTOL1","BGN","ANXA2","COL18A1","GAS6","WT1","FLRT2","OLFML2A","TMPRSS6","COL4A1","THSD4","COL14A1","COL4A5","AGRN","OPTC","MMP23B","NOXO1","SULF2","LAMA2","MMP1","NF1","COL27A1","CD47","LAMB3","PDGFA","COL13A1","ELANE","COL4A6","MFAP5","DPP4","ADAMTSL2","ERO1A","MMP17","EGFL6","COL11A2","COL5A2","COL15A1","DDR1","PRSS1","VIT","SERPINB5","COLQ","ITGA1","COL28A1","ATXN1L","TNF","MARCOL","CAPNS2","ITGB3","SCX","PECAM1","MMP12","MMP28","MIR98","PRSS2","MIR29B1","MIR29B2")
glutamate.uptake.go <- unique(c(go.pathways.msigdb$GO_GLUTAMATE_BINDING, go.pathways.msigdb$GO_GLUTAMATE_BIOSYNTHETIC_PROCESS, go.pathways.msigdb$GO_GLUTAMATE_GATED_CALCIUM_ION_CHANNEL_ACTIVITY, go.pathways.msigdb$GO_GLUTAMATE_METABOLIC_PROCESS, gprofiler_database$`glutamatergic synapse`, gprofiler_database$`glutamate reuptake`))
phagocytosis.go <- unique(c(go.pathways.msigdb$GO_REGULATION_OF_PHAGOCYTOSIS, go.pathways.msigdb$GO_PHAGOCYTOSIS))
all.reactivity.go <- c(focal.adhesion.go, collagen.formation.go, extracellular.matrix.go, go.pathways.msigdb$GO_RESPONSE_TO_WOUNDING,
                       RNA_BINDING_PROTEINS,
                       go.pathways.msigdb$GO_IMMUNE_EFFECTOR_PROCESS, go.pathways.msigdb$GO_CYTOKINE_MEDIATED_SIGNALING_PATHWAY, go.pathways.msigdb$GO_RESPONSE_TO_CYTOKINE, 
                       go.pathways.msigdb$GO_RESPONSE_TO_INTERFERON_ALPHA, go.pathways.msigdb$GO_RESPONSE_TO_INTERFERON_BETA, go.pathways.msigdb$GO_RESPONSE_TO_INTERFERON_GAMMA,
                       gprofiler_database$`Cellular responses to stress`, all.astrocyte.markers, 
                       go.pathways.msigdb$GO_PHAGOCYTOSIS, glutamate.uptake.go, go.pathways.msigdb$GO_REGULATION_OF_PHAGOCYTOSIS,
                       go.pathways.msigdb$GO_EXCITATORY_SYNAPSE,go.pathways.msigdb$GO_NEUROFILAMENT, go.pathways.msigdb$GO_NEUROFILAMENT_CYTOSKELETON_ORGANIZATION,
                       go.pathways.msigdb$GO_NEUROINFLAMMATORY_RESPONSE, go.pathways.msigdb$GO_GLIAL_CELL_ACTIVATION, 
                       go.pathways.msigdb$GO_NEURON_APOPTOTIC_PROCESS, go.pathways.msigdb$GO_NEURON_CELL_CELL_ADHESION, go.pathways.msigdb$GO_NEURON_DEATH,
                       go.pathways.msigdb$GO_RESPONSE_TO_OXIDATIVE_STRESS,go.pathways.msigdb$GO_RESPONSE_TO_REACTIVE_OXYGEN_SPECIES, 
                       go.pathways.msigdb$GO_REACTIVE_OXYGEN_SPECIES_METABOLIC_PROCESS, go.pathways.msigdb$GO_NEURON_DEATH_IN_RESPONSE_TO_OXIDATIVE_STRESS, go.pathways.msigdb$GO_RESPONSE_TO_REACTIVE_OXYGEN_SPECIES,  
                       go.pathways.msigdb$GO_INTERMEDIATE_FILAMENT_CYTOSKELETON,go.pathways.msigdb$GO_INTERMEDIATE_FILAMENT_BASED_PROCESS, go.pathways.msigdb$GO_ACTIN_FILAMENT, go.pathways.msigdb$GO_CYTOSKELETON_ORGANIZATION)

astrocyte.reactivity.splicing.labels = c("VIM", "COL1A1", "MYL6", "DDX5", "FN1", "FLNA", "FUS", "SFPQ","OSMR", "NONO", "HLA-A", "HLA-B", "HLA-C", "STAT3", "RND3", "SLC26AA6","PIMREG", "RRP8", "RECK", "LOXL3","UBN1", "COL7A1","MX1",
                                "TTBK1", "TSC1", "MBD1", "EMP1", "PRPF4", "PIF1", "AARS2", "PLSCR1", "NARPT", "CSRP1", "IRF7", "ID4", "AGER", "GBP2", "ANKZF1", "AGRN", "GTBP2", "TTBK1", "CCS", "TGFB1I1","COL1A1A", "TNC", 
                                "RPL10", "NUP199", "LAMB2", "IDH1", "TP53", "LRP1")
astrocyte.reactivity.labels = c("VIM", "COL1A1", "MYL6", "FN1", "FLNA", "OSMR", "HLA-A", "HLA-B", "HLA-C", "STAT3", "RND3", "SLC26AA6","PIMREG", "RRP8", "RECK", "LOXL3","UBN1", "COL7A1","MX1",
                              "TTBK1", "TSC1", "MBD1", "EMP1", "PRPF4", "PIF1", "AARS2", "PLSCR1", "NARPT", "CSRP1", "IRF7", "ID4", "AGER", "GBP2", "ANKZF1", "AGRN", "GTBP2", "TTBK1", "CCS", "TGFB1I1","COL1A1A", "TNC", 
                              "RPL10", "NUP199", "LAMB2", "IDH1", "TP53", "LRP1",
                              "TINF2", "ILK", "MCAM", "FN1", "FBLN5", "CAPN2", "TGFB1I1", "CDKN1A", "ADAM19", "SERPINE1", "ITGA6", "ADAMTSL4", "ITGA5", "PDLIM7", "FLOT1", "DIXDC1", "BAG3", "PVR")  # overlap increased exp & decreased IR
astrocyte.reactivity.labels = c(astrocyte.reactivity.labels, astrocyte.reactive.genes)

astrocyte.reactivity.labels.ac_nuc <- c("VIM", "COL1A1", "MYL6", "FN1", "FLNA", "OSMR", "HLA-A", "HLA-B", "HLA-C", "STAT3", "RND3", "SLC26AA6",
                                         "PIMREG", "RRP8", "RECK", "MX1", "TTBK1", "MBD1", "EMP1", "PRPF4", "PLSCR1", "NARPT", "CSRP1", "IRF7", "ANKZF1", "GTBP2",   
                                         "CCS", "TGFB1I1", "COL1A1A", "RPL10", "NUP199", "LAMB2", "IDH1", "LRP1", "TINF2", "ILK" , "MCAM", "CAPN2", "ADAM19", "SERPINE1", "ITGA6", "PDLIM7", "FLOT1",  
                                         "BAG3", "LCN2", "STEAP4", "S1PR3", "TIMP1", "HSPB1", "CXCL10", "CD44", "CP", "GFAP", "SERPINA3", "ASPG", "C3", "HLA-F", "MICA", "H2-T23", "H2-D1", "AMIGO2", "SERPING1", "GGTA1P", "GGTA1", "UGT1A1",  
                                         "FKBP5", "PSMB8", "SRGN", "IIGP1", "S100A10", "CLCF1", "TGM1", "PTX3", "SPHK1", "CD109", "PTGS2", "SLC10A6", "TM4SF1", "B3GNT5", "CD14")  

astrocyte.reactivity.labels.ac_cyt <- c("VIM", "COL1A1", "FN1", "FLNA", "OSMR", "HLA-A", "HLA-B", "HLA-C", "STAT3", "RND3", "SLC26AA6","PIMREG", "RRP8", "RECK", "LOXL3", "COL7A1", "MX1", "EMP1",
  "NARPT", "CSRP1", "ID4", "AGER", "GBP2", "ANKZF1", "AGRN", "GTBP2",   "CCS", "TGFB1I1", "COL1A1A", "TNC", "NUP199", "LAMB2", "IDH1", "TP53", "ILK" ,    
  "MCAM", "FBLN5", "CAPN2", "CDKN1A", "ADAM19", "SERPINE1", "ADAMTSL4", "ITGA5", "PDLIM7", "FLOT1", "DIXDC1",  "BAG3", "LCN2", "STEAP4", "S1PR3", "TIMP1", "HSPB1", "CXCL10", "CD44", "CP", "GFAP", "SERPINA3",
  "ASPG", "C3", "HLA-E", "MICA", "H2-T23", "H2-D1", "AMIGO2", "SERPING1", "GGTA1P", "GGTA1", "UGT1A1",  "FKBP5", "SRGN", "IIGP1", "S100A10", "CLCF1", "TGM1", "PTX3", "SPHK1", "CD109", "PTGS2", "SLC10A6", "TM4SF1", "B3GNT5", "CD14")  

astrocyte.reactivity.labels.ac_nucIR_cytDGE <- c("VIM", "COL1A1", "FN1", "HLA-A", "HLA-B", "HLA-C", "SLC26AA6","PIMREG", "RRP8", "RECK", "MX1", "TTBK1", "EMP1", "PRPF4",   
  "NARPT", "CSRP1", "IRF7", "GBP2", "ANKZF1", "GTBP2","CCS", "TGFB1I1", "COL1A1A", "RPL10", "NUP199", "TP53", "TINF2", "ILK",    
  "MCAM", "CAPN2", "ADAM19", "SERPINE1", "ITGA6", "ADAMTSL4", "FLOT1","BAG3", "PVR", "LCN2", "STEAP4", "S1PR3", "TIMP1", "HSPB1", "CXCL10", "CD44", "CP", "GFAP", "SERPINA3",
  "ASPG", "C3", "HLA-E", "HLA-F", "MICA", "H2-T23", "H2-D1", "AMIGO2", "SERPING1", "GGTA1P", "GGTA1", "UGT1A1",  "FKBP5", "PSMB8", "SRGN", "IIGP1", "S100A10", "CLCF1", "TGM1", "PTX3", "SPHK1", "CD109", "PTGS2", "SLC10A6", "TM4SF1", "B3GNT5", "CD14")  

astrocyte.reactivity.labels.ac_c9orf72 <- c("VIM", "COL1A1", "FN1", "FLNA", "OSMR", "HLA-A", "HLA-B", "HLA-C", "STAT3", "RND3", "SLC26AA6","PIMREG", "RRP8", "RECK", "UBN1", "MX1", "TTBK1", "TSC1", "MBD1", "EMP1", "PRPF4",   
  "PIF1", "AARS2", "PLSCR1", "NARPT", "CSRP1", "IRF7", "ID4", "AGER", "GBP2", "ANKZF1", "AGRN", "GTBP2","CCS", "TGFB1I1", "COL1A1A", "TNC", "RPL10", "NUP199", "LAMB2", "IDH1", "TP53", "LRP1", "TINF2", "ILK" ,    
  "MCAM", "FBLN5", "CAPN2", "CDKN1A", "ADAM19", "SERPINE1", "ITGA6", "ADAMTSL4", "ITGA5", "PDLIM7", "FLOT1", "DIXDC1", "BAG3", "PVR", "LCN2", "STEAP4", "S1PR3", "TIMP1", "HSPB1", "CXCL10", "CD44", "CP", "GFAP", "SERPINA3",
  "ASPG", "C3", "HLA-E", "HLA-F", "MICA", "H2-T23", "H2-D1", "AMIGO2", "SERPING1", "GGTA1P", "GGTA1", "UGT1A1",  "FKBP5", "PSMB8", "SRGN", "IIGP1", "S100A10", "CLCF1", "TGM1", "SPHK1", "CD109", "PTGS2", "SLC10A6", 
  "TM4SF1", "B3GNT5", "CD14")  

astrocyte.reactivity.labels.ac_sod1 <-  c("COL1A1", "MYL6", "FN1", "FLNA", "OSMR", "HLA-A", "HLA-B", "HLA-C", "STAT3", "RND3", "SLC26AA6",
  "PIMREG", "RRP8", "RECK", "LOXL3", "UBN1", "COL7A1", "MX1", "TTBK1", "TSC1", "MBD1", "EMP1", "PRPF4",   
  "PIF1", "AARS2", "PLSCR1", "NARPT", "CSRP1", "IRF7", "ID4", "AGER", "GBP2", "ANKZF1", "AGRN", "GTBP2",   
  "CCS", "TGFB1I1", "COL1A1A", "TNC", "RPL10", "NUP199", "LAMB2", "IDH1", "TP53", "LRP1", "TINF2", "ILK" ,    
  "FBLN5", "CAPN2", "CDKN1A", "ADAM19", "SERPINE1", "ITGA6", "ADAMTSL4", "ITGA5", "PDLIM7", "DIXDC1",  
  "BAG3", "PVR", "LCN2", "STEAP4", "S1PR3", "TIMP1", "HSPB1", "CXCL10", "CD44", "CP", "GFAP", "SERPINA3",
  "ASPG", "C3", "HLA-E", "HLA-F","H2-T23", "H2-D1", "AMIGO2", "SERPING1", "GGTA1P", "GGTA1", "UGT1A1",  
  "FKBP5", "PSMB8", "SRGN", "IIGP1", "S100A10", "CLCF1", "TGM1", "PTX3", "SPHK1", "CD109", "PTGS2", "SLC10A6", 
  "TM4SF1", "B3GNT5", "CD14")  
astrocyte.reactivity.labels.ac_a1 <- c("FN1", "ILK", "OSMR", "TINF2", "IRF7", "MCAM", "SOD2", "PSMB7", "PSMB10", "DNAJC2", "KMT2E", "COL27A1", "ERF")

astrocyte.reactivity.labels.ac_tdp43 <- c("FLOT1", "MCAM", "MBD1", "ANKZF1", "AARS2", "PRRT1", "RAB34", "DUSP1", "SPG11", "ABCD1", "CFP", "TRIP6", "ADAMTS20", "PTPN6", "ITGA10", "ADAMTS9", "BCAN", "NSMF", "MAP7", "PLOD3")

astrocyte.reactivity.labels.ac_protein <- c("FN1", "ILK", "OSMR", "TINF2", "IRF7", "MCAM", "SOD2", "PSMB7", "PSMB10", "DNAJC2", "KMT2E", "COL27A1", "ERF", "ATP6V1F", "RFTN1", "IDH1", "HNRNPU", "VCP", "ACTB", "COL5A2", "HSPB1", "PDLIM2", "PURA", "SVIL", "RAP1A", "TGFB1I1", "HLA-E", "ITGA3", "RELA", "TUBA1C")
RNA.capping <-c("POLR2J","MNAT1","POLR2B","NCBP3","POLR2E","POLR2F","RNMT","POLR2C","ERCC2","POLR2I","GTF2H1","GTF2H3","RNGTT","NCBP2","GTF2F1","CDK7","CCNH","NCBP1","CMTR1","TGS1","POLR2D","GTF2H2","POLR2K","ERCC3","POLR2H","POLR2G","RAMAC","POLR2L","CMTR2","POLR2A","GTF2F2","SUPT5H","GTF2H4","RAMACL","GTF2H5")
exon.junction.complex <- c("UPF1","THRAP3","SMG6","TDRD3","PNN","R3HCC1","CASC3","MAGOHB","UPF3B","SRSF1","EIF4A3","SAP18","UPF2","MAGOH","R3HCC1L","UPF3A","PYM1","RNPS1","RBM8A","SRRM1","PRPF8","EIF4A3","MLN51","MAGOH","Y14","CASC3","EIF4A3","MAGOH","RBM8A","PYM1","NXF1","DDX39B")
RNA.polyadenylation <- c("PAF1","ZC3H3","PABPC1","CPSF1","SNRPA","YTHDC1","PAPOLA","ZC3H14","PABPC1L","CSTF2","RNF40","CPEB3","MTPAP","CPSF6","TENT4A","PAPOLG","CPSF3","TENT4B","SYMPK","CCNT1","GRSF1","CDC73","WDR33","CDK9","PNPT1","APP","FIP1L1","TUT1","CPSF7","RNF20","SSU72","TENT2","VIRMA","PCF11","CPSF2","LEO1","NUDT21","AHCYL1","CLP1","CSTF3","HSF1","SUPT5H","CTR9","NELFE","SCAF8","CPEB1","PAPOLB","SSU72P8","AP002990.1","SSU72P4","SSU72P5","SSU72P2","SSU72P7","SSU72P3","NCBP2","NCBP1","","CPSF4","NUDT21","CSTF2T","CPSF4")
# Transport of Mature Transcript to Cytoplasm	NUP160	TPR	THOC3	ZC3H11A	NDC1	U2AF2	NUP133	CPSF1	NUP37	THOC1	NUP50	AAAS	NUP188	POLDIP3	THOC5	SRSF5	RAE1	LUZP4	NUP93	CASC3	NUP88	MAGOHB	NUP107	SRSF9	SRSF3	NUP155	NCBP2	SRSF7	SRSF4	SRSF11	CPSF3	GLE1	NUP43	FYTTD1	DDX39A	SRSF6	NUP153	UPF3B	NUP85	THOC2	SYMPK	NUP214	THOC6	NUP210	NXT1	SRRM1	NUP42	SRSF1	WDR33	NCBP1	NUP54	DHX38	EIF4A3	FIP1L1	EIF4E	RANBP2	NUP205	SEC13	U2AF1	CHTOP	CPSF4	U2AF1L4	SRSF2	NXF1	MAGOH	NUP35	THOC7	SLBP	SLU7	CPSF2	CDC40	POM121	DDX39B	SARNP	RNPS1	NUP62	RBM8A	NXF2	NXF2B	POM121C
minor.spliceosome.go <- c("RNU4ATAC","RNU6ATAC","PHF5A","SF3B1","SF3B2","SF3B3","SF3B4","SF3B5","SF3B6","SNRPB","SNRPD1","SNRPD2","SNRPD3","SNRPE","SNRPF","SNRPG","SNRPN","DHX15","PDCD7","RNPC3","RNU11","RNU12","SNRNP25","SNRNP35","SNRNP48","YBX1","ZCRB1","ZMAT5","ZRSR2","CD2BP2","DDX23","EFTUD2","PRPF6","PRPF8","RNU5A-1","SNRNP200","SNRNP40","TXNL4A")
minor.splicing.go <- unique(c(gprofiler_database$`mRNA Splicing - Minor Pathway`, minor.spliceosome.go))
RNA.capping <-c("POLR2J","MNAT1","POLR2B","NCBP3","POLR2E","POLR2F","RNMT","POLR2C","ERCC2","POLR2I","GTF2H1","GTF2H3","RNGTT","NCBP2","GTF2F1","CDK7","CCNH","NCBP1","CMTR1","TGS1","POLR2D","GTF2H2","POLR2K","ERCC3","POLR2H","POLR2G","RAMAC","POLR2L","CMTR2","POLR2A","GTF2F2","SUPT5H","GTF2H4","RAMACL","GTF2H5")
exon.junction.complex <- c("UPF1","THRAP3","SMG6","TDRD3","PNN","R3HCC1","CASC3","MAGOHB","UPF3B","SRSF1","EIF4A3","SAP18","UPF2","MAGOH","R3HCC1L","UPF3A","PYM1","RNPS1","RBM8A","SRRM1","PRPF8","EIF4A3","MLN51","MAGOH","Y14","CASC3","EIF4A3","MAGOH","RBM8A","PYM1","NXF1","DDX39B")
RNA.polyadenylation <- c("PAF1","ZC3H3","PABPC1","CPSF1","SNRPA","YTHDC1","PAPOLA","ZC3H14","PABPC1L","CSTF2","RNF40","CPEB3","MTPAP","CPSF6","TENT4A","PAPOLG","CPSF3","TENT4B","SYMPK","CCNT1","GRSF1","CDC73","WDR33","CDK9","PNPT1","APP","FIP1L1","TUT1","CPSF7","RNF20","SSU72","TENT2","VIRMA","PCF11","CPSF2","LEO1","NUDT21","AHCYL1","CLP1","CSTF3","HSF1","SUPT5H","CTR9","NELFE","SCAF8","CPEB1","PAPOLB","SSU72P8","AP002990.1","SSU72P4","SSU72P5","SSU72P2","SSU72P7","SSU72P3","NCBP2","NCBP1","","CPSF4","NUDT21","CSTF2T","CPSF4")
exosome.complex <- c("MTREX","EXOSC7","EXOSC5","DIS3","KHSRP","GTPBP1","EXOSC3","EXOSC8","EXOSC9","EXOSC2","WDR74","ZFC3H1","MPHOSPH6","PNPT1","NVL","DIS3L2","CARHSP1","SUPV3L1","DIS3L","EXOSC1","EXOSC10","EXOSC4","C1D","EXOSC6")
# response.to.typeIinterferon <- c("SP100","IFI35","IP6K2","UBE2K","MAVS","OAS1","MUL1","HSP90AB1","TTLL12","SAMHD1","CACTIN","TYK2","CDC37","OAS3","OAS2","PTPN6","WNT5A","STAT1","IRF6","IFIT3","IFIT2","IFNA6","IFNA8","EGR1","TRIM6","ZBP1","IRF1","IRF3","IFI6","IRF5","BST2","SHFL","XAF1","RSAD2","OASL","RNASEL","IFNA21","IRF4","NLRC5","IRF8","IFITM3","IFNAR1","IFNA5","IFNA16","MX1","IFNAR2","ADAR","JAK1","GBP2","DCST1","ABCE1","METTL3","IFI27","FADD","IRF2","TRIM56","STAT2","IFNB1","ISG20","MYD88","PTPN2","PTPN11","SETD2","SHMT2","MX2","TBK1","IRAK1","USP18","IFITM2","IRF7","YTHDF3","IFIT1","IFITM1","IFNA10","ISG15","IFNA2","PTPN1","IFNA1","YTHDF2","CNOT7","PSMB8","HLA-C","HLA-E","HLA-G","HLA-F","HLA-A","TREX1","IRF9","IFNA7","IFNA14","IFNA13","HLA-B","IFNA17","IFNA4","LSM14A","MMP12","IKBKE","MIR21")	
# response.to.typeIinterferon <- c("SP100","IFI35","IP6K2","UBE2K","MAVS","OAS1","MUL1","HSP90AB1","TTLL12","SAMHD1","CACTIN","TYK2","CDC37","OAS3","OAS2","PTPN6","WNT5A","STAT1","IRF6","IFIT3","IFIT2","IFNA6","IFNA8","EGR1","TRIM6","ZBP1","IRF1","IRF3","IFI6","IRF5","BST2","SHFL","XAF1","RSAD2","OASL","RNASEL","IFNA21","IRF4","NLRC5","IRF8","IFITM3","IFNAR1","IFNA5","IFNA16","MX1","IFNAR2","ADAR","JAK1","GBP2","DCST1","ABCE1","METTL3","IFI27","FADD","IRF2","TRIM56","STAT2","IFNB1","ISG20","MYD88","PTPN2","PTPN11","SETD2","SHMT2","MX2","TBK1","IRAK1","USP18","IFITM2","IRF7","YTHDF3","IFIT1","IFITM1","IFNA10","ISG15","IFNA2","PTPN1","IFNA1","YTHDF2","CNOT7","PSMB8","HLA-C","HLA-E","HLA-G","HLA-F","HLA-A","TREX1","IRF9","IFNA7","IFNA14","IFNA13","HLA-B","IFNA17","IFNA4","LSM14A","MMP12","IKBKE","MIR21")	
# response.interferon.gamma	CX3CL1	CCL26	NOS2	NUB1	WAS	SLC11A1	NR1H3	VIM	CD44	IFNGR1	PIAS1	VAMP3	SEC61A1	CAMK2B	SP100	CAMK2A	CDC42	TXK	STXBP2	EDN1	OAS1	KIF16B	ICAM1	HSP90AB1	JAK2	IL12RB1	CD40	HCK	MID1	ACOD1	CCL22	CCL17	MEFV	STX4	CDC37	CCL24	CCL7	CCL2	CCL8	CCL1	HPX	CYP27B1	OAS3	OAS2	IFNG	GAPDH	TRIM38	IL12B	WNT5A	CCL20	ACTR3	STAT1	KYNU	SUMO1	STXBP3	TRIM62	CD58	GBP3	GBP1	FASLG	IRF6	ARG1	TRIM25	NMI	MT2A	IRF1	MED1	CITED1	IRF3	IRF5	BST2	ASS1	SHFL	CCL25	NR1H2	GCH1	TRIM21	PPARG	TRIM5	TRIM22	OASL	SYNCRIP	FLNB	CALCOCO2	EPRS1	STXBP1	TLR4	MYC	CCL21	IRF4	TLR2	CASP1	ACTR2	PARP9	RAB20	PML	NLRC5	IRF8	TP53	IFITM3	RPL13A	XCL1	XCL2	SNCA	CAMK2D	STAR	GSN	CAMK2G	NCAM1	CDC42EP2	FCGR1A	VPS26B	GBP5	IFNGR2	ZYX	ADAMTS13	CXCL16	JAK1	IL23R	GBP2	GBP4	VCAM1	CLDN1	PRKCD	OTOP1	TLR3	SLC30A8	STXBP4	B2M	TRIM68	DAPK3	IRF2	LGALS9	PTAFR	STX8	KIF5B	TRIM8	AQP4	CCL11	CCL19	RAB43	PARP14	PDE12	PTPN2	CIITA	CDC42EP4	CCL13	GBP6	ACTG1	SOCS3	IFITM2	SOCS1	IRF7	IFITM1	HLA-DRB1	EVL	DAPK1	HLA-DQA1	CD47	MYO1C	FCGR1B	SIRPA	HLA-DRB5	HLA-DRA	AIF1	HLA-C	HLA-E	TRIM31	HLA-G	HLA-F	RAB12	HLA-A	GBP7	UBD	IRF9	IFI30	HLA-DPB1	SLC26A6	HLA-DPA1	HLA-DQB2	TRIM26	HLA-B	HLA-DQA2	IRGM	TDGF1	TRIM34	MRC1	CCL5	CCL23	CCL16	CCL4	CCL18	CCL15-CCL14	CCL15	CCL4L2	CCL3L1	CCL14	RAB7B	CCL3
# response.interferon.alpha	TPR	EIF2AK2	KLHL20	LAMP3	GATA3	RO60	IFIT3	IFIT2	BST2	IFITM3	IFNAR1	STAR	IFNAR2	ADAR	PYHIN1	AXL	PDE12	GAS6	MX2	IFITM2	IFITM1
# response.interferon.beta	UBE2K	CDC34	ACOD1	HTRA2	STAT1	TRIM6	IRF1	BST2	SHFL	XAF1	PNPT1	IFITM3	IFNAR2	CAPN2	MNDA	PYHIN1	IFI16	AIM2	TLR3	IFNB1	DDX41	STING1	UBE2G2	IFITM2	IFITM1	NDUFA13	PLSCR1	TREX1	IRGM	IKBKE
# regulation of response to interferon-gamma	NR1H3	IFNGR1	PIAS1	TXK	HSP90AB1	JAK2	CDC37	HPX	IFNG	STAT1	SUMO1	ARG1	MED1	NR1H2	PPARG	PARP9	NLRC5	IFNGR2	JAK1	OTOP1	PARP14	PTPN2	SOCS3	SOCS1	IRGM
# GO:0060331	negative regulation of response to interferon-gamma	NR1H3	ARG1	NR1H2	PPARG	OTOP1	PARP14	PTPN2
# GO:0060332	positive regulation of response to interferon-gamma	TXK	HPX	MED1	PARP9	NLRC5	IRGM
# GO:0060333	interferon-gamma-mediated signaling pathway	NR1H3	CD44	IFNGR1	PIAS1	CAMK2B	SP100	CAMK2A	TXK	OAS1	ICAM1	HSP90AB1	JAK2	HCK	MID1	CDC37	HPX	OAS3	OAS2	IFNG	TRIM38	STAT1	SUMO1	TRIM62	GBP1	IRF6	ARG1	TRIM25	NMI	MT2A	IRF1	MED1	IRF3	IRF5	NR1H2	TRIM21	PPARG	TRIM5	TRIM22	OASL	IRF4	PARP9	PML	NLRC5	IRF8	TP53	CAMK2D	CAMK2G	NCAM1	FCGR1A	IFNGR2	JAK1	GBP2	VCAM1	PRKCD	OTOP1	B2M	TRIM68	IRF2	PTAFR	TRIM8	PARP14	PTPN2	CIITA	SOCS3	SOCS1	IRF7	HLA-DRB1	HLA-DQA1	FCGR1B	HLA-DRB5	HLA-DRA	HLA-C	HLA-E	TRIM31	HLA-G	HLA-F	HLA-A	IRF9	IFI30	HLA-DPB1	HLA-DPA1	HLA-DQB2	TRIM26	HLA-B	HLA-DQA2	IRGM	TRIM34
# GO:0060334	regulation of interferon-gamma-mediated signaling pathway	NR1H3	IFNGR1	PIAS1	TXK	HSP90AB1	JAK2	CDC37	HPX	IFNG	STAT1	SUMO1	ARG1	MED1	NR1H2	PPARG	PARP9	NLRC5	IFNGR2	JAK1	OTOP1	PARP14	PTPN2	SOCS3	SOCS1	IRGM
# GO:0060335	positive regulation of interferon-gamma-mediated signaling pathway	TXK	HPX	MED1	PARP9	NLRC5	IRGM
# GO:0060336	negative regulation of interferon-gamma-mediated signaling pathway	NR1H3	ARG1	NR1H2	PPARG	OTOP1	PARP14	PTPN2
# GO:0060337	type I interferon signaling pathway	SP100	IFI35	IP6K2	UBE2K	MAVS	OAS1	MUL1	HSP90AB1	TTLL12	SAMHD1	CACTIN	TYK2	CDC37	OAS3	OAS2	PTPN6	WNT5A	STAT1	IRF6	IFIT3	IFIT2	IFNA6	IFNA8	EGR1	TRIM6	ZBP1	IRF1	IRF3	IFI6	IRF5	BST2	XAF1	RSAD2	OASL	RNASEL	IFNA21	IRF4	NLRC5	IRF8	IFITM3	IFNAR1	IFNA5	IFNA16	MX1	IFNAR2	ADAR	JAK1	GBP2	DCST1	ABCE1	METTL3	IFI27	FADD	IRF2	STAT2	IFNB1	ISG20	MYD88	PTPN2	PTPN11	MX2	TBK1	IRAK1	USP18	IFITM2	IRF7	YTHDF3	IFIT1	IFITM1	IFNA10	ISG15	IFNA2	PTPN1	IFNA1	YTHDF2	CNOT7	PSMB8	HLA-C	HLA-E	HLA-G	HLA-F	HLA-A	TREX1	IRF9	IFNA7	IFNA14	IFNA13	HLA-B	IFNA17	IFNA4	LSM14A	MMP12	IKBKE	MIR21
# GO:0060338	regulation of type I interferon-mediated signaling pathway	UBE2K	MAVS	MUL1	HSP90AB1	TTLL12	SAMHD1	CACTIN	CDC37	PTPN6	WNT5A	TRIM6	ZBP1	IRF3	RNASEL	NLRC5	IFNAR2	ADAR	DCST1	ABCE1	METTL3	FADD	PTPN2	PTPN11	TBK1	USP18	IRF7	YTHDF3	PTPN1	YTHDF2	CNOT7	TREX1	LSM14A	MMP12	IKBKE	MIR21
# negative regulation of type I interferon-mediated signaling pathway	MUL1	TTLL12	SAMHD1	CACTIN	NLRC5	ADAR	DCST1	METTL3	PTPN2	YTHDF3	YTHDF2	CNOT7	TREX1	MMP12	MIR21
# positive regulation of type I interferon-mediated signaling pathway	UBE2K	MAVS	WNT5A	TRIM6	ZBP1	IRF3	NLRC5	FADD	TBK1	IRF7	LSM14A	MMP12	IKBKE

astrocyte.markers <- unique(c("S100B", "SOX9", "FGFR3", "GFAP", "AQP4", "VIM", "EDNRB","NFIX","PLA2G7","SLC39A12","MLC1","AU067665","DIO2","SLC14A1","ALDH1L1","CYP4F3","ALDOC","TTPA","ACSBG1","CHRDL1","GM266","SLC4A4","SLC1A2","SLC25A18","SLC1A3","F3","PPP1R3G","CYP4F12","LOC388419","FZD2","2900019G14RIK","MERTK","EZR","EVA1A","GJB6","HAPLN1","RFX4","PAPSS2","SLC15A2","PPP1R3C","TLR3","ACOT11","ATP1A2","BMPR1B","DEPDC2","PRODH","GLI3","TMEM47","SLC9A3R1","CTH","NTSR2","SLC7A10","VCAM1","FGFR3","CCDC80","ENTPD2","CYBRD1","KCNE5","FAM20A","KIAA1161","EG328479","TNC","TLCD1","S1PR1","CBS","PBXIP1","GRIN2C","A730056I06RIK","ADHFE1","AGT","GLDC","SLC7A2","BC055107","GJA1","PDK4","EGFR","SOX9","CLDN10","PLCD4","ID4","FMO1","EMP2","LONRF3","HTRA1","MGST1", "THRSP"))
motor.neuron.markers <- unique(c("OLIG2", "NEUROG2", "ISL1", "ISL2", "CHAT", "FOXP1", "LHX3", "HOXA5", "POU3F1", "TSHZ1", "PHOX2B", "OTX2", "MAFB", "NEFH", "TUBB3", "MAP2", "NCAM1", "NEUROD6","GLRA2","CCN3","PRDM8","SLA","CRH","GABRG2","HTR2C","HS3ST2","MAL2","9130024F11RIK","NECAB1","STMN2","GABRA5","NTS","A130090K04RIK","GABRA1","SATB2","GPR88","SYT1","4930544G21RIK","GDA","MYT1L","SLC17A6","A930034L06RIK","A930009L07RIK","NPAS4","CALB1","SLC12A5","EPHA7","VIP","MEF2C","SSTR2","TENM2","PCSK2","SNAP25","SCG2","PGM2L1","PLCXD3","DLX1","VSNL1","SYT4","NRG3","SCN2A","ICAM5","KCNF1","AI838057","CCK","VGF","TMEM130","CAMK4","ASPH","C030017B01RIK","SLC6A7","OLFR1344","ICA1L","MYO5B","NELL1","NEFM","NEFL","TTR","6330527O06RIK","CDH8","SV2B","GAP43","TRHDE","D11BWG0517E","CAMK2B","PENK","RGS4","LPL","CACNA1B","KCNC2","TTC9","L1CAM","9130430L19RIK","CLSTN2","NAPB","CYB561","NKX6-1","NKX6-2","CXADR","GLI2","NKX2-2","FOXA2","CHAT","MNX1","ISL1","FOXP1","LHX3","ETV4","LMO4","NOTCH1","SIM1"))
ac.mn.gene.markers <- c(astrocyte.markers, motor.neuron.markers)

fractionation_labels = c("MALAT1", "GAPDH", "NEAT", "MALAT1", "GFAP", "TUBB", "ASH2L", "PAF49", "CENPA", "CHOP", "EIF6", "ERS1", "H2AX", "H2AZ1", "H4F3", "H3F3B", "LSD1", "LMNA", "LMNC", "MYC", "NOP2", "NUP98", "SIR1", "SRSF2")
fractionation_labels_filt = c("MLXIPL", "MALAT1", "GAPDH", "NEAT1", "GLUL", "ACTB", "ASS1", " SLC2A2", "APOE", "PCK1", "MALAT1", "GFAP", "TUBB", "PAF49", "CHOP", "ERS1", "H4F3", "LSD1", "HIST", "LMNC", "NOP2", "NUP98", "SIR1", "HIST4H4", "HIST1H2APS3", "HIST1H2APS5", "HIST2H2BD", "HIST2H2BF", "HIST2H2BK", "HIST2H2AK", "HIST1H2PS3") #NLRP6 GCGR


phast <- phastCons100way.UCSC.hg38
gtf <- readRDS("/camp/home/ziffo/home/genomes/ensembl/Homo_sapiens.GRCh38.99.gtf.coords.RDS")
# gene_coordinates <- gtf %>% dplyr::mutate(gene_id = ensemblID, gene_coords = paste(.$seqnames, ":", .$start, "-", .$end, sep = "")) %>% dplyr::select(gene_coords, gene_id)

# RNAInter database
# sbatch -N 1 -c 4 --mem=10G -t 1:00:00 --wrap="wget 'http://www.rna-society.org/rnainter/download/RNA-Protein.zip'" --mail-type=ALL,ARRAY_TASKS --mail-user=oliver.ziff@crick.ac.uk --job-name=RNAinter
# sbatch -N 1 -c 4 --mem=10G -t 2:00:00 --wrap="wget 'http://www.rna-society.org/rnainter/download/Full%20Version.zip'" --mail-type=ALL,ARRAY_TASKS --mail-user=oliver.ziff@crick.ac.uk --job-name=RNAinterFULL
# sbatch -N 1 -c 4 --mem=10G -t 1:00:00 --wrap="unzip 'Full Version.zip'" --mail-type=ALL,ARRAY_TASKS --mail-user=oliver.ziff@crick.ac.uk --job-name=RNAinter
# sbatch -N 1 -c 4 --mem=10G -t 1:00:00 --wrap="unzip 'RNA-Protein.zip'" --mail-type=ALL,ARRAY_TASKS --mail-user=oliver.ziff@crick.ac.uk --job-name=RNAinter
# RNAInter_interaction_full <- read.table("/camp/lab/luscomben/home/users/ziffo/genomes/rbp-databases/rnainter/RNAInter_interaction_full.txt", header=FALSE, sep = "\t")
# RNAInter_RNA_Protein.txt <- read.table("/camp/lab/luscomben/home/users/ziffo/genomes/rbp-databases/rnainter/RNA-Protein.txt", header=FALSE, sep = "\t") %>% as_tibble() %>% select(gene_name = V2, gene_id = V3, biotype = V4, gene.organism = V5, protein_name = V6, protein_id = V7, protein_type = V8, protein.organism = V9, score = V10)
# RNAInter_RNA_Protein.human <- RNAInter_RNA_Protein.txt %>% filter(gene.organism == "Homo sapiens", protein.organism == "Homo sapiens")
# saveRDS(RNAInter_RNA_Protein.human, "/camp/lab/luscomben/home/users/ziffo/genomes/rbp-databases/rnainter/rna_protein.human.rds")
# RNAInter_RNA_Protein.human <- readRDS("/camp/lab/luscomben/home/users/ziffo/genomes/rbp-databases/rnainter/rna_protein.human.rds")
# Sys.setenv("http_proxy" = "http://my.proxy.org:9999") # options(RCurlOptions = list(proxy="uscache.kcc.com:80",proxyuserpwd="------:-------"))
# httr::set_config(httr::config(ssl_verifypeer = FALSE)) # httr::set_config(httr::config(ssl_cipher_list = "DEFAULT@SECLEVEL=1"))
# ensembl99 <- useEnsembl(biomart = 'genes', dataset = 'hsapiens_gene_ensembl', version = 99) # Select biomaRt database (genes) and dataset (hsapiens)
# ens2entrez <- getBM(attributes =  c("ensembl_gene_id",  "entrezgene_id"),values = TRUE,mart = ensembl99,useCache = FALSE) %>% select(ensemblID = ensembl_gene_id, entrezID = entrezgene_id) # # NCBI entrez gene ID to Ensembl Gene ID
# gtf <- gtf %>% left_join(ens2entrez)
# RNAInter_RNA_Protein <- RNAInter_RNA_Protein.human %>%  filter(!gene_name %in% gtf$gene_name | gene_id %in% gtf$ensemblID | gene_id %in% gtf$entrezID) # 1,990,466 interactions can be matched based on gene identifiers (many missing e.g. ENSG00000124835)
# RNAInter_RNA_RBP <- RNAInter_RNA_Protein.human %>%  filter(protein_name %in% RNA_BINDING_PROTEINS) %>% # 22,501,074 interactions, 1,471 RBPs
  # # filter(protein_type == "RBP")  # some important RBPs are annotated as TFs eg FUS and SFPQ so cant rely on protein_type == 192 RBPs and 36,280 genes
  # filter(!gene_name %in% gtf$gene_name | gene_id %in% gtf$ensemblID | gene_id %in% gtf$entrezID) # 1,990,466 interactions can be matched based on gene identifiers (many missing e.g. ENSG00000124835)
# saveRDS(RNAInter_RNA_RBP, "/camp/lab/luscomben/home/users/ziffo/genomes/rbp-databases/rnainter/rna_rbp.human.rds")
RNAInter_RNA_RBP <- readRDS("/camp/lab/luscomben/home/users/ziffo/genomes/rbp-databases/rnainter/rna_rbp.human.rds")


# Generic Functions ------------------
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}

transpose_tibble <- function(tibble,old_rownames,new_rownames){
  tibble.t <- tibble %>% column_to_rownames(old_rownames) %>% t %>% as_tibble(rownames = new_rownames)
  tibble.t
}
# transpose_tibble(vsd_counts.phag_glut, old_rownames = "gene_name", new_rownames = "group")

scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)


left_join_NA <- function(x, y, ...) {
  left_join(x = x, y = y, by = ...) %>% 
    mutate_each(funs(replace(., which(is.na(.)), 0)))
}


# Tx2ens and gene2ens annotation - for HTSeq, Salmon, Kallisto
# mus_musculus.gtf <- import("/camp/home/ziffo/home/genomes/ensembl/Mus_musculus.GRCm38.101.gtf") %>% as_tibble
# saveRDS(mus_musculus.gtf, "/camp/home/ziffo/home/genomes/ensembl/Mus_musculus.GRCm38.101.gtf.RDS")
# mus_musculus.gtf <- readRDS("/camp/home/ziffo/home/genomes/ensembl/Mus_musculus.GRCm38.101.gtf.RDS")
# mus.musculus.gene2ens <- mus_musculus.gtf %>% dplyr::select(gene_id, gene_name) %>% unique()
# saveRDS(mus.musculus.gene2ens, "/camp/home/ziffo/home/genomes/ensembl/Mus_musculus.GRCm38.101.gene2ens.RDS")
mus.musculus.gene2ens <- readRDS("/camp/home/ziffo/home/genomes/ensembl/Mus_musculus.GRCm38.101.gene2ens.RDS")
# mus.musculus.tx2gene <- mus_musculus.gtf %>% dplyr::select(transcript_id, gene_id) %>% unique() %>% drop_na()
# saveRDS(mus.musculus.tx2gene, "/camp/home/ziffo/home/genomes/ensembl/Mus_musculus.GRCm38.101.tx2gene.RDS")
mus.musculus.tx2gene <- readRDS("/camp/home/ziffo/home/genomes/ensembl/Mus_musculus.GRCm38.101.tx2gene.RDS")

# homo_sapiens.gtf <- import("/camp/home/ziffo/home/genomes/ensembl/Homo_sapiens.GRCh38.99.gtf") %>% as_tibble
# transcript.gtf <- import("/camp/home/ziffo/home/genomes/ensembl/Homo_sapiens.GRCh38.99.gtf") %>% as_tibble %>% distinct(transcript_id, .keep_all = TRUE) %>% mutate(coord = paste0("chr", seqnames, ":", start, "-", end))

# gene2ens <- homo_sapiens.gtf %>% dplyr::select(gene_id, gene_name)
gene2ens <- readRDS("/camp/home/ziffo/home/genomes/ensembl/Homo_sapiens.GRCh38.99.gene2ens.RDS") %>% unique
tx2gene <- readRDS("/camp/home/ziffo/home/genomes/ensembl/Homo_sapiens.GRCh38.99.chr_patch_hapl_scaff.tx2gene.rds") %>% rename("transcript_id" = TXNAME, "gene_id" = GENEID) %>% unique # data.frame with 2 columns: ENST ID & ENSG ID
# homo_sapiens.gtf %>% dplyr::select(transcript_id, gene_name) %>% drop_na %>% unique %>% saveRDS("/camp/home/ziffo/home/genomes/ensembl/Homo_sapiens.GRCh38.99.enst2gene.rds")
enst2gene <- readRDS("/camp/home/ziffo/home/genomes/ensembl/Homo_sapiens.GRCh38.99.enst2gene.rds") # data.frame with 2 columns: ENST ID & ENSG ID

ensembl_transcript_biotype <- readRDS("/camp/home/ziffo/home/genomes/ensembl/Homo_sapiens.GRCh38.99.transcript_biotype.rds")
ensembl_gene_biotype <- readRDS("/camp/home/ziffo/home/genomes/ensembl/Homo_sapiens.GRCh38.99.gene_biotype.rds")


ens_to_gene <- function(ens) {
  gene <- replace(ens, ens %in% as.character(gene2ens$gene_id), as.character(gene2ens$gene_name[gene2ens$gene_id %in% ens]))
  # gene <- gene2ens %>% dplyr::filter(gene_id == ens) %>% dplyr::pull(gene_name) %>% as.character()
  return(gene)
}
gene_to_ens <- function(gene) {
  ens <- gene2ens %>% dplyr::filter(gene_name == gene) %>% dplyr::pull(gene_id) %>% as.character()
  print(ens)
}
ens_to_gene_list <- function(ens) {
  gene <- list()
  for (i in seq_along(ens)) {
    gene[[i]] <- replace(ens[[i]], ens[[i]] %in% as.character(gene2ens$gene_id), as.character(gene2ens$gene_name[gene2ens$gene_id %in% ens[[i]]]))
  }
  return(gene)
}

# Sys.setenv("http_proxy" = "http://my.proxy.org:9999") # Configure biomaRt. NB this will stop gprofiler2 working 
# Sys.unsetenv("http_proxy") # unset proxy to fix gprofiler2 
# httr::set_config(httr::config(ssl_verifypeer = FALSE)) # httr::set_config(httr::config(ssl_cipher_list = "DEFAULT@SECLEVEL=1"))

# human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
# mouse2human = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = mus_musculus.gtf$gene_name, mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
mouse2human = readRDS("/camp/home/ziffo/home/genomes/ensembl/mouse2human.rds")

# # Convert mouse to human gene names https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/
# convertMouseGeneList <- function(x){
#   require("biomaRt")
#   human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#   mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
#   genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = mus_musculus.gtf$gene_name, mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
#   humanx <- unique(genesV2[, 2])
#   # Print the first 6 genes found to the screen
#   print(head(humanx))
#   return(humanx)
# }


sinaplot <- 
  function(data = ac_who_ir_dge, groups = "p0.05_reliable_IRdirection", continuous = "log2FoldChange", stat.y.position = 1, tip.length = 0.01, stats.test = "t_test", x.position = NULL, dashed.line = 0, width = 0.2, flip = "no", colours = 4){
    data2 <- data %>% dplyr::select(x = groups, y = continuous)
    print(paste("Stats test used = ", stats.test))
    if(stats.test == "t_test"){
      stats_test <- data2 %>% t_test(y ~ x, var.equal = TRUE) %>% adjust_pvalue() %>%  add_significance() %>% mutate(y.position = stat.y.position)
      print(stats_test)
    } else{
      stats_test <- data2 %>% wilcox_test(y ~ x) %>% adjust_pvalue() %>%  add_significance() %>% mutate(y.position = stat.y.position)
      print(stats_test)
    }
    if(flip == "yes"){
      if(dashed.line == "NA"){
        ggplot(data2, aes(x = x, y= y, color = x)) + 
          geom_violin() + geom_sina(size = 0.5) + geom_boxplot(data = data2, aes(fill = x), colour = "black", width=width, outlier.shape = NA) +
          scale_colour_manual(values = get_palette("npg",colours)) + scale_fill_manual(values = get_palette("npg",colours)) +  
          theme_classic() + theme(legend.position = "none", axis.text.x=element_text(size=10), axis.text.y=element_text(size=8), axis.title = element_text(size=10)) + coord_flip() +
          # geom_hline(yintercept = dashed.line, linetype = 2) + 
          stat_pvalue_manual(stats_test, label = "p.adj.signif", tip.length = tip.length, x = x.position, hide.ns = TRUE, size = 3)
      } else{
        ggplot(data2, aes(x = x, y= y, color = x)) + 
          geom_violin() + geom_sina(size = 0.5) + geom_boxplot(data = data2, aes(fill = x), colour = "black", width=width, outlier.shape = NA) +
          scale_colour_manual(values = get_palette("npg",colours)) + scale_fill_manual(values = get_palette("npg",colours)) +  
          theme_classic() + theme(legend.position = "none", axis.text.x=element_text(size=10), axis.text.y=element_text(size=8), axis.title = element_text(size=10)) + coord_flip() +
          geom_hline(yintercept = dashed.line, linetype = 2) + coord_flip() +
          stat_pvalue_manual(stats_test, label = "p.adj.signif", tip.length = tip.length, x = x.position, hide.ns = TRUE, size = 3) # add manually p values from binomial test
      }
    } else{
      if(dashed.line == "NA"){
        ggplot(data2, aes(x = x, y= y, color = x)) + 
          geom_violin() + geom_sina(size = 0.5) + geom_boxplot(data = data2, aes(fill = x), colour = "black", width=width, outlier.shape = NA) +
          # scale_colour_manual(values = get_palette("npg",colours)) + scale_fill_manual(values = get_palette("npg",colours)) +
          scale_colour_manual(values = c("dodgerblue2", "firebrick2")) + scale_fill_manual(values = c("dodgerblue2", "firebrick2")) +
          theme_classic() + theme(legend.position = "none", axis.text.x=element_text(size=10), axis.text.y=element_text(size=8), axis.title = element_text(size=10)) + #coord_flip() +
          # geom_hline(yintercept = dashed.line, linetype = 2) + 
          stat_pvalue_manual(stats_test, label = "p.adj.signif", tip.length = tip.length, x = x.position, hide.ns = TRUE, size = 3)
      } else{
        ggplot(data2, aes(x = x, y= y, color = x)) + 
          geom_violin() + geom_sina(size = 0.5) + geom_boxplot(data = data2, aes(fill = x), colour = "black", width=width, outlier.shape = NA) +
          # scale_colour_manual(values = get_palette("npg",colours)) + scale_fill_manual(values = get_palette("npg",colours)) +
          scale_colour_manual(values = c("dodgerblue2", "firebrick2")) + scale_fill_manual(values = c("dodgerblue2", "firebrick2")) +
          theme_classic() + theme(legend.position = "none", axis.text.x=element_text(size=10), axis.text.y=element_text(size=8), axis.title = element_text(size=10)) + #coord_flip() +
          geom_hline(yintercept = dashed.line, linetype = 2) + #coord_flip() +
          stat_pvalue_manual(stats_test, label = "p.adj.signif", tip.length = tip.length, x = x.position, hide.ns = TRUE, size = 3) # add manually p values from binomial test
      }
    }    
  }


# sinaplot(data = ac_who_ir_dge, groups = "p0.05_reliable_IRdirection", continuous = "log2FoldChange") + ylim(c(-1,1)) + ylab(label = expression(Gene~expression~log[2]~fold~change~VCP:CTRL) ) + xlab("")


sinaplot.one.sample <- 
  function(data = boxplot.lfc.merged, groups = "category", continuous = "log2FoldChange", flip = "yes", tip.length = 0.01, labels = "gene_name", colours = 3, stat.y.position = 1.00, label_number = 10, stats.test = "t_test", expected = 0, width = 0.2, dashed.line = 0){
    data2 <- data %>% dplyr::select(group = groups, lfc = continuous, label = labels)
    print(data2)
    print(paste("Stats test used = ", stats.test))
    if(stats.test == "t_test"){
      stats_test <- data2 %>% dplyr::select(group, lfc) %>% group_by(group) %>% filter(lfc != Inf) %>% filter(lfc != -Inf) %>% t_test(lfc ~ 1, mu = expected) %>% add_significance() %>% mutate(y.position = stat.y.position)
      print(stats_test)
    } else{
      stats_test <- data2 %>% dplyr::select(group, lfc) %>% group_by(group) %>% wilcox_test(lfc ~ 1, mu = expected) %>% add_significance() %>% mutate(y.position = stat.y.position)
      print(stats_test)
    }
    if(flip == "yes"){
      print(paste("Flipping coords"))
      if(colours < 5){
        ggplot(data2, aes(x = group, y= lfc)) + 
          geom_violin(aes(colour = group)) + 
          geom_sina(size = 0.5, aes(colour = group)) + 
          geom_boxplot(data = data2, aes(fill = group), colour = "black", width=width, outlier.shape = NA) +
          scale_colour_manual(values = c("firebrick2", "dodgerblue2", "forestgreen", "gold2")) + scale_fill_manual(values = c("firebrick2", "dodgerblue2", "forestgreen", "gold2")) +
          theme_classic() + theme(legend.position = "none", axis.text =element_text(size=10), axis.title = element_text(size=10)) +
          geom_hline(yintercept =dashed.line, linetype = 2) + coord_flip() +
          geom_text_repel(data = top_n(data2, label_number, abs(lfc)), aes(x = group, y = lfc, label = label), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0, size=2.7, max.overlaps = Inf) +
          stat_pvalue_manual(stats_test, label = "p.signif", tip.length = tip.length, xmin = "group", xmax = NULL, hide.ns = TRUE, size = 3) #, position_dodge(0.8)) # add manually p values
      } else{
        ggplot(data2, aes(x = group, y= lfc, color = group)) + 
          geom_violin(aes(colour = group)) + 
          geom_sina(size = 0.5, aes(colour = group)) + 
          geom_boxplot(data = data2, aes(fill = group), colour = "black", width=width, outlier.shape = NA) +
          scale_colour_manual(values = get_palette("npg",colours)) + scale_fill_manual(values = get_palette("npg",colours)) +
          theme_classic() + theme(legend.position = "none", axis.text =element_text(size=10), axis.title = element_text(size=10)) +
          geom_hline(yintercept =dashed.line, linetype = 2) + coord_flip() +
          geom_text_repel(data = top_n(data2, label_number, abs(lfc)), aes(x = group, y = lfc, label = label), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0, size=2.7, max.overlaps = Inf) +
          stat_pvalue_manual(stats_test, label = "p.signif", tip.length = tip.length, xmin = "group", xmax = NULL, hide.ns = TRUE, size = 3) #, position_dodge(0.8)) # add manually p values
      }
    } else{
      print(paste("Not flipping coords"))
      if(colours < 5){
        ggplot(data2, aes(x = group, y= lfc, color = group)) + 
          geom_violin() +  #aes(colour = group)
          geom_sina(size = 0.5) + # , aes(colour = group)
          geom_boxplot(data = data2, aes(fill = group), colour = "black", width=width, outlier.shape = NA) + 
          scale_colour_manual(values = c("firebrick2", "dodgerblue2", "forestgreen", "gold2")) + scale_fill_manual(values = c("firebrick2", "dodgerblue2", "forestgreen", "gold2")) +
          theme_classic() + theme(legend.position = "none", axis.text = element_text(size=10), axis.title = element_text(size=10)) +
          geom_hline(yintercept = dashed.line, linetype = 2) + #coord_flip() +
          geom_text_repel(data = top_n(data2, label_number, abs(lfc)), aes(x = group, y = lfc, label = label), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0, size=2.7, max.overlaps = Inf) +
          stat_pvalue_manual(stats_test, label = "p.signif", tip.length = tip.length, xmin = "group", xmax = NULL, hide.ns = TRUE, size = 3) #, position_dodge(0.8)) # add manually p values
      } else{
        ggplot(data2, aes(x = group, y= lfc, color = group)) + 
          geom_violin(aes(colour = group)) + 
          geom_sina(size = 0.5, aes(colour = group)) + 
          geom_boxplot(data = data2, aes(fill = group), colour = "black", width=width, outlier.shape = NA) +
          scale_colour_manual(values = get_palette("npg",colours)) + scale_fill_manual(values = get_palette("npg",colours)) +
          theme_classic() + theme(legend.position = "none", axis.text = element_text(size=10), axis.title = element_text(size=10)) +
          geom_hline(yintercept = dashed.line, linetype = 2) + #coord_flip() +
          geom_text_repel(data = top_n(data2, label_number, abs(lfc)), aes(x = group, y = lfc, label = label), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0, size=2.7, max.overlaps = Inf) +
          stat_pvalue_manual(stats_test, label = "p.signif", tip.length = tip.length, xmin = "group", xmax = NULL, hide.ns = TRUE, size = 3) #, position_dodge(0.8)) # add manually p values
      }
    }
  }



# sinaplot.one.sample(data = boxplot.lfc.merged, groups = "category", continuous = "log2FoldChange") + ylim(c(-1,1)) + ylab(label = expression(Gene~expression~log[2]~fold~change~VCP:CTRL) ) + xlab("")

display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

# DESeq2 Functions ------------------
DESeq.analysis <- function(metadata = treated.who.metadata, design = ~ 1, contrast = "NA", species = "human", transcript.level = TRUE){
  metadata.files <- metadata$file_salmon
  names(metadata.files) <- metadata$sample
  rownames(metadata) <- metadata$sample
  if(length(contrast) == 1){
    cat(red(paste(contrast, "result contrast specified")))
    if(species == "human"){
      cat(blue("\nusing human gene names\n"))
      dds <- tximport(metadata.files, type="salmon", tx2gene=tx2gene) %>% DESeqDataSetFromTximport(colData = metadata, design = design) %>% DESeq()
      print(resultsNames(dds))
      vsd <- vst(dds, blind=FALSE)
      vsd.counts <- as_tibble(assay(vsd), rownames = "gene_id")  %>% left_join(gene2ens)
      if(contrast != "NA"){ res <- DESeq2::results(dds, name = contrast) %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(pvalue) }
      if(transcript.level == TRUE){
        cat(blue("transcript level analysis\n"))
        dds.tx <- tximport(metadata.files, type="salmon", txOut = TRUE) %>% DESeqDataSetFromTximport(colData = metadata, design = design) %>% DESeq()
        if(contrast != "NA"){ res.tx <- DESeq2::results(dds.tx, name = contrast) %>% as_tibble(rownames = "transcript_id") %>% left_join(tx2gene, by="transcript_id") %>% arrange(pvalue) }
        }
      } else if(species == "mouse"){
      cat(blue("\nusing mouse gene names\n"))
      dds <- tximport(metadata.files, type="salmon", tx2gene=mus.musculus.tx2gene) %>% DESeqDataSetFromTximport(colData = metadata, design = design) %>% DESeq()
      print(resultsNames(dds))
      vsd <- vst(dds, blind=FALSE)
      vsd.counts <- as_tibble(assay(vsd), rownames = "gene_id")  %>% left_join(mus.musculus.gene2ens)
      if(contrast != "NA"){ res <- DESeq2::results(dds, name = contrast) %>% as_tibble(rownames = "gene_id") %>% left_join(mus.musculus.gene2ens) %>% mutate(gene_name.mouse = gene_name, gene_name = toupper(gene_name)) %>% arrange(pvalue) }
      if(transcript.level == TRUE){
        cat(blue("transcript level analysis\n"))
        dds.tx <- tximport(metadata.files, type="salmon", txOut = TRUE) %>% DESeqDataSetFromTximport(colData = metadata, design = design) %>% DESeq()
        if(contrast != "NA"){ res.tx <- DESeq2::results(dds.tx, name = contrast) %>% as_tibble(rownames = "transcript_id") %>% left_join(mus.musculus.tx2gene, by="transcript_id") %>% mutate(gene_name.mouse = gene_name, gene_name = toupper(gene_name)) %>% arrange(pvalue) }
        }
      }
    if(contrast != "NA"){ cat(blue(paste("Significant events: padj < 0.05 = ", nrow(filter(res, padj < 0.05)), ", pvalue < 0.05 = ", nrow(filter(res, pvalue < 0.05))))) }
  } else if(length(contrast) > 1){
    cat(green(paste(length(contrast), "contrasts specified")))
    if(species == "human"){
      cat(blue("\nusing human gene names\n"))
      dds <- tximport(metadata.files, type="salmon", tx2gene=tx2gene) %>% DESeqDataSetFromTximport(colData = metadata, design = design) %>% DESeq()
      print(resultsNames(dds))
      res <- list()
      for(i in seq_along(contrast)) {
        res[[i]] <- DESeq2::results(dds, name = contrast[[i]]) %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(pvalue) 
        names(res)[i] <- contrast[[i]]
      }
      vsd <- vst(dds, blind=FALSE)
      vsd.counts <- as_tibble(assay(vsd), rownames = "gene_id")  %>% left_join(gene2ens)
      if(transcript.level == TRUE){
        cat(blue("\ntranscript level analysis\n"))
        dds.tx <- tximport(metadata.files, type="salmon", txOut = TRUE) %>% DESeqDataSetFromTximport(colData = metadata, design = design) %>% DESeq()
        res.tx <- list()
        for(i in seq_along(contrast)) {
          res.tx[[i]] <- DESeq2::results(dds.tx, name = contrast[[i]]) %>% as_tibble(rownames = "transcript_id") %>% left_join(tx2gene, by="transcript_id") %>% arrange(pvalue) 
          names(res.tx)[i] <- contrast[[i]]
        }
      }
      } else if(species == "mouse") {
        cat(blue("\nusing mouse gene names\n"))
        dds <- tximport(metadata.files, type="salmon", tx2gene=mus.musculus.tx2gene) %>% DESeqDataSetFromTximport(colData = metadata, design = design) %>% DESeq()
        print(resultsNames(dds))
        res <- list()
        for(i in seq_along(contrast)) {
          res[[i]] <- DESeq2::results(dds, name = contrast[[i]]) %>% as_tibble(rownames = "gene_id") %>% left_join(mus.musculus.gene2ens) %>% mutate(gene_name.mouse = gene_name, gene_name = toupper(gene_name)) %>% arrange(pvalue)
          names(res)[i] <- contrast[[i]]
        }
        vsd <- vst(dds, blind=FALSE)
        vsd.counts <- as_tibble(assay(vsd), rownames = "gene_id")  %>% left_join(mus.musculus.gene2ens)
        if(transcript.level == TRUE){
          cat(blue("\ntranscript level analysis\n"))
          dds.tx <- tximport(metadata.files, type="salmon", txOut = TRUE) %>% DESeqDataSetFromTximport(colData = metadata, design = design) %>% DESeq()
          res.tx <- list()
          for(i in seq_along(contrast)) {
            res.tx[[i]] <- DESeq2::results(dds.tx, name = contrast[[i]]) %>% as_tibble(rownames = "transcript_id") %>% left_join(mus.musculus.tx2gene, by="transcript_id") %>% mutate(gene_name.mouse = gene_name, gene_name = toupper(gene_name)) %>% arrange(pvalue)
            names(res.tx)[i] <- contrast[[i]]
          }
        }
      }
  }
  if(transcript.level == TRUE & contrast != "NA"){ return(list(dds=dds,vsd=vsd,vsd.counts=vsd.counts,res=res,res.tx=res.tx)) } 
  else if(transcript.level == FALSE & contrast != "NA"){ return(list(dds=dds,vsd=vsd,vsd.counts=vsd.counts,res=res)) }
  else if(contrast == "NA"){ return(list(dds=dds,vsd=vsd,vsd.counts=vsd.counts)) }
}




curated_profiles <- 
  function(gprofiles = IRevents_gprofiles_curated, colours = 4, cols =  c("firebrick2", "dodgerblue2", "forestgreen", "gold2")){
    gprofiles <- gprofiles %>% distinct(term_name, .keep_all = TRUE) %>%  arrange( -log10(p_value) ) %>%  mutate(term_name = factor(term_name, levels = term_name)) #%>%
    # filter(term_name %in% curated_terms) %>% separate(col = query, into = c("group", "direction"),sep ="_", remove = FALSE) %>% filter( group == group_name )  %>% 
    if(colours < 5){
      ggplot(gprofiles, aes( x = term_name, y = -log10(p_value), fill = query) ) + 
        geom_col(aes(y = ( (-log10(p_value))/sum(-log10(p_value)) * 100.0) )) + # makes more space to axis labels
        scale_fill_manual(values = cols) +
        facet_grid(query ~ ., space = "free", scales = "free") +
        coord_flip() + theme_classic() + guides(fill=FALSE) +
        ylab(label = expression(log[10]~(P~value)) ) + xlab("") +
        theme(#strip.background = element_blank(), 
              strip.text.y = element_text(size=8, face = "bold", angle = 90), axis.text= element_text(size = 8), axis.title = element_text(size = 9), legend.position = "none", 
              plot.subtitle = element_text(size = 10)) + # legend.position = c(0.8, 0.2), 
        scale_y_continuous(expand = c(0, 0))
    } else{
      ggplot(gprofiles, aes( x = term_name, y = -log10(p_value), fill = query) ) + 
        geom_col(aes(y = ( (-log10(p_value))/sum(-log10(p_value)) * 100.0) )) + # makes more space to axis labels
        scale_fill_manual(values = get_palette("npg", 10)) +
        facet_grid(query ~ ., space = "free", scales = "free") +
        coord_flip() + theme_classic() + guides(fill=FALSE) +
        ylab(label = expression(log[10]~(P~value)) ) + xlab("") +
        theme(strip.background = element_blank(), strip.text.y = element_text(size=8, face = "bold",angle = 90), axis.text= element_text(size = 8), axis.title = element_text(size = 9), legend.position = "none", 
              plot.subtitle = element_text(size = 10)) + # legend.position = c(0.8, 0.2), 
        scale_y_continuous(expand = c(0, 0))
    }
  }


signed_z <- function(log2fc, pvalue){
  z <- ifelse(test = log2fc > 0,
              yes = qnorm(1 - (pvalue / 2)),
              no = qnorm(pvalue / 2) )
  return(z)
}


nest_batch <- function(inner, within, set_to=".") {
  n_instances <- apply(table(inner, within)!=0, 1, sum)
  if (any(n_instances>1)) {
    stop(paste(names(n_instances)[n_instances>1], collapse=", "), " appear in multiple parents")
  }
  levels_in <- levels(inner)
  # each inner level has in a unique group, find it
  corresponding_group <-within[match(levels_in, inner)]
  # for each group, find index of first inner level
  which_inner <- match(levels(within), corresponding_group)
  # and set it to the common value
  if (any(levels(inner)[-which_inner]==set_to )) { # we're going to duplicate an existing level so best stop
    stop("One of the batches is already called ", set_to, ". Please use a non-existing level")
  }
  levels(inner)[which_inner] <- set_to
  inner
}

res_top <- function(r, i=1) {
  ind <- which(r$padj<0.05) # indexes of significant genes
  ind[order(abs(r$log2FoldChange)[ind])][i]
}

ma_plot_deseq <- 
  function(data = mn_d25.res, labels_list = ma.labels){
    labels_maplot <- data %>% filter(gene_name %in% labels_list)
    ggplot(data = data, aes(x = baseMean, y = log2FoldChange, colour = padj < 0.05)) +  # NB if pvalue == NA then will not show
      geom_point(size = 0.5) +
      scale_colour_manual( values = c("darkgray", "firebrick2") ) +  theme_bw() +  theme(panel.grid = element_blank(), legend.title = element_blank(), legend.position = "none") +
      xlab(expression( Log[2]~mean~expression )) + ylab(expression(Log[2]~fold~change~VCP:CTRL)) + 
      geom_text_repel(data = labels_maplot, aes(x = baseMean, y = log2FoldChange, label = gene_name), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0, size=3, max.overlaps = Inf) +
      geom_hline(yintercept = 0, linetype = 2)
  }

ma_plot_deseq.pvalue <- 
  function(data = mn_d25.res, labels_list = ma.labels){
    labels_maplot <- data %>% filter(gene_name %in% labels_list) %>% arrange(pvalue) %>% distinct(gene_name, .keep_all = TRUE)
    ggplot(data = data, aes(x = baseMean, y = log2FoldChange, colour = pvalue < 0.05)) +  # NB if pvalue == NA then will not show
      geom_point(size = 0.5) +
      scale_colour_manual( values = c("darkgray", "firebrick2") ) +  theme_bw() +  theme(panel.grid = element_blank(), legend.title = element_blank(), legend.position = "none") +
      xlab(expression( Log[2]~mean~expression )) + ylab(expression(Log[2]~fold~change~VCP:CTRL)) + 
      geom_text_repel(data = labels_maplot, aes(x = baseMean, y = log2FoldChange, label = gene_name), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0, size=3, max.overlaps = Inf) +
      geom_hline(yintercept = 0, linetype = 2)
  }


volcano_plot_deseq <- 
  function(data = mn_d25.res, labels_list = volcano.labels){
    labels_plot <- data %>% filter(gene_name %in% labels_list, padj < 0.05)
    data2 <- data %>% mutate(direction = case_when(padj < 0.05 & log2FoldChange > 0 ~ "upregulated", padj < 0.05 & log2FoldChange < 0 ~ "downregulated", TRUE ~ "unchanged"))
    ggplot(data = data2, aes(x = log2FoldChange, y =  -log10(padj), colour = direction)) +  geom_point(size = 0.5) + #xlim(1,5) + ylim(-3,3) +
      scale_colour_manual( values = c(unchanged = "darkgray", upregulated = "firebrick2", downregulated = "dodgerblue2") ) +  theme_bw() +  theme(panel.grid = element_blank(), legend.title = element_blank(), legend.position = "none") + 
      ylab(expression( -log[10]~FDR )) + 
      geom_text_repel(data = labels_plot, aes(x = log2FoldChange, y = -log10(padj), label = gene_name), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0.5, size=2.6, max.overlaps = Inf) +
      geom_hline(yintercept = -log10(0.05), linetype = 2) + geom_vline(xintercept =0, linetype = 2)
  }

# volcano plotting undajusted p value
volcano_plot_deseq.pvalue <-  
  function(data = mn_d25.res, labels_list = volcano.labels){
    labels_plot <- data %>% filter(gene_name %in% labels_list, pvalue < 0.05) %>% arrange(pvalue) %>% distinct(gene_name, .keep_all = TRUE)
    data2 <- data %>% mutate(direction = case_when(pvalue < 0.05 & log2FoldChange > 0 ~ "upregulated", pvalue < 0.05 & log2FoldChange < 0 ~ "downregulated", TRUE ~ "unchanged"))
    ggplot(data = data2, aes(x = log2FoldChange, y =  -log10(pvalue), colour = direction)) +  geom_point(size = 0.5) + #xlim(1,5) + ylim(-3,3) +
      scale_colour_manual( values = c(unchanged = "darkgray", upregulated = "firebrick2", downregulated = "dodgerblue2") ) +  
      theme_bw() +  theme(panel.grid = element_blank(), legend.title = element_blank(), legend.position = "none", axis.text=element_text(size=8), axis.title=element_text(size=10)) + 
      ylab(expression( -log[10]~P~value )) + 
      geom_text_repel(data = labels_plot, aes(x = log2FoldChange, y = -log10(pvalue), label = gene_name), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0.5, size=2.6, max.overlaps = Inf) +
      geom_hline(yintercept = -log10(0.05), linetype = 2) + geom_vline(xintercept =0, linetype = 2)
  }



# IR Functions ----------------------
IRFinder.analysis <- function(metadata = treated.who.metadata, sample.names = "sample", variable.name = "stimulated", ctrl = NULL, mut = NULL, batch = "NA", file.var = "Sample limsid", ge.res = treated.who.a1_vs_a0$res, animal = "human", 
                              design = ~ Condition + Condition:IRFinder, contrast = NULL,
                              irfinder.dir = "/camp/lab/luscomben/home/shared/projects/patani-collab/astrocyte-stimulated-bulk-rnaseq/splicing/irfinder"){
  if(length(variable.name) == 1){
  if(batch == "NA"){
    experiment = metadata %>% select(SampleNames = sample.names, Condition = variable.name, file.name = file.var) %>% 
      mutate(Condition = factor(Condition,levels=c(ctrl, mut)), file_paths.dir = file.path(irfinder.dir,file.name,"IRFinder-IR-dir.txt"), file_paths.nondir = file.path(irfinder.dir,file.name,"IRFinder-IR-nondir.txt"))
    if(FALSE %in% file.exists(experiment$file_paths.dir) == FALSE){
      cat(blue("Importing stranded IRFinder output: IRFinder-IR-dir.txt"))
      metaList = DESeqDataSetFromIRFinder(filePaths=experiment$file_paths.dir, designMatrix=experiment, designFormula = design)
    }else{
      cat(blue("Importing non-stranded IRFinder output: IRFinder-IR-nondir.txt"))
      metaList = DESeqDataSetFromIRFinder(filePaths=experiment$file_paths.nondir, designMatrix=experiment, designFormula = design)
    }
  } else{
    experiment = metadata %>% select(SampleNames = sample.names, Condition = variable.name, file.name = file.var, Batch = batch) %>% 
      mutate(Condition = factor(Condition,levels=c(ctrl, mut)), file_paths.dir = file.path(irfinder.dir,file.name,"IRFinder-IR-dir.txt"), file_paths.nondir = file.path(irfinder.dir,file.name,"IRFinder-IR-nondir.txt"))
    cat(blue(paste0("Accounting for batch: ", batch)))
    if(FALSE %in% file.exists(experiment$file_paths.dir) == FALSE){
      cat(blue("Importing stranded IRFinder output: IRFinder-IR-dir.txt"))
      metaList = DESeqDataSetFromIRFinder(filePaths=experiment$file_paths.dir, designMatrix=experiment, designFormula = ~ Batch + Condition + Condition:IRFinder)
    }else{
      cat(blue("Importing non-stranded IRFinder output: IRFinder-IR-nondir.txt"))
      metaList = DESeqDataSetFromIRFinder(filePaths=experiment$file_paths.nondir, designMatrix=experiment, designFormula = ~ Batch + Condition + Condition:IRFinder)
    }
  }
  dds = DESeq(metaList$DESeq2Object)
  res.ctrl = DESeq2::results(dds, name = paste0("Condition",ctrl,".IRFinderIR")) # tests IR reads vs spliced reads in ctrl
  res.mut = DESeq2::results(dds, name = paste0("Condition",mut,".IRFinderIR")) # get IR ratio in the mut samples
  res = DESeq2::results(dds, contrast=list(paste0("Condition",mut,".IRFinderIR"),  paste0("Condition",ctrl,".IRFinderIR")))
  ir.res = cleanIRFinder(mutant_vs_ctrl = res, mutant = res.mut, ctrl = res.ctrl, species = animal)
  } else if(length(variable.name) == 2){
    cat(red(paste("Multi-factor design examining ", variable.name[1], "and ", variable.name[2])))
    if(batch == "NA"){
      experiment = metadata %>% select(SampleNames = sample.names, Condition1 = variable.name[1], Condition2 = variable.name[2], file.name = file.var) %>% 
        mutate(Condition1 = as.factor(Condition1), Condition2 = as.factor(Condition2), # NB this is not explicit and will factor based on alphabetical order i.e. =c("cytoplasmic", "nuclear") & c("ctrl", "vcp")
               file_paths.dir = file.path(irfinder.dir,file.name,"IRFinder-IR-dir.txt"), file_paths.nondir = file.path(irfinder.dir,file.name,"IRFinder-IR-nondir.txt"))
      if(FALSE %in% file.exists(experiment$file_paths.dir) == FALSE){
        cat(blue("Importing stranded IRFinder output: IRFinder-IR-dir.txt"))
        metaList = DESeqDataSetFromIRFinder(filePaths=experiment$file_paths.dir, designMatrix=experiment, designFormula = ~ Condition1 * Condition2 * IRFinder)
      }else{
        cat(blue("Importing non-stranded IRFinder output: IRFinder-IR-nondir.txt"))
        metaList = DESeqDataSetFromIRFinder(filePaths=experiment$file_paths.nondir, designMatrix=experiment, designFormula = ~ Condition1 * Condition2 * IRFinder)
      }
    } else{
      experiment = metadata %>% select(SampleNames = sample.names, Condition1 = variable.name[1], Condition2 = variable.name[2], file.name = file.var, Batch = batch) %>% 
        mutate(Condition1 = as.factor(Condition1), Condition2 = as.factor(Condition2), # NB this is not explicit and will factor based on alphabetical order i.e. =c("cytoplasmic", "nuclear") & c("ctrl", "vcp")
               file_paths.dir = file.path(irfinder.dir,file.name,"IRFinder-IR-dir.txt"), file_paths.nondir = file.path(irfinder.dir,file.name,"IRFinder-IR-nondir.txt"))
      cat(blue(paste0("Accounting for batch: ", batch)))
      if(FALSE %in% file.exists(experiment$file_paths.dir) == FALSE){
        cat(blue("Importing stranded IRFinder output: IRFinder-IR-dir.txt"))
        metaList = DESeqDataSetFromIRFinder(filePaths=experiment$file_paths.dir, designMatrix=experiment, designFormula = ~ Batch + Condition1 * Condition2 * IRFinder)
      }else{
        cat(blue("Importing non-stranded IRFinder output: IRFinder-IR-nondir.txt"))
        metaList = DESeqDataSetFromIRFinder(filePaths=experiment$file_paths.nondir, designMatrix=experiment, designFormula = ~ Batch + Condition1 * Condition2 * IRFinder)
      }
    }
    dds = DESeq(metaList$DESeq2Object)
    resultsNames(dds)
    res = DESeq2::results(dds, name = contrast)
    ir.res = cleanIRFinder.multidesign(IRFinder = res, species = animal)
  }
  ir.ge.res = summariseIRjoinGE(cleanIRFinder = ir.res, GE = ge.res)
  return(list(irfinder=ir.res,ir.ge=ir.ge.res))
}


# using output from IRFinder DESeq2 GLM with replicates https://github.com/williamritchie/IRFinder/wiki/Generalized-Linear-Model-with-Replicates
DESeqDataSetFromIRFinder = function(filePaths,designMatrix,designFormula){
  irfinder.tsv = filePaths %>% map(read_tsv)
  names(irfinder.tsv) = designMatrix$SampleNames
  irtab <- map_dfr(irfinder.tsv, bind_rows, .id = "SampleNames") %>% mutate(IntronDepth = round(IntronDepth), SpliceExact = round(SpliceExact), MaxSplice = round(pmax(SpliceLeft, SpliceRight)), irnames = paste0(Name,"/",Chr,":",Start,"-",End,":",Strand))
  IntronDepth = irtab %>% pivot_wider(irnames, names_from = SampleNames, values_from = IntronDepth, names_glue = "intronDepth.{SampleNames}") %>% column_to_rownames("irnames")
  SpliceExact = irtab %>% pivot_wider(irnames, names_from = SampleNames, values_from = SpliceExact, names_glue = "totalSplice.{SampleNames}") %>% column_to_rownames("irnames")
  MaxSplice = irtab %>% pivot_wider(irnames, names_from = SampleNames, values_from = MaxSplice, names_glue = "maxSplice.{SampleNames}") %>% column_to_rownames("irnames")
  group = bind_rows(mutate(designMatrix, IRFinder = "IR"), mutate(designMatrix, IRFinder = "Splice")) %>% mutate(IRFinder = factor(IRFinder, levels=c("Splice","IR")))
  counts.IRFinder = bind_cols(IntronDepth,MaxSplice) %>% drop_na()
  
  dd = DESeqDataSetFromMatrix(countData = counts.IRFinder, colData = group, design = designFormula)
  sizeFactors(dd)=rep(1,nrow(group))
  final=list(dd,IntronDepth,SpliceExact,MaxSplice)
  names(final)=c("DESeq2Object","IntronDepth","SpliceDepth","MaxSplice")
  return(final)
}


cleanIRFinder <- function(mutant_vs_ctrl = irfinder.res.ac_nuc_vcp_vs_ctrl, mutant = irfinder.res_ac_nuc.vcp, ctrl = irfinder.res_ac_nuc.ctrl, deseq_res = NULL, mass_spec = NULL, species = "human"){
  irfinder.mutant <- mutant %>% as_tibble(rownames = "Intron.GeneName.GeneID.Coords") %>% 
    mutate(IR_vs_Splice.mutant = 2^.$log2FoldChange, IRratio.mutant = IR_vs_Splice.mutant/(1+IR_vs_Splice.mutant), baseMean.mutant = baseMean) %>% dplyr::select(Intron.GeneName.GeneID.Coords, IRratio.mutant, baseMean.mutant) 
  irfinder.ctrl <- ctrl %>% as_tibble(rownames = "Intron.GeneName.GeneID.Coords") %>% 
    mutate(IR_vs_Splice.ctrl = 2^.$log2FoldChange, IRratio.ctrl = IR_vs_Splice.ctrl/(1+IR_vs_Splice.ctrl), baseMean.ctrl = baseMean) %>% dplyr::select(Intron.GeneName.GeneID.Coords, IRratio.ctrl, baseMean.ctrl) 
  IRFinder <- mutant_vs_ctrl %>% as_tibble(rownames = "Intron.GeneName.GeneID.Coords") %>% 
    left_join(irfinder.mutant, by = "Intron.GeneName.GeneID.Coords") %>% left_join(irfinder.ctrl, by = "Intron.GeneName.GeneID.Coords") %>% 
    mutate(gene_name = str_split_fixed(.$Intron.GeneName.GeneID.Coords, "/", 4)[,1], gene_id = str_split_fixed(.$Intron.GeneName.GeneID.Coords, "/", 4)[,2], 
           
           IRratio.diff = IRratio.mutant - IRratio.ctrl, 
           IRdirection = case_when(IRratio.diff > 0 ~ "up", IRratio.diff < 0 ~ "down", TRUE ~ "unchanged"), 
           
           reliable = case_when(baseMean < 10 ~ "unreliable", TRUE ~ "reliable"), 
           retained = case_when((IRratio.mutant >= 0.1 | IRratio.ctrl >= 0.1) ~ "retained", TRUE ~ "spliced"),
           reliable_retained = case_when(reliable == "reliable" & retained == "retained" ~ "retained", TRUE ~ "spliced"),
           
           p0.05 = case_when(pvalue < 0.05 ~ "significant", TRUE ~ "none_significant"), p0.05 = factor(p0.05, levels = c("significant", "none_significant")),
           p0.05_IRdirection = case_when(IRdirection == "up" & pvalue < 0.05 ~ "up", IRdirection == "down" & pvalue < 0.05 ~ "down", TRUE ~ "none_significant"),
           p0.05_reliable = case_when(pvalue < 0.05 & reliable == "reliable" ~ "significant", TRUE ~ "none_significant"),
           p0.05_reliable_IRdirection = case_when(IRdirection == "up" & p0.05_reliable == "significant" ~ "IR up", IRdirection == "down" & p0.05_reliable  == "significant" ~ "IR down", TRUE ~ "none_significant"),
           
           padj0.05 = case_when(padj < 0.05 ~ "significant", TRUE ~ "none_significant"), padj0.05 = factor(padj0.05, levels = c("significant", "none_significant")),
           padj0.05_IRdirection = case_when(IRdirection == "up" & padj < 0.05 ~ "up", IRdirection == "down" & padj < 0.05 ~ "down", TRUE ~ "none_significant"),
           padj0.05_reliable = case_when(padj < 0.05 & reliable == "reliable" ~ "significant", TRUE ~ "none_significant"),
           padj0.05_reliable_IRdirection = case_when(IRdirection == "up" & padj0.05_reliable == "significant" ~ "IR up", IRdirection == "down" & padj0.05_reliable  == "significant" ~ "IR down", TRUE ~ "none_significant"),
           
           intron_type = str_split_fixed(.$Intron.GeneName.GeneID.Coords, "/", 4)[,3], intron_length = GenomicRanges::width(as.GenomicRange(Intron.GeneName.GeneID.Coords)), 
           coords = str_split_fixed(.$Intron.GeneName.GeneID.Coords, "/",4)[,4], Chr = str_split_fixed(.$Intron.GeneName.GeneID.Coords, "/|:", 6)[,4], Start = as.numeric(str_split_fixed(.$Intron.GeneName.GeneID.Coords, "/|:|-", 7)[,5]), 
           End = as.numeric(str_split_fixed(.$Intron.GeneName.GeneID.Coords, "/|:|-", 7)[,6]), Direction = str_split_fixed(.$Intron.GeneName.GeneID.Coords, ":", 3)[,3])
  print(paste(nrow(filter(IRFinder, pvalue < 0.05)), "significant IR events P < 0.05"))
  print(table(IRFinder$p0.05_IRdirection))
  print(paste(nrow(filter(IRFinder, padj < 0.05)), "significant IR events Padj < 0.05"))
  print(table(IRFinder$padj0.05_IRdirection))
  gr_introns <- as.GenomicRange(IRFinder$Intron.GeneName.GeneID.Coords)
  seqlevelsStyle(gr_introns) <- "UCSC"
  if(species == "human"){
    print("GC content based on human")
    IRFinder$gc_content <- GCcontent(Hsapiens, gr_introns)
  }else if(species == "mouse"){
    print("GC content based on mouse")
    IRFinder$gc_content <- GCcontent(Mmusculus, gr_introns)
  }else if(species == "rat"){
    print("GC content based on rat")
    IRFinder$gc_content <- GCcontent(Rnorvegicus, gr_introns)
  }
  if(species != "human"){
    IRFinder <- IRFinder %>% mutate(gene_name = toupper(gene_name))  # covert gene_names to upper case for matching to Deseq results & Human gene names
  }
  # IRFinder$phast_cons <- gscores(phast, gr_introns) # too slow to keep in function
  # if(!is.null(deseq_res)){
  #   deseq_res <- deseq_res %>% rename(ge.log2FoldChange = log2FoldChange, ge.padj = padj, ge.pvalue = pvalue, ge.baseMean = baseMean) %>% dplyr::select(-lfcSE, -stat)
  #   IRFinder <- IRFinder %>% left_join(deseq_res, by = c("gene_id", "gene_name")) %>% dplyr::mutate(ge.direction = case_when(ge.log2FoldChange > 0 ~ "up", TRUE ~ "down"),
  #                                          direction = case_when(ge.log2FoldChange > 0 & IRratio.diff < 0 ~ "IR down & expression up",
  #                                                                ge.log2FoldChange > 0 & IRratio.diff > 0 ~ "IR up & expression up",
  #                                                                ge.log2FoldChange < 0 & IRratio.diff < 0 ~ "IR down & expression down",
  #                                                                ge.log2FoldChange < 0 & IRratio.diff > 0 ~ "IR up & expression down"))
  # }
  # if(!is.null(mass_spec)){
  #   IRFinder <- IRFinder %>%  left_join(mass_spec, by = c("gene_name"="name"))
  #   IRFinder <- IRFinder %>% dplyr::mutate(massspec_direction = case_when(vcp_vs_ctrl_ratio > 0 ~ "up", vcp_vs_ctrl_ratio < 0 ~ "down"),
  #                                          massspec_ir_direction = case_when(vcp_vs_ctrl_ratio > 0 & IRratio.diff < 0 ~ "IR down & protein up", vcp_vs_ctrl_ratio > 0 & IRratio.diff > 0 ~ "IR up & protein up",
  #                                                                            vcp_vs_ctrl_ratio < 0 & IRratio.diff < 0 ~ "IR down & protein down", vcp_vs_ctrl_ratio < 0 & IRratio.diff > 0 ~ "IR up & protein down"),
  #                                          massspec_expression_ir_direction = case_when(direction == "IR down & expression up" & massspec_direction == "up" ~ "IR down & expression up & protein up",
  #                                                                                       direction == "IR down & expression up" & massspec_direction == "down" ~ "IR down & expression up & protein down", 
  #                                                                                       direction == "IR down & expression down" & massspec_direction == "up" ~ "IR down & expression up & protein down", 
  #                                                                                       direction == "IR down & expression down" & massspec_direction == "down" ~ "IR down & expression down & protein down",
  #                                                                                       direction == "IR up & expression down" & massspec_direction == "up" ~ "IR up & expression down & protein up",
  #                                                                                       direction == "IR up & expression up" & massspec_direction == "up" ~ "IR up & expression up & protein up", 
  #                                                                                       direction == "IR up & expression up" & massspec_direction == "down" ~ "IR up & expression up & protein down",
  #                                                                                       direction == "IR up & expression down" & massspec_direction == "down" ~ "IR up & expression down & protein up",
  #                                                                                       TRUE ~ "other"))
  # }
  return(IRFinder)
}


# normalise IR per gene and left_join to GE results
summariseIRjoinGE <- function(cleanIRFinder = irfinder.ac_nuc_vcp_vs_ctrl, GE = ac_cyt_vcp_vs_ctrl.res){
  cleanIRFinder.summarised <- cleanIRFinder %>% group_by(gene_id, gene_name) %>% 
    filter(reliable == "reliable") %>% 
    summarise(all_intron_lengths = sum(intron_length), introns = dplyr::n(), 
              IR.log2FoldChange.norm = sum( log2FoldChange * (intron_length / all_intron_lengths)), 
              # IRratio.diff.norm = sum( IRratio.diff * (intron_length / all_intron_lengths)), 
              IR.pvalue = min(pvalue, na.rm = TRUE), 
              IR.baseMean = sum( baseMean * (intron_length / all_intron_lengths))) %>% ungroup %>% 
    mutate(IR.direction = case_when(IR.log2FoldChange.norm > 0 ~ "IR up", IR.log2FoldChange.norm <= 0 ~ "IR down")) %>% 
    dplyr::select(gene_id, gene_name, IR.log2FoldChange.norm, IR.pvalue, IR.baseMean, IR.direction) # , IRratio.diff.norm
  cleanIRFinder.summarised.GE <- GE %>% left_join(cleanIRFinder.summarised, by = c("gene_id", "gene_name")) %>% 
    mutate(GE.direction = case_when(log2FoldChange > 0 ~ "GE up", log2FoldChange < 0 ~ "GE down")) %>% 
    rename(GE.log2FoldChange = log2FoldChange, GE.baseMean = baseMean, GE.padj = padj, GE.pvalue = pvalue) %>%
    mutate(significant = case_when(GE.pvalue < 0.05 & IR.pvalue < 0.05 ~ "both", GE.pvalue < 0.05 ~ "GE only", IR.pvalue < 0.05 ~ "IR only", TRUE ~ "None"), 
           direction = case_when(IR.direction == "IR up" & GE.direction == "GE down" ~ "IR up mRNA down", IR.direction == "IR up" & GE.direction == "GE up" ~ "IR up mRNA up",
                                 IR.direction == "IR down" & GE.direction == "GE down" ~ "IR down mRNA down", IR.direction == "IR down" & GE.direction == "GE up" ~ "IR down mRNA up")) %>%
    dplyr::select(gene_id, gene_name, GE.log2FoldChange, GE.padj, GE.pvalue, GE.baseMean, IR.log2FoldChange.norm, IR.pvalue, IR.baseMean, direction, significant, IR.direction, GE.direction) # , IRratio.diff.norm
  print(paste("Normalising", nrow(cleanIRFinder), "IR events in", nrow(GE), "GE events"))
  return(cleanIRFinder.summarised.GE)
}




cleanIRFinder.multidesign <- function(IRFinder = irfinder.ac.nuc_vs_cyt.vcp_vs_ctrl, species = "human"){
  IRFinder <- IRFinder %>% as_tibble(rownames = "Intron.GeneName.GeneID.Coords") %>%
    mutate(gene_name = str_split_fixed(.$Intron.GeneName.GeneID.Coords, "/", 4)[,1], gene_id = str_split_fixed(.$Intron.GeneName.GeneID.Coords, "/", 4)[,2], 
           IRdirection = case_when(log2FoldChange > 0 ~ "up", log2FoldChange < 0 ~ "down", TRUE ~ "unchanged"), 
           reliable = case_when(baseMean < 10 ~ "unreliable", TRUE ~ "reliable"), 
           p0.05 = case_when(pvalue < 0.05 ~ "significant", TRUE ~ "none_significant"), p0.05 = factor(p0.05, levels = c("significant", "none_significant")),
           p0.05_IRdirection = case_when(IRdirection == "up" & pvalue < 0.05 ~ "up", IRdirection == "down" & pvalue < 0.05 ~ "down", TRUE ~ "none_significant"),
           p0.05_reliable = case_when(pvalue < 0.05 & reliable == "reliable" ~ "significant", TRUE ~ "none_significant"),
           p0.05_reliable_IRdirection = case_when(IRdirection == "up" & p0.05_reliable == "significant" ~ "up", IRdirection == "down" & p0.05_reliable  == "significant" ~ "down", TRUE ~ "none_significant"),
           padj0.05 = case_when(padj < 0.05 ~ "significant", TRUE ~ "none_significant"), padj0.05 = factor(padj0.05, levels = c("significant", "none_significant")),
           padj0.05_IRdirection = case_when(IRdirection == "up" & padj < 0.05 ~ "up", IRdirection == "down" & padj < 0.05 ~ "down", TRUE ~ "none_significant"),
           padj0.05_reliable = case_when(padj < 0.05 & reliable == "reliable" ~ "significant", TRUE ~ "none_significant"),
           padj0.05_reliable_IRdirection = case_when(IRdirection == "up" & padj0.05_reliable == "significant" ~ "up", IRdirection == "down" & padj0.05_reliable  == "significant" ~ "down", TRUE ~ "none_significant"),
           
           intron_type = str_split_fixed(.$Intron.GeneName.GeneID.Coords, "/", 4)[,3], intron_length = GenomicRanges::width(as.GenomicRange(Intron.GeneName.GeneID.Coords)), 
           coords = str_split_fixed(.$Intron.GeneName.GeneID.Coords, "/",4)[,4], Chr = str_split_fixed(.$Intron.GeneName.GeneID.Coords, "/|:", 6)[,4], Start = as.numeric(str_split_fixed(.$Intron.GeneName.GeneID.Coords, "/|:|-", 7)[,5]), 
           End = as.numeric(str_split_fixed(.$Intron.GeneName.GeneID.Coords, "/|:|-", 7)[,6]), Direction = str_split_fixed(.$Intron.GeneName.GeneID.Coords, ":", 3)[,3])
  print(paste(nrow(filter(IRFinder, pvalue < 0.05)), "significant IR events P < 0.05"))
  print(table(IRFinder$p0.05_IRdirection))
  print(paste(nrow(filter(IRFinder, padj < 0.05)), "significant IR events Padj < 0.05"))
  print(table(IRFinder$padj0.05_IRdirection))
  gr_introns <- as.GenomicRange(IRFinder$Intron.GeneName.GeneID.Coords)
  seqlevelsStyle(gr_introns) <- "UCSC"
  if(species == "human"){
    print("GC content based on human")
    IRFinder$gc_content <- GCcontent(Hsapiens, gr_introns)
  }else if(species == "mouse"){
    print("GC content based on mouse")
    IRFinder$gc_content <- GCcontent(Mmusculus, gr_introns)
  }else if(species == "rat"){
    print("GC content based on rat")
    IRFinder$gc_content <- GCcontent(Rnorvegicus, gr_introns)
  }
  if(species != "human"){
    IRFinder <- IRFinder %>% mutate(gene_name = toupper(gene_name))  # covert gene_names to upper case for matching to Deseq results & Human gene names
  }
  # IRFinder$phast_cons <- gscores(phast, gr_introns) # too slow to keep in function
  return(IRFinder)
}


# using output from Small Amount of Replicates via Audic and Claverie Test https://github.com/williamritchie/IRFinder/wiki/Small-Amounts-of-Replicates-via-Audic-and-Claverie-Test
cleanIRFinder.lowreplicates <- function(IRFinder = ac_who_vcp_vs_ctrl.IRFinder, deseq_res = NULL, mass_spec = NULL, species = "Hsapiens"){
  IRFinder <- IRFinder %>% mutate(gene_name = str_split_fixed(IRFinder$Intron.GeneName.GeneID, "/", 3)[,1], gene_id = str_split_fixed(IRFinder$Intron.GeneName.GeneID, "/", 3)[,2], intron_type = str_split_fixed(IRFinder$Intron.GeneName.GeneID, "/", 3)[,3],
                        coords = paste(.$Chr, ":", .$Start, "-", .$End, sep = ""), Intron.GeneName.GeneID.Coords = paste(Intron.GeneName.GeneID, "/", coords, ":", Direction, sep = ""),
                        intron_length = GenomicRanges::width(as.GenomicRange(Intron.GeneName.GeneID.Coords)),
                        # reliable / stable intron expression
                        A.IR.coverage = A.SplicesMax + A.IntronDepth,
                        B.IR.coverage = B.SplicesMax + B.IntronDepth,
                        IR.coverage = (A.IR.coverage + B.IR.coverage)/2,
                        # reliable = case_when(A.IRok != "-" ~ "unreliable", B.IRok != "-" ~ "unreliable",  TRUE ~ "reliable"),
                        reliable = case_when(A.IRok %in% c("LowCover","LowSplicing") ~ "unreliable", B.IRok %in% c("LowCover","LowSplicing") ~ "unreliable",  TRUE ~ "reliable"),
                        # reliable = case_when(((A.SplicesMax + A.IntronDepth) > 100) & ((B.SplicesMax + B.IntronDepth) > 100) ~ "reliable", TRUE ~ "unreliable"), # based on https://github.com/lbroseus/SIRFindeR/blob/master/vignettes/SIRFindeR.pdf
                        # reliable_lenient = case_when(((A.SplicesMax + A.IntronDepth) > 10) & ((B.SplicesMax + B.IntronDepth) > 10) ~ "reliable", TRUE ~ "unreliable"),
                        retained = case_when((A.IRratio >= 0.1 | B.IRratio >= 0.1) ~ "retained", TRUE ~ "spliced"),
                        reliable_retained = case_when(reliable == "reliable" & retained == "retained" ~ "retained", TRUE ~ "spliced"),
                        # differential IR
                        ## effect size
                        IRratio.diff = A.IRratio - B.IRratio, # VCP - CTRL as per Dadi at IRFinder https://github.com/williamritchie/IRFinder/issues/99
                        IR.lfc = log(IRFinder$A.IRratio/(1-IRFinder$A.IRratio)) -  log(IRFinder$B.IRratio/(1-IRFinder$B.IRratio)) %>% replace(., is.na(.), 0), # as per lucile broseus https://github.com/lbroseus/SIRFindeR/issues/1#issuecomment-668611143
                        IR.l2fc = log2(IRFinder$A.IRratio / IRFinder$B.IRratio) %>% replace(., is.na(.), 0),
                        IRdirection = case_when(IRratio.diff > 0 ~ "up", TRUE ~ "down"), # if positive then IRratio greater in VCP, negative then IRratio greater in CTRL
                        ## significance
                        p0.05 = case_when(p.diff < 0.05 ~ "significant", TRUE ~ "none_significant"), # classify significant IR events by p.diff < 0.05
                        p0.05_IRdirection = case_when(IRdirection == "up" & p.diff < 0.05 ~ "IR up", IRdirection == "down" & p.diff < 0.05 ~ "IR down", TRUE ~ "none_significant"),
                        p0.05_reliable = case_when(p.diff < 0.05 & reliable == "reliable" ~ "significant", TRUE ~ "none_significant"),
                        p0.05_reliable_IRdirection = case_when(IRdirection == "up" & p0.05_reliable == "significant" ~ "IR up",
                                                               IRdirection == "down" & p0.05_reliable  == "significant" ~ "IR down", TRUE ~ "none_significant"))
                        # # adjust for multiple comparisons
                        # padj = p.adjust(p.diff),
                        # padj0.05 = case_when(padj < 0.05 ~ "significant", TRUE ~ "none_significant"), # classify significant IR events by padj < 0.05
                        # padj0.05_IRdirection = case_when(IRdirection == "up" & padj < 0.05 ~ "IR up", IRdirection == "down" & padj < 0.05 ~ "IR down", TRUE ~ "none_significant"),
                        # padj0.05_reliable = case_when(padj < 0.05 & reliable == "reliable" ~ "significant", TRUE ~ "none_significant"),
                        # padj0.05_reliable_IRdirection = case_when(IRdirection == "up" & padj0.05_reliable == "significant" ~ "IR up",
                        #                                           IRdirection == "down" & padj0.05_reliable  == "significant" ~ "IR down",TRUE ~ "none_significant"))
  gr_introns <- as.GenomicRange(IRFinder$Intron.GeneName.GeneID.Coords)
  seqlevelsStyle(gr_introns) <- "UCSC"
  if(species == "Hsapiens"){
    print("GC content based on Hsapiens")
    IRFinder$gc_content <- GCcontent(Hsapiens, gr_introns)
  }else if(species == "Mmusculus"){
    print("GC content based on Mmusculus")
    IRFinder$gc_content <- GCcontent(Mmusculus, gr_introns)
  }else if(species == "Rnorvegicus"){
    print("GC content based on Rnorvegicus")
    IRFinder$gc_content <- GCcontent(Rnorvegicus, gr_introns)
  }
  if(species != "Hsapiens"){
  IRFinder <- IRFinder %>% mutate(gene_name = toupper(gene_name))  # covert gene_names to upper case for matching to Deseq results & Human gene names
  }
  # IRFinder$phast_cons <- gscores(phast, gr_introns) # too slow to keep in function
  if(!is.null(deseq_res)){
    IRFinder <- IRFinder %>% left_join(deseq_res, by = c("gene_id", "gene_name")) # left join deseq2 results
    IRFinder <- IRFinder %>% dplyr::mutate(DGE.direction = case_when(log2FoldChange > 0 ~ "up", TRUE ~ "down"),
                                                           direction = case_when(log2FoldChange > 0 & IRratio.diff < 0 ~ "IR down & expression up",
                                                                                 log2FoldChange > 0 & IRratio.diff > 0 ~ "IR up & expression up",
                                                                                 log2FoldChange < 0 & IRratio.diff < 0 ~ "IR down & expression down",
                                                                                 log2FoldChange < 0 & IRratio.diff > 0 ~ "IR up & expression down"))
  }
  if(!is.null(mass_spec)){
    IRFinder <- IRFinder %>%  left_join(mass_spec, by = c("gene_name"="name"))
    IRFinder <- IRFinder %>% dplyr::mutate(massspec_direction = case_when(vcp_vs_ctrl_ratio > 0 ~ "up", vcp_vs_ctrl_ratio < 0 ~ "down"),
                                           massspec_ir_direction = case_when(vcp_vs_ctrl_ratio > 0 & IRratio.diff < 0 ~ "IR down & protein up", vcp_vs_ctrl_ratio > 0 & IRratio.diff > 0 ~ "IR up & protein up",
                                                                          vcp_vs_ctrl_ratio < 0 & IRratio.diff < 0 ~ "IR down & protein down", vcp_vs_ctrl_ratio < 0 & IRratio.diff > 0 ~ "IR up & protein down"),
                                           massspec_expression_ir_direction = case_when(direction == "IR down & expression up" & massspec_direction == "up" ~ "IR down & expression up & protein up",
                                                                                        direction == "IR down & expression up" & massspec_direction == "down" ~ "IR down & expression up & protein down", 
                                                                                        direction == "IR down & expression down" & massspec_direction == "up" ~ "IR down & expression up & protein down", 
                                                                                        direction == "IR down & expression down" & massspec_direction == "down" ~ "IR down & expression down & protein down",
                                                                                        direction == "IR up & expression down" & massspec_direction == "up" ~ "IR up & expression down & protein up",
                                                                                        direction == "IR up & expression up" & massspec_direction == "up" ~ "IR up & expression up & protein up", 
                                                                                        direction == "IR up & expression up" & massspec_direction == "down" ~ "IR up & expression up & protein down",
                                                                                        direction == "IR up & expression down" & massspec_direction == "down" ~ "IR up & expression down & protein up",
                                                                                        TRUE ~ "other"))
  }
  return(IRFinder)
}


# clean the pooled replicate quantified IRFinder output (not differential IR). Needed to calculate absolute IR per condition.
cleanIRFinderQuant <- function(IRFinder = ac_who_ctrl.IRFinder, species = "Hsapiens"){
  IRFinder <- IRFinder %>% mutate(gene_name = str_split_fixed(IRFinder$Name, "/", 3)[,1], gene_id = str_split_fixed(IRFinder$Name, "/", 3)[,2], intron_type = str_split_fixed(IRFinder$Name, "/", 3)[,3],
                                  coords = paste(.$Chr, ":", .$Start, "-", .$End, sep = ""), Name.Coords = paste(Name, "/", coords, ":", Strand, sep = ""),
                                  intron_length = GenomicRanges::width(as.GenomicRange(Name.Coords)),
                                  retained = case_when(IRratio >= 0.1 ~ "retained", TRUE ~ "spliced"),
                                  ExonDepth = case_when(SpliceRight > SpliceLeft ~ SpliceRight, TRUE ~ SpliceLeft), # aka SplicesMax
                                  # IR.coverage = ExonDepth + IntronDepth,
                                  # reliable = case_when(Warnings == "-" ~ "reliable", TRUE ~ "unreliable"), # case_when((ExonDepth + IntronDepth) > 100 ~ "reliable", TRUE ~ "unreliable"),
                                  reliable = case_when(Warnings %in% c("LowCover","LowSplicing") ~ "unreliable", TRUE ~ "reliable"), # case_when((ExonDepth + IntronDepth) > 100 ~ "reliable", TRUE ~ "unreliable"),
                                  # reliable_lenient = case_when((ExonDepth + IntronDepth) > 10 ~ "reliable", TRUE ~ "unreliable"),
                                  reliable_retained = case_when(reliable == "reliable" & retained == "retained" ~ "reliable-retained", TRUE ~ "unreliable-or-spliced"))
                                  # reliable_lenient_retained = case_when(reliable_lenient == "reliable" & retained == "retained" ~ "reliable-retained", TRUE ~ "unreliable-or-spliced"))
  gr_introns <- as.GenomicRange(IRFinder$Name.Coords)
  seqlevelsStyle(gr_introns) <- "UCSC"
  if(species == "Hsapiens"){
    print("GC content based on Hsapiens")
    IRFinder$gc_content <- GCcontent(Hsapiens, gr_introns)
  }else if(species == "Mmusculus"){
    print("GC content based on Mmusculus")
    IRFinder$gc_content <- GCcontent(Mmusculus, gr_introns)
  }else if(species == "Rnorvegicus"){
    print("GC content based on Rnorvegicus")
    IRFinder$gc_content <- GCcontent(Rnorvegicus, gr_introns)
  }
  if(species != "Hsapiens"){
    IRFinder <- IRFinder %>% mutate(gene_name = toupper(gene_name))  # covert gene_names to upper case for matching to Deseq results & Human gene names
  }
  return(IRFinder)
}

# Compare the mutant and control pooled replicate quantified IRFinder output (not differential IR)
cleanIRFinderQuantCompared <- function(mutant = ac_who_vcp.IR, control = ac_who_ctrl.IR, differentialIR= NULL, deseq_res = NULL, mass_spec = NULL){
  IRFinder.control <- control %>% dplyr::select(gene_id, gene_name, coords, Chr, Start, End, Name, Strand, intron_length, gc_content,
                                                "B.IRratio" = IRratio, "B.retained" = retained, "B.reliable" = reliable, "B.reliable_retained" = reliable_retained, "B.Coverage" = Coverage, "B.IntronDepth" = IntronDepth, 
                                                "B.SpliceLeft" = SpliceLeft, "B.SpliceRight" = SpliceRight, "B.SpliceExact" = SpliceExact, "B.Warnings" = Warnings, "B.intron_type" = intron_type)
  IRFinder.mutant <- mutant %>% dplyr::select(gene_id, gene_name, coords, Chr, Start, End, Name, Strand, intron_length, gc_content,
                                              "A.IRratio" = IRratio, "A.retained" = retained, "A.reliable" = reliable, "A.reliable_retained" = reliable_retained, "A.Coverage" = Coverage, "A.IntronDepth" = IntronDepth, 
                                              "A.SpliceLeft" = SpliceLeft, "A.SpliceRight" = SpliceRight, "A.SpliceExact" = SpliceExact, "A.Warnings" = Warnings, "A.intron_type" = intron_type)
  IRFinder <- IRFinder.control %>% full_join(IRFinder.mutant, by = c("gene_id", "gene_name", "coords", "Chr", "Start", "End", "Name", "Strand", "intron_length", "gc_content"))
  IRFinder <- IRFinder %>% mutate(
    # reliable event
    A.IR.coverage = max(A.SpliceLeft, A.SpliceRight) + A.IntronDepth,
    B.IR.coverage = max(B.SpliceLeft, B.SpliceRight) + B.IntronDepth,
    IR.coverage = (A.IR.coverage + B.IR.coverage)/2,
    reliable = case_when(A.Warnings %in% c("LowCover","LowSplicing") ~ "unreliable", B.Warnings %in% c("LowCover","LowSplicing") ~ "unreliable",  TRUE ~ "reliable"),
    retained = case_when((A.IRratio >= 0.1 | B.IRratio >= 0.1) ~ "retained", TRUE ~ "spliced"),
    reliable_retained = case_when(reliable == "reliable" & retained == "retained" ~ "retained", TRUE ~ "spliced"),
    # compare IR between conditions
    IRratio.diff = A.IRratio - B.IRratio, # Mutant - CTRL as per Dadi at IRFinder https://github.com/williamritchie/IRFinder/issues/99
    IR.lfc = log(IRFinder$A.IRratio/(1-IRFinder$A.IRratio)) -  log(IRFinder$B.IRratio/(1-IRFinder$B.IRratio)) %>% replace(., is.na(.), 0), # as per lucile broseus https://github.com/lbroseus/SIRFindeR/issues/1#issuecomment-668611143
    IR.l2fc = log2(IRFinder$A.IRratio / IRFinder$B.IRratio) %>% replace(., is.na(.), 0),
    IR.direction = case_when(IRratio.diff > 0 ~ "up", TRUE ~ "down")) # if positive then IRratio greater in VCP, negative then IRratio greater in CTRL
  if(!is.null(differentialIR)){
    IRFinder.differential <- differentialIR %>% dplyr::select(Chr, Start, End, "Name" = Intron.GeneName.GeneID, "Strand" = Direction, p.diff)
    IRFinder <- IRFinder %>% left_join(IRFinder.differential, by = c("Chr", "Start", "End", "Name", "Strand")) # left join Differential IR results
    IRFinder <- IRFinder %>% dplyr::mutate(p0.05 = case_when(p.diff < 0.05 ~ "significant", TRUE ~ "none_significant"), # classify significant IR events by p.diff < 0.05
    # p0.05_IRr0.1 = case_when(p.diff < 0.05 & (A.IRratio >= 0.1 | B.IRratio >= 0.1) ~ "significant", TRUE ~ "none_significant"), # significant IR events p.diff < 0.05 and IRratio 0.1 in either sample
    # p0.05_dIRr0.1 = case_when(p.diff < 0.05 & (abs(IRratio.diff) >= 0.1) ~ "significant", TRUE ~ "none_significant"), # significant IR events p.diff < 0.05 and delta IRratio of 0.1 between samples
    p0.05_IRdirection = case_when(IR.direction == "up" & p.diff < 0.05 ~ "IR up",
                                  IR.direction == "down" & p.diff < 0.05 ~ "IR down",
                                  TRUE ~ "none_significant"),
    # p0.05_dIRr0.1_IRdirection = case_when(IRdirection == "up" & p.diff < 0.05 & (abs(IRratio.diff) >= 0.1) ~ "IR up",
    # IRdirection == "down" & p.diff < 0.05 & (abs(IRratio.diff) >= 0.1) ~ "IR down",
    # TRUE ~ "none_significant"),
    p0.05_reliable = case_when(p.diff < 0.05 & reliable == "reliable" ~ "significant", TRUE ~ "none_significant"),
    # p0.05_reliable_lenient = case_when(p.diff < 0.05 & reliable_lenient == "reliable" ~ "significant", TRUE ~ "none_significant"),
    p0.05_reliable_IR.direction = case_when(IR.direction == "up" & p0.05_reliable == "significant" ~ "IR up",
                                           IR.direction == "down" & p0.05_reliable  == "significant" ~ "IR down",
                                           TRUE ~ "none_significant"))
  }
  if(!is.null(deseq_res)){
    IRFinder <- IRFinder %>% left_join(deseq_res, by = c("gene_id", "gene_name")) # left join deseq2 results
    IRFinder <- IRFinder %>% dplyr::mutate(DGE.direction = case_when(log2FoldChange > 0 ~ "up", TRUE ~ "down"),
                                           direction = case_when(log2FoldChange > 0 & IRratio.diff < 0 ~ "IR down & expression up",
                                                                 log2FoldChange > 0 & IRratio.diff > 0 ~ "IR up & expression up",
                                                                 log2FoldChange < 0 & IRratio.diff < 0 ~ "IR down & expression down",
                                                                 log2FoldChange < 0 & IRratio.diff > 0 ~ "IR up & expression down"))
  }
  if(!is.null(mass_spec)){
    IRFinder <- IRFinder %>%  left_join(mass_spec, by = c("gene_name"="name"))
    IRFinder <- IRFinder %>% dplyr::mutate(massspec_direction = case_when(vcp_vs_ctrl_ratio > 0 ~ "up", TRUE ~ "down"),
                                           massspec_ir_direction = case_when(vcp_vs_ctrl_ratio >= 0 & IRratio.diff <= 0 ~ "IR down & protein up", vcp_vs_ctrl_ratio > 0 & IRratio.diff > 0 ~ "IR up & protein up",
                                                                             vcp_vs_ctrl_ratio < 0 & IRratio.diff < 0 ~ "IR down & protein down", vcp_vs_ctrl_ratio < 0 & IRratio.diff > 0 ~ "IR up & protein down"),
                                           massspec_expression_ir_direction = case_when(direction == "IR down & expression up" & massspec_direction == "up" ~ "IR down & expression up & protein up",
                                                                                        direction == "IR down & expression up" & massspec_direction == "down" ~ "IR down & expression up & protein down", 
                                                                                        direction == "IR down & expression down" & massspec_direction == "up" ~ "IR down & expression up & protein down", 
                                                                                        direction == "IR down & expression down" & massspec_direction == "down" ~ "IR down & expression down & protein down",
                                                                                        direction == "IR up & expression down" & massspec_direction == "up" ~ "IR up & expression down & protein up",
                                                                                        direction == "IR up & expression up" & massspec_direction == "up" ~ "IR up & expression up & protein up", 
                                                                                        direction == "IR up & expression up" & massspec_direction == "down" ~ "IR up & expression up & protein down",
                                                                                        direction == "IR up & expression down" & massspec_direction == "down" ~ "IR up & expression down & protein up",
                                                                                        TRUE ~ "other"))
  }
  return(IRFinder)
}




# IRFinder.analysis <- function(metadata = treated.who.metadata, sample.names = "sample", condition = "stimulated", ctrl = "A0", mut = "A1", batch = "NA", file.var = "Sample limsid", ge.res = treated.who.a1_vs_a0$res, 
#                               irfinder.dir = "/camp/lab/luscomben/home/shared/projects/patani-collab/astrocyte-stimulated-bulk-rnaseq/splicing/irfinder"){
#   if(batch == "NA"){
#     experiment = metadata %>% select(SampleNames = sample.names, Condition = condition, file.name = file.var) %>% remove_rownames() %>% mutate(Condition = factor(Condition,levels=c(ctrl, mut)))
#     file_paths = experiment %>% mutate(file_paths = case_when(file.exists(file.path(irfinder.dir,file.name,"IRFinder-IR-dir.txt")) ~ file.path(irfinder.dir,file.name,"IRFinder-IR-dir.txt"), TRUE ~ file.path(irfinder.dir,file.name,"IRFinder-IR-nondir.txt"))) %>% pull(file_paths)
#     print("Importing stranded IRFinder output: IRFinder-IR-dir.txt or IRFinder-IR-nondir.txt")
#       metaList = DESeqDataSetFromIRFinder(filePaths=file_paths, designMatrix=experiment, designFormula =~ Condition + Condition:IRFinder)
#   } else{
#     print(paste0("Accounting for batch: ", batch))
#     experiment = metadata %>% select(SampleNames = sample.names, Condition = condition, Batch = batch, file.name = file.var) %>% remove_rownames() %>% mutate(Condition = factor(Condition,levels=c(ctrl, mut)))
#     file_paths = experiment %>% mutate(file_paths = case_when(file.exists(file.path(irfinder.dir,file.name,"IRFinder-IR-dir.txt")) ~ file.path(irfinder.dir,file.name,"IRFinder-IR-dir.txt"), TRUE ~ file.path(irfinder.dir,file.name,"IRFinder-IR-nondir.txt"))) %>% pull(file_paths)
#     print("Importing stranded IRFinder output: IRFinder-IR-dir.txt or IRFinder-IR-nondir.txt")
#     metaList = DESeqDataSetFromIRFinder(filePaths=file_paths, designMatrix=experiment, designFormula = ~ Batch + Condition + Condition:IRFinder)
#   }
#   dds = DESeq(metaList$DESeq2Object)
#   res.ctrl = DESeq2::results(dds, name = paste0("Condition",ctrl,".IRFinderIR")) # tests IR reads vs spliced reads in ctrl
#   res.mut = DESeq2::results(dds, name = paste0("Condition",mut,".IRFinderIR")) # get IR ratio in the mut samples
#   res = DESeq2::results(dds, contrast=list(paste0("Condition",mut,".IRFinderIR"),  paste0("Condition",ctrl,".IRFinderIR")))
#   ir.res = cleanIRFinder(mutant_vs_ctrl = res, mutant = res.mut, ctrl = res.ctrl)
#   ir.ge.res = summariseIRjoinGE(cleanIRFinder = ir.res, GE = ge.res)
#   return(list(irfinder=ir.res,ir.ge=ir.ge.res))
# }



mergeIRFinderQuantDESeq2 <- function(cleanIRFinderQuantCompared = NULL, deseq_res = NULL, mass_spec = NULL){
  cleanIRFinderQuantCompared.GE <- cleanIRFinderQuantCompared  %>% filter(!A.Warnings %in% c("LowCover","LowSplicing"), !B.Warnings %in% c("LowCover","LowSplicing")) %>%
    group_by(gene_name) %>% dplyr::summarise(all_intron_lengths = sum(intron_length), introns = dplyr::n(), 
                                             IRratio.diff.norm = sum( IRratio.diff * (intron_length / all_intron_lengths)), 
                                             IRratio.diff.mean = mean(IRratio.diff, na.rm = TRUE),
                                             IRratio.diff.max = IRratio.diff[which.max(abs(IRratio.diff))],
                                             # IRl2fc.norm = sum( IR.l2fc * (intron_length / all_intron_lengths)),
                                             # IRl2fc.mean = mean(IR.l2fc, na.rm = TRUE),
                                             # IRl2fc.max = IR.l2fc[which.max(abs(IR.l2fc))],
                                             IR.pval = min(p.diff, na.rm = TRUE), IR.p0.05_reliable = max(p0.05_reliable, na.rm = TRUE)) %>% ungroup %>% mutate(IR.pval = case_when(IR.pval == Inf ~ NaN, TRUE ~ IR.pval))
  if(!is.null(deseq_res)){
    deseq_res <- deseq_res %>% group_by(gene_name) %>% slice_max(baseMean, n = 1, with_ties = FALSE)
    cleanIRFinderQuantCompared.GE <- cleanIRFinderQuantCompared.GE %>% 
      full_join(filter(dplyr::select(deseq_res, gene_name, "mRNA.lfc" = log2FoldChange, "mRNA.FDR" = padj, "mRNA.baseMean" = baseMean)))
  }
  if(!is.null(mass_spec)){
    cleanIRFinderQuantCompared.GE <- cleanIRFinderQuantCompared.GE %>% 
      full_join(dplyr::select(mass_spec, "gene_name" = name, protein.lfc = vcp_vs_ctrl_ratio, "protein.pval" = vcp_vs_ctrl_p.val), by = "gene_name")
  }
  return(cleanIRFinderQuantCompared.GE)
}


cleanVASTTOOLS <- function(diff = ac_who_vcp_vs_ctrl_diff, comp = ac_who_vcp_vs_ctrl_compare){
  txt <- diff %>% left_join(dplyr::select(comp, -starts_with("nuclear_"),-starts_with("cytoplasmic_")), by = c("GENE", "EVENT"))
  txt <- txt %>% mutate(type = case_when(grepl("HsaEX", EVENT) ~ "exon skipping",
                                         grepl("HsaINT", EVENT) ~ "intron retention",
                                         grepl("HsaALTA", EVENT) ~ "alternate 3'",
                                         grepl("HsaALTD", EVENT) ~ "alternate 5'"))
  txt <- txt %>% mutate(compare.direction = case_when(dPSI > 0 ~ "up", TRUE ~ "down"),
                        dPSI0.1 = case_when(dPSI > 0.1 ~ "significant", TRUE ~ "none_significant"), 
                        dPSI0.1_direction = case_when(compare.direction == "up" & dPSI0.1 > 0.1 ~ "up",
                                                      compare.direction == "down" & dPSI0.1 < -0.1 ~ "down",
                                                      TRUE ~ "none_significant"),
                        diff.direction = case_when(E.dPsi. > 0 ~ "up", TRUE ~ "down"), # if positive then IRratio greater in VCP, negative then IRratio greater in CTRL
                        MV0.1 = case_when(MV.dPsi._at_0.95 > 0.1 ~ "significant", TRUE ~ "none_significant"), 
                        MV0.1_dPSI0.1 = case_when(MV.dPsi._at_0.95 > 0.1 & abs(E.dPsi.) >= 0.1 ~ "significant", TRUE ~ "none_significant"), # significant IR events p.diff < 0.05 and delta IRratio of 0.1 between samples
                        MV0.1_direction = case_when(diff.direction == "up" & MV.dPsi._at_0.95 > 0.1 ~ "up",
                                                    diff.direction == "down" & MV.dPsi._at_0.95 > 0.1 ~ "down",
                                                    TRUE ~ "none_significant"),
                        MV0.1_dPSI0.1_direction = case_when(diff.direction == "up" & MV.dPsi._at_0.95 > 0.1 & (abs(E.dPsi.) >= 0.1) ~ "up",
                                                            diff.direction == "down" & MV.dPsi._at_0.95 > 0.1 & (abs(E.dPsi.) >= 0.1) ~ "down",
                                                            TRUE ~ "none_significant"))
  fromSplit <- str_split_fixed(txt$COORD, ":", 2)
  txt$Chr <- fromSplit[,1] %>% gsub("chr", "", .)
  coord <- str_split_fixed(fromSplit[,2], "-", 2)
  txt$Start <- coord[,1] %>% as.numeric
  txt$End <- coord[,2] %>% as.numeric
  return(txt)
}


# labels = c("VIM", "COL1A1", "MYL6", "DDX5", "FN1", "FLNA", "HNRNPK", "FUS", "HNRNPL", "SFPQ", "EMP1", "HNRNPH1", "OSMR", "NONO", "HLA-A", "HLA-B", "HLA-C", "STAT3", "ANXA2", "SPARC", "HNRNPU", "HNRNPA2B1", "B2M", "PPIA", "LRP1","SRRM2", "SRSF5", "TNC", "SLC44A2", "CXCL2")

histogram_plot_irfinder <- 
  function(data = irfinder.ac_nuc_vcp_vs_ctrl, reliabilityThr = 10, legend = TRUE){
    data <- data %>% mutate(reliable.threshold = case_when(baseMean < reliabilityThr ~ "unreliable", TRUE ~ "reliable"))
    if(legend == TRUE){
      ggplot(filter(data, abs(IRratio.diff) > 0), aes(x = IRratio.diff, fill = reliable.threshold)) + 
        geom_histogram(position="dodge",binwidth=0.03) +  theme_bw() + 
        theme(strip.background = element_rect(colour="black", fill="white"), panel.grid = element_blank(), legend.position = "top", legend.box = "horizontal", legend.margin=margin(0.2,0.2,0.2,0.2), legend.box.margin=margin(-5,-5,-5,-5)) + labs(fill = expression(Coverage)) + 
        geom_vline(xintercept = 0.1, linetype = 3, colour = "darkgrey")  + geom_vline(xintercept = -0.1, linetype = 3, colour = "darkgrey")  + geom_vline(xintercept = 0, linetype = 1, colour = "darkgrey")  + 
        geom_hline(yintercept = 0, linetype = 1, colour = "black") +
        scale_fill_manual( values = c("dodgerblue2", "firebrick2") ) + scale_colour_manual( values = c("dodgerblue2", "firebrick2") ) +
        coord_cartesian(xlim=c(-0.5,0.5)) + xlab(label = expression(Delta~IR~ratio) ) + ylab(label = expression(IR~event~count) )
        } else{
      ggplot(filter(data, abs(IRratio.diff) > 0), aes(x = IRratio.diff, fill = reliable.threshold)) + 
        geom_histogram(position="dodge",binwidth=0.03) +  theme_bw() + 
        theme(strip.background = element_rect(colour="black", fill="white"), panel.grid = element_blank(), legend.position = "none") + 
        geom_vline(xintercept = 0.1, linetype = 3, colour = "darkgrey")  + geom_vline(xintercept = -0.1, linetype = 3, colour = "darkgrey")  + geom_hline(yintercept = 0, linetype = 1, colour = "black") +
        scale_fill_manual( values = c("dodgerblue2", "firebrick2") ) + scale_colour_manual( values = c("dodgerblue2", "firebrick2") )+
        coord_cartesian(xlim=c(-0.5,0.5)) + xlab(label = expression(Delta~IR~ratio) ) + ylab(label = expression(IR~event~count) )
          }
  }


histogram_plot_irfinder.log2FoldChange <- 
  function(data = irfinder.ac_nuc_vcp_vs_ctrl, reliabilityThr = 10, legend = TRUE, y.position = 100000){
    data <- data %>% mutate(reliable.threshold = case_when(baseMean < reliabilityThr ~ "unreliable", TRUE ~ "reliable"))
    if(legend == TRUE){
      ggplot(filter(data, abs(log2FoldChange) > 0), aes(x = log2FoldChange, fill = reliable.threshold)) + 
        geom_histogram(position="dodge",binwidth=0.1) +  theme_bw() + 
        theme(strip.background = element_rect(colour="black", fill="white"), panel.grid = element_blank(), legend.position = "top", legend.box = "horizontal", legend.margin=margin(0.2,0.2,0.2,0.2), legend.box.margin=margin(-5,-5,-5,-5)) + labs(fill = expression(Coverage)) + 
        geom_vline(xintercept = 0.1, linetype = 3, colour = "darkgrey")  + geom_vline(xintercept = -0.1, linetype = 3, colour = "darkgrey")  + geom_vline(xintercept = 0, linetype = 1, colour = "darkgrey")  + 
        geom_hline(yintercept = 0, linetype = 1, colour = "black") +
        scale_fill_manual( values = c("dodgerblue2", "firebrick2") ) + scale_colour_manual( values = c("dodgerblue2", "firebrick2") ) +
        coord_cartesian(xlim=c(-3,3)) + xlab(label = expression(LFC~IR) ) + ylab(label = expression(IR~event~count) ) #+
      # geom_segment(aes(x = 1, xend = 6, y= y.position, yend= y.position), arrow=arrow(length=unit(0.3,"cm")), color = "black") + geom_segment(aes(x = 1, xend = -6, y= y.position, yend= y.position),arrow=arrow(length=unit(0.3,"cm")), color = "black") +
      # annotate("text", x = 5, y = y.position + 10000, label = "IR Up", colour = "black") +   annotate("text", x = -5, y =  y.position + 10000, label = "IR Down", colour = "black")
    } else{
      ggplot(filter(data, abs(log2FoldChange) > 0), aes(x = log2FoldChange, fill = reliable.threshold)) + 
        geom_histogram(position="dodge",binwidth=0.1) +  theme_bw() + 
        theme(strip.background = element_rect(colour="black", fill="white"), panel.grid = element_blank(), legend.position = "none") + 
        geom_vline(xintercept = 0.1, linetype = 3, colour = "darkgrey")  + geom_vline(xintercept = -0.1, linetype = 3, colour = "darkgrey")  + geom_hline(yintercept = 0, linetype = 1, colour = "black") +
        scale_fill_manual( values = c("dodgerblue2", "firebrick2") ) + scale_colour_manual( values = c("dodgerblue2", "firebrick2") )+
        coord_cartesian(xlim=c(-3,3)) + xlab(label = expression(LFC~IR) ) + ylab(label = expression(IR~event~count) ) +
        geom_segment(aes(x = 1, xend = 6, y= y.position, yend= y.position), arrow=arrow(length=unit(0.3,"cm")), color = "black") + geom_segment(aes(x = 1, xend = -6, y= y.position, yend= y.position),arrow=arrow(length=unit(0.3,"cm")), color = "black") +
        annotate("text", x = 5, y = y.position + 10000, label = "IR Up", colour = "black") +   annotate("text", x = -5, y =  y.position + 10000, label = "IR Down", colour = "black")
    }
  }


histogram_plot_irfinder.lowreplicates <- 
  function(ctrl = ac_who_ctrl.IR, mutant = ac_who_vcp.IR, ctrl.name = "ctrl", mutant.name = "VCP", legend = TRUE){
    ctrl.filt <- ctrl %>% dplyr::select(Chr, Start, End, Name, Strand, IRratio, gene_name, reliable) %>% mutate(condition = ctrl.name)
    mutant.filt <- mutant %>% dplyr::select(Chr, Start, End, Name, Strand, IRratio, gene_name, reliable) %>% mutate(condition = mutant.name)
    ctrl_and_mutant.filt <- ctrl.filt %>% bind_rows(mutant.filt) %>% mutate(condition = factor(condition, levels = c(ctrl.name, mutant.name)))
    if(legend == TRUE){
    ggplot(filter(ctrl_and_mutant.filt, IRratio > 0), aes(x = IRratio, fill = reliable)) + 
      geom_histogram(position="dodge",binwidth=0.03) +  theme_bw() + 
      theme(strip.background = element_rect(colour="black", fill="white"), panel.grid = element_blank(), legend.position = "top", legend.box = "horizontal", legend.margin=margin(0.2,0.2,0.2,0.2), legend.box.margin=margin(-5,-5,-5,-5)) + labs(fill = expression(Coverage)) + 
      geom_vline(xintercept = 0.1, linetype = 3, colour = "darkgrey")  + geom_hline(yintercept = 0, linetype = 1, colour = "black") +
      facet_grid(. ~ condition) + scale_fill_manual( values = c("dodgerblue2", "firebrick2") ) + scale_colour_manual( values = c("dodgerblue2", "firebrick2") )
    } else{
    ggplot(filter(ctrl_and_mutant.filt, IRratio > 0), aes(x = IRratio, fill = reliable)) + 
      geom_histogram(position="dodge",binwidth=0.03) +  theme_bw() + 
      theme(strip.background = element_rect(colour="black", fill="white"), panel.grid = element_blank(), legend.position = "none") + 
      geom_vline(xintercept = 0.1, linetype = 3, colour = "darkgrey")  + geom_hline(yintercept = 0, linetype = 1, colour = "black") +
        facet_grid(. ~ condition)  + scale_fill_manual( values = c("dodgerblue2", "firebrick2") ) + scale_colour_manual( values = c("dodgerblue2", "firebrick2") )
    }
  }



IRratio.diff_histogram_plot_irfinder.lowreplicates <- 
  function(data = ac_who_vcp_vs_ctrl, reliabilityThr = 100){
    if(reliabilityThr == 100){
      ggplot(data, aes(x = IRratio.diff, fill = reliable)) + 
        geom_histogram(colour="black", position="dodge",binwidth=0.03) +  theme_bw() + 
        theme(strip.background = element_rect(colour="black", fill="white"), panel.grid = element_blank(), legend.position = "top") + labs(fill = expression(Coverage)) + 
        geom_vline(xintercept = 0.1, linetype = 2, colour = "darkred") + geom_vline(xintercept = -0.1, linetype = 2, colour = "darkred")  + geom_vline(xintercept = 0, linetype = 1, colour = "black") +
        scale_fill_manual( values = c("dodgerblue", "firebrick2") )
    } else{
      ggplot(data, aes(x = IRratio.diff, fill = reliable_lenient)) + 
        geom_histogram(colour="black", position="dodge",binwidth=0.01) +  theme_bw() + 
        theme(strip.background = element_rect(colour="black", fill="white"), panel.grid = element_blank(), legend.position = "top") + labs(fill = expression(Coverage)) + 
        geom_vline(xintercept = 0.1, linetype = 2, colour = "darkred") + geom_vline(xintercept = -0.1, linetype = 2, colour = "darkred")  + geom_vline(xintercept = 0, linetype = 1, colour = "black")  +
        scale_fill_manual( values = c("dodgerblue", "firebrick2") )
    }
  }


reliabile_plot_irfinder.lowreplicates <- 
  function(data = ac_who_ctrl, reliabilityThr = 100){
    ggplot(filter(data, IRratio > 0), aes( x = IRratio, y = (ExonDepth + IntronDepth),  label = gene_name)) + 
    theme(strip.background = element_rect(colour="black", fill="white"), panel.grid = element_blank()) + theme_bw() +
    geom_point(colour = "purple", alpha = 0.5, size = 0.5) +  theme_bw() +  xlab( expression(IR~ratio)  ) +  ylab( expression(log[2]~Intron~Exon~abundance) ) + scale_y_continuous(trans='log2') +
    geom_hline(yintercept = log2(reliabilityThr), linetype = 3, colour = "darkred") + geom_vline(xintercept = 0.1, linetype = 3, colour = "darkred")
  }


ma_plot_irfinder <- 
  function(data = irfinder.mn.nuc_vs_cyt.vcp_vs_ctrl, labels_list = irfinder.mn.nuc_vs_cyt.vcp_vs_ctrl.ma.labels, padj = TRUE){
    data <- data %>% mutate(p0.05_reliable_IRdirection = factor(p0.05_reliable_IRdirection, levels = c("up", "down", "none_significant")))
    if(padj == TRUE){
    labels_maplot <- data %>% filter(gene_name %in% labels_list, padj0.05_reliable == "significant") %>% arrange(-(abs(log2FoldChange))) %>% distinct(gene_name, .keep_all = TRUE)
    ggplot(data = data, aes(x = log10(baseMean), y = log2FoldChange, colour = padj0.05_reliable_IRdirection)) +  # NB if padj== NA then will not show
      geom_point(size = 0.5) +
      scale_colour_manual(values = c("up"="firebrick2", "down"="dodgerblue2", "none_significant"="darkgray") ) +  
      theme_bw() +  theme(strip.background = element_rect(colour="black", fill="white"), panel.grid = element_blank(), legend.position = "none") +
      xlab(expression( Log[10]~IR~event~coverage )) + ylab(expression(log[2]~Fold~Change)) + 
      geom_text_repel(data = labels_maplot, aes(x = log10(baseMean), y = log2FoldChange, label = gene_name), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0, size=2.7, max.overlaps = Inf) +
      geom_hline(yintercept = 0, linetype = 1) + geom_vline(xintercept = log10(10), linetype = 3, colour = "darkgrey")
    }
    else{
      labels_maplot <- data %>% filter(gene_name %in% labels_list, p0.05_reliable == "significant") %>% arrange(-(abs(log2FoldChange))) %>% distinct(gene_name, .keep_all = TRUE)
      ggplot(data = data, aes(x = log10(baseMean), y = log2FoldChange, colour = p0.05_reliable_IRdirection)) +  # NB if pvalue == NA then will not show
        geom_point(size = 0.5) +
        scale_colour_manual(values = c("up"="firebrick2", "down"="dodgerblue2", "none_significant"="darkgray") ) +  
        theme_bw() +  theme(strip.background = element_rect(colour="black", fill="white"), panel.grid = element_blank(), legend.position = "none") +
        xlab(expression( Log[10]~IR~event~coverage )) + ylab(expression(log[2]~Fold~Change)) + 
        geom_text_repel(data = labels_maplot, aes(x = log10(baseMean), y = log2FoldChange, label = gene_name), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0, size=2.7, max.overlaps = Inf) +
        geom_hline(yintercept = 0, linetype = 1) + geom_vline(xintercept = log10(10), linetype = 3, colour = "darkgrey")
    }
  }

ma_plot_irfinder.IRratio.diff <- 
  function(data = irfinder.ac_nuc_vcp_vs_ctrl, labels_list = irfinder.ac_nuc_vcp_vs_ctrl.ma.labels){
    labels_maplot <- data %>% filter(gene_name %in% labels_list) %>% filter(p0.05_reliable == "significant") %>% arrange(-(abs(IRratio.diff))) %>% distinct(gene_name, .keep_all = TRUE)
    ggplot(data = data, aes(x = log10(baseMean), y = IRratio.diff, colour = p0.05_reliable_IRdirection)) +  # NB if pvalue == NA then will not show
      geom_point(size = 0.5) +
      scale_colour_manual( values = c("firebrick2", "dodgerblue2", "darkgray") ) +  theme_bw() +  
      theme(panel.grid = element_blank(), legend.title = element_blank(), legend.position = "none") +
      xlab(expression( Log[10]~IR~event~coverage )) + ylab(expression(Delta~IR~ratio)) + 
      geom_text_repel(data = labels_maplot, aes(x = log10(baseMean), y = IRratio.diff, label = gene_name), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0, size=2.7, max.overlaps = Inf) +
      geom_hline(yintercept = 0, linetype = 1) + geom_vline(xintercept = log10(10), linetype = 3, colour = "darkgrey")
  }


ma_plot_irfinder.lowreplicates <- 
  function(data = ac_ctrl_cyt_vs_nuc, xlabs = "Log[10] exon + intron coverage", labels_list = astrocyte.reactivity.labels){
    ggmaplot <- data %>% mutate(baseMean = log10(IR.coverage))
    labels_maplot <- data %>% mutate(baseMean = log10(IR.coverage)) %>% filter(p0.05_reliable == "significant") %>% filter(gene_name %in% labels_list) %>% arrange( -(abs(IRratio.diff)) ) %>% distinct(gene_name, .keep_all = TRUE)
    ggplot(data = ggmaplot, aes(x = baseMean, y = IRratio.diff, colour = p0.05_reliable_IRdirection)) +  geom_point(size = 0.5) + #xlim(1,5) + ylim(-3,3) +
        scale_colour_manual( values = c("firebrick2", "dodgerblue2", "darkgray") ) +  theme_bw() +  
        theme(panel.grid = element_blank(), axis.title.y = element_text(size = 10), legend.position = "none") +
        xlab(expression( Log[10]~coverage )) + 
        geom_text_repel(data = labels_maplot, aes(x = baseMean, y = IRratio.diff, label = gene_name), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0,
                        size=2.5, max.overlaps = Inf) +
        geom_hline(yintercept = 0, linetype = 2) + geom_vline(xintercept = 1, linetype = 2)
  }


volcano_plot_irfinder <- 
  function(data = ac_ctrl_cyt_vs_nuc, xlabs = "Log[10] exon + intron coverage", labels_list = astrocyte.reactivity.labels, reliability = "reliable"){
    if(reliability == "reliable"){
      labels_plot <- data %>% filter(p0.05_reliable == "significant") %>% filter(gene_name %in% labels_list) %>% arrange( -(abs(IRratio.diff)) ) %>% distinct(gene_name, .keep_all = TRUE)
      ggplot(data = data, aes(x = IRratio.diff, y =  -log10(p.diff), colour = p0.05_reliable_IRdirection)) +  geom_point(size = 0.5) + #xlim(1,5) + ylim(-3,3) +
        scale_colour_manual( values = c("firebrick2", "dodgerblue", "darkgray") ) +  theme_bw() +  theme(panel.grid = element_blank(), legend.title = element_blank(), legend.position = "none") +
        ylab(expression( -log[10]~FDR )) + #xlab( expression( IR~log[2]~fold~change~CTRL~cytoplasmic:nuclear)) +
        geom_text_repel(data = labels_plot, aes(x = IRratio.diff, y = -log10(p.diff), label = gene_name), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0, max.overlaps = Inf) +
        geom_hline(yintercept =1.30, linetype = 2) + geom_vline(xintercept =0, linetype = 2)
    } else if(reliability == "reliable_lenient"){
      labels_plot <- data %>% filter(p0.05_reliable_lenient == "significant") %>% filter(gene_name %in% labels_list) %>% arrange( -(abs(IRratio.diff)) ) %>% distinct(gene_name, .keep_all = TRUE)
      ggplot(data = data, aes(x = IRratio.diff, y =  -log10(p.diff), colour = p0.05_reliable_lenient_IRdirection)) +  geom_point(size = 0.5) + #xlim(1,5) + ylim(-3,3) +
        scale_colour_manual( values = c("firebrick2", "dodgerblue", "darkgray") ) +  theme_bw() +  theme(panel.grid = element_blank(), legend.title = element_blank(), legend.position = "none") +
        ylab(expression( -log[10]~FDR )) + #xlab( expression( IR~log[2]~fold~change~CTRL~cytoplasmic:nuclear)) +
        geom_text_repel(data = labels_plot, aes(x = IRratio.diff, y = -log10(p.diff), label = gene_name), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0, max.overlaps = Inf) +
        geom_hline(yintercept =1.30, linetype = 2) + geom_vline(xintercept =0, linetype = 2)
    } else{
      labels_plot <- data %>% filter(p0.05 == "significant") %>% filter(gene_name %in% labels_list) %>% arrange( -(abs(IRratio.diff)) ) %>% distinct(gene_name, .keep_all = TRUE)
      ggplot(data = data, aes(x = IRratio.diff, y =  -log10(p.diff), colour = p0.05_reliable_IRdirection)) +  geom_point(size = 0.5) + #xlim(1,5) + ylim(-3,3) +
      scale_colour_manual( values = c("firebrick2", "dodgerblue", "darkgray") ) +  theme_bw() +  theme(panel.grid = element_blank(), legend.title = element_blank(), legend.position = "none") +
      ylab(expression( -log[10]~FDR )) + #xlab( expression( IR~log[2]~fold~change~CTRL~cytoplasmic:nuclear)) +
      geom_text_repel(data = labels_plot, aes(x = IRratio.diff, y = -log10(p.diff), label = gene_name), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0, max.overlaps = Inf) +
      geom_hline(yintercept =1.30, linetype = 2) + geom_vline(xintercept =0, linetype = 2)
  }
  }

volcano_plot_irfinder_p0.05 <- 
  function(data = ac_ctrl_cyt_vs_nuc, xlabs = "Log[10] exon + intron coverage", labels_list = labels){
    labels_plot <- data %>% filter(p0.05 == "significant") %>% filter(gene_name %in% labels_list) %>% arrange( -(abs(IRratio.diff)) ) %>% distinct(gene_name, .keep_all = TRUE)
    ggplot(data = data, aes(x = IRratio.diff, y =  -log10(p.diff), colour = p0.05_IRdirection)) +  geom_point(size = 0.5) + #xlim(1,5) + ylim(-3,3) +
      scale_colour_manual( values = c("firebrick2", "dodgerblue", "darkgray") ) +  theme_bw() +  theme(panel.grid = element_blank(), legend.title = element_blank(), legend.position = "none") +
      ylab(expression( -log[10]~FDR )) + #xlab( expression( IR~log[2]~fold~change~CTRL~cytoplasmic:nuclear)) +
      geom_text_repel(data = labels_plot, aes(x = IRratio.diff, y = -log10(p.diff), label = gene_name), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0, max.overlaps = Inf) +
      geom_hline(yintercept =1.30, linetype = 2) + geom_vline(xintercept =0, linetype = 2)
  }


# volcano_plot_irfinder(data = ac_ctrl_cyt_vs_nuc, labels_list = labels) + xlab( expression( IR~log[2]~fold~change~CTRL~cytoplasmic:nuclear )) + #xlim(1,5) + ylim(-3,3)


scatter_plot_irfinder <-
  function(data = ctrl_vcp_scatter, labels_list = labels, x = IR.lfc_ctrl, y = IR.lfc_vcp){
    # labels_plot <- data %>% filter(gene_name %in% labels) %>% arrange( -(abs(x)) ) %>% distinct(gene_name, .keep_all = TRUE)
    ggplot(data, aes( x = x, y = y)) + geom_point(colour = "purple", alpha = 0.5, size = 0.5) +  theme_bw() +  
    # geom_text_repel(data = labels_plot, aes(x = x, y = y, label = gene_name), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0.5, max.overlaps = Inf) +
    theme(strip.background = element_rect(colour="black", fill="white"), panel.grid = element_blank()) +
    geom_hline(yintercept = 0, linetype = 3, colour = "darkgrey") + geom_vline(xintercept = 0, linetype = 3, colour = "darkgrey") + 
    geom_smooth(data = data, aes(x = x, y = y), method = "lm", se = FALSE, show.legend = TRUE, colour = "black") # + geom_abline(slope = -0.11, intercept = 0, linetype = 1)
  }


# Coverage Plots ---------------
# ggbio http://www.sthda.com/english/wiki/ggbio-visualize-genomic-data   https://www.bioconductor.org/packages/release/bioc/vignettes/ggbio/inst/doc/ggbio.pdf
# Gviz is better # https://blog.liang2.tw/posts/2016/01/plot-seq-depth-gviz/ http://www.sthda.com/english/wiki/gviz-visualize-genomic-data # Dave Tang blog: https://davetang.org/muse/2013/10/03/using-gviz/
options(ucscChromosomeNames=FALSE) # crucial to ensure AlignmentsTrack plot works
gtrack <- GenomeAxisTrack() # genome axis track
# itrack <- IdeogramTrack(genome = "hg38", chromosome = 1) # add chromosome ideogram # Error in readHTMLTable(url)[[1L]] : subscript out of bounds
# Transcript annotations
### TURN ON ONLY WHEN NEEDED
# txdb <- makeTxDbFromGFF("/camp/home/ziffo/home/genomes/ensembl/Homo_sapiens.GRCh38.99.chr_patch_hapl_scaff.gtf", format="gtf") # make txdb from GTF - turn on only when needed
# seqlevels(txdb, pruning.mode="coarse") <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y") # removes all the other chromsome annotations & scaffolds which are also in the BAM

getGeneIDsFromTxDb <- function(gr, txdb){
  stopifnot(is(gr, "GRanges"))
  stopifnot(length(gr)>0)
  stopifnot(is(txdb, "TxDb"))
  if(length(gr)>1){
    warning("The length of gr is greater than 1. Only first genomic location will be used.")
    gr <- gr[1]
  }
  genes <- genes(txdb, columns="gene_id")
  genes <- subsetByOverlaps(genes, gr)
  return(genes$gene_id)
}

# Track viewer https://www.bioconductor.org/packages/release/bioc/vignettes/trackViewer/inst/doc/trackViewer.html#introduction
track_viewer_nuc_cyt <- function(gene.name = "UBA6", mut.nuc.bam = ac_nuc_vcp.bam, ctrl.nuc.bam = ac_nuc_ctrl.bam, mut.cyt.bam = ac_cyt_vcp.bam, ctrl.cyt.bam = ac_cyt_ctrl.bam,
                                 plot_transcript_number = 6, print2screen = TRUE, save = TRUE, limit.max.ctrl = "NA", limit.min.ctrl = "NA", limit.max.vcp = "NA", limit.min.vcp = "NA"){
  gtf.gene <- gtf %>% filter(gene_name == gene.name) %>% group_by(ensemblID, gene_name, strand, seqnames) %>% summarise(start = min(start), end = max(end), width = max(width))
  print(gtf.gene)
  # cat(blue(gtf.gene))
  if (nrow(gtf.gene) == 0){
    stop("ensembl ID not found in GTF")
  }
  gene_name <-  gtf.gene$gene_name
  thechr <-  gtf.gene$seqnames  # without chr before number
  st <- gtf.gene$start
  en <- gtf.gene$end
  strand <- gtf.gene$strand
  gr <- GRanges(thechr, IRanges(st - 100, en + 100), strand = strand) # Granges for retained intron with padding either side
  
  CTRL <- importBam(ctrl.nuc.bam, ctrl.cyt.bam, ranges=gr, pairs = TRUE) # import BAMs
  VCP <- importBam(mut.nuc.bam, mut.cyt.bam, ranges=gr, pairs = TRUE)
  CTRL$dat <- coverageGR(CTRL$dat)    # calculate coverage
  VCP$dat <- coverageGR(VCP$dat)
  
  gene <- geneModelFromTxdb(txdb, org.Hs.eg.db, gr = gr, strand = as.character(strand))   # Build Gene model
  
  # Set track arguments
  optSty <- optimizeStyle(trackList(VCP, CTRL, gene))
  trackList <- optSty$tracks
  viewerStyle <- optSty$style  # viewerStyle <- trackViewerStyle()
  setTrackViewerStyleParam(viewerStyle, "xaxis", FALSE)
  setTrackViewerStyleParam(viewerStyle, "margin", c(.01, .10, -0.1, .01))
  
  # optimise coverage plots
  for(i in 1:2) {
    setTrackXscaleParam(trackList[[i]], "gp", list(cex=0.9))
    setTrackStyleParam(trackList[[i]], "height", .36)
    setTrackStyleParam(trackList[[i]], "color", c("firebrick2", "dodgerblue2"))
    setTrackStyleParam(trackList[[i]], "ylabgp", list(cex=0.8))
  }
  if (length(trackList) > (plot_transcript_number + 2) ){
        print(paste("There are ", length(trackList), " ensembl transcripts. Plotting only the first ", plot_transcript_number)) 
        for(i in (plot_transcript_number + 3):length(trackList)) {
          print(trackList[[i]]$name) 
        }
        idx = c((plot_transcript_number + 3):length(trackList))
        trackList <- trackList[-idx] # only plot first selected transcript tracks - remove remaining transcript tracks
      }
      # modify the transcript tracks
      for(i in 3:length(trackList)) {
        setTrackStyleParam(trackList[[i]], "height", .28/(plot_transcript_number+2))
        setTrackStyleParam(trackList[[i]], "ylabgp", list(cex=0))
        setTrackStyleParam(trackList[[i]], "color", "black")
      }
  
  if (limit.max.ctrl != "NA"){ 
  setTrackStyleParam(trackList[[2]], "ylim", c(limit.min.ctrl, limit.max.ctrl))  ## Adjust the limit of CTRL y-axis
  }
  if (limit.max.vcp != "NA"){
    setTrackStyleParam(trackList[[1]], "ylim", c(limit.min.vcp, limit.max.vcp))  ## Adjust the limit of VCP y-axis
  }

  if (print2screen == TRUE ){
  print("Printing to screen.")
  vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
  grid.text("Nuclear reads", x=.85, y=.76, just="bottom", gp=gpar(fontsize=6, col = "firebrick2"))
  grid.text("Cytoplasm reads", x=.85, y=.41, just="bottom", gp=gpar(fontsize=6, col = "dodgerblue2"))
  grid.text("Nuclear reads", x=.85, y=.36, just="bottom", gp=gpar(fontsize=6, col = "firebrick2"))
  grid.text("Cytoplasm reads", x=0.85, y=.01, just="bottom", gp=gpar(fontsize=6, col = "dodgerblue2"))
  grid.text(gene_name, x=.06, y=.88, just="bottom", gp=gpar(fontsize=10))#, rot=90)
  dev.off()
  }
  if (save == TRUE ){
  print(paste0("Saving figure as /camp/lab/luscomben/home/users/ziffo/projects/motor-neuron-mislocalisation/expression/coverage/", gene_name, ".nuc_cyt_coverage_plot.png"))
  png(file = paste0("/camp/lab/luscomben/home/users/ziffo/projects/motor-neuron-mislocalisation/expression/coverage/", gene_name, ".nuc_cyt_coverage_plot.png"), height = 3.5, width = 6.5, units ="in", res = 300)
  vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
  grid.text("Nuclear reads", x=.85, y=.76, just="bottom", gp=gpar(fontsize=6, col = "firebrick2"))
  grid.text("Cytoplasm reads", x=.85, y=.41, just="bottom", gp=gpar(fontsize=6, col = "dodgerblue2"))
  grid.text("Nuclear reads", x=.85, y=.36, just="bottom", gp=gpar(fontsize=6, col = "firebrick2"))
  grid.text("Cytoplasm reads", x=0.85, y=.01, just="bottom", gp=gpar(fontsize=6, col = "dodgerblue2"))
  grid.text(gene_name, x=.06, y=.88, just="bottom", gp=gpar(fontsize=10))#, rot=90)
  dev.off()
}}

# homo_sapiens.gtf <- import("/camp/home/ziffo/home/genomes/ensembl/Homo_sapiens.GRCh38.99.gtf") %>% as_tibble
# gtf.3utr <- homo_sapiens.gtf %>% filter(type == "three_prime_utr", transcript_biotype == "protein_coding")  %>% group_by(gene_id) %>% top_n(end, n = 2) %>% select(transcript_id, gene_id, gene_name, seqnames, start, end, width, strand, ccds_id) %>% unique()
# gtf.exon <- homo_sapiens.gtf %>% filter(type == "exon", transcript_biotype == "protein_coding") %>% group_by(gene_id) %>% top_n(end, n = 2) %>% select(transcript_id, gene_id, gene_name, seqnames, start, end, width, strand, ccds_id) %>% unique()
# gtf.apa <- homo_sapiens.gtf %>% filter(type %in% c("three_prime_utr", "exon"), transcript_biotype == "protein_coding") %>% group_by(gene_id) %>% top_n(end, n = 2) %>% select(transcript_id, gene_id, gene_name, seqnames, start, end, width, strand, ccds_id) %>% unique()
# saveRDS(gtf.apa, "/camp/home/ziffo/home/genomes/ensembl/Homo_sapiens.GRCh38.99.gtf.apa.coords.rds")
gtf.apa <- readRDS("/camp/home/ziffo/home/genomes/ensembl/Homo_sapiens.GRCh38.99.gtf.apa.coords.rds")

track_viewer_nuc_cyt_apa <- function(gene.name = "ZSWIM7", mut.nuc.bam = mn_nuc_vcp.bam, ctrl.nuc.bam = mn_nuc_ctrl.bam, mut.cyt.bam = mn_cyt_vcp.bam, ctrl.cyt.bam = mn_cyt_ctrl.bam,
                                 plot_transcript_number = 6, print2screen = TRUE, save = TRUE, limit.max.ctrl = "NA", limit.min.ctrl = "NA", limit.max.vcp = "NA", limit.min.vcp = "NA"){
  gtf.gene <- gtf.apa %>% filter(gene_name == gene.name) %>% group_by(gene_id, gene_name, strand, seqnames) %>% summarise(start = min(start), end = max(end), width = max(width))
  print(gtf.gene)
  # cat(blue(gtf.gene))
  if (nrow(gtf.gene) == 0){
    stop("ensembl ID not found in GTF")
  }
  gene_name <-  gtf.gene$gene_name
  thechr <-  gtf.gene$seqnames  # without chr before number
  st <- gtf.gene$start
  en <- gtf.gene$end
  strand <- gtf.gene$strand
  utr_length = gtf.gene$width
  # xpad <- utr_length[1]*2.0
  gr <- GRanges(thechr, IRanges(st - utr_length[1]*0.1, en + utr_length[1]*0.1), strand = strand)
  # gr <- GRanges(thechr, IRanges(st - 50, en + 400), strand = strand) # Granges for retained intron with padding either side
  
  CTRL <- importBam(ctrl.nuc.bam, ctrl.cyt.bam, ranges=gr, pairs = TRUE) # import BAMs
  VCP <- importBam(mut.nuc.bam, mut.cyt.bam, ranges=gr, pairs = TRUE)
  CTRL$dat <- coverageGR(CTRL$dat)    # calculate coverage
  VCP$dat <- coverageGR(VCP$dat)
  
  gene <- geneModelFromTxdb(txdb, org.Hs.eg.db, gr = gr, strand = as.character(strand))   # Build Gene model
  # gene[(elementMetadata(gene)[, "feature"] == "utr3")]
  # entrezIDforFMR1 <- BiocGenerics::get("METTL22", org.Hs.egSYMBOL2EG)
  # three_prime_utr <- geneTrack(ids = entrezIDforFMR1, txdb = txdb)
  # 
  # geneTrack(ids = entrezIDforFMR1, txdb = txdb)
  # , symbols = "METTL22", type = "transcript"
  # columns(org.Hs.eg.db)
  # keytypes(org.Hs.eg.db)
  # txby <- transcriptsBy(txdb, by="three_prime_utr")
  # txby
  # threeUTRsByTranscript(txdb)
  
  # Set track arguments
  optSty <- optimizeStyle(trackList(VCP, CTRL, gene))
  trackList <- optSty$tracks
  viewerStyle <- optSty$style  # viewerStyle <- trackViewerStyle()
  setTrackViewerStyleParam(viewerStyle, "xaxis", FALSE)
  setTrackViewerStyleParam(viewerStyle, "margin", c(.01, .10, -0.1, .01))
  
  # optimise coverage plots
  for(i in 1:2) {
    setTrackXscaleParam(trackList[[i]], "gp", list(cex=0.9))
    setTrackStyleParam(trackList[[i]], "height", .36)
    setTrackStyleParam(trackList[[i]], "color", c("firebrick2", "dodgerblue2"))
    setTrackStyleParam(trackList[[i]], "ylabgp", list(cex=0.8))
  }
  if (length(trackList) > (plot_transcript_number + 2) ){
    print(paste("There are ", length(trackList), " ensembl transcripts. Plotting only the first ", plot_transcript_number)) 
    for(i in (plot_transcript_number + 3):length(trackList)) {
      print(trackList[[i]]$name) 
    }
    idx = c((plot_transcript_number + 3):length(trackList))
    trackList <- trackList[-idx] # only plot first selected transcript tracks - remove remaining transcript tracks
  }
  # modify the transcript tracks
  for(i in 3:length(trackList)) {
    setTrackStyleParam(trackList[[i]], "height", .28/(plot_transcript_number+2))
    setTrackStyleParam(trackList[[i]], "ylabgp", list(cex=0))
    setTrackStyleParam(trackList[[i]], "color", "black")
  }
  
  if (limit.max.ctrl != "NA"){ 
    setTrackStyleParam(trackList[[2]], "ylim", c(limit.min.ctrl, limit.max.ctrl))  ## Adjust the limit of CTRL y-axis
  }
  if (limit.max.vcp != "NA"){
    setTrackStyleParam(trackList[[1]], "ylim", c(limit.min.vcp, limit.max.vcp))  ## Adjust the limit of VCP y-axis
  }
  
  if (print2screen == TRUE ){
    print("Printing to screen.")
    vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
    grid.text("Nuclear reads", x=.85, y=.76, just="bottom", gp=gpar(fontsize=6, col = "firebrick2"))
    grid.text("Cytoplasm reads", x=.85, y=.41, just="bottom", gp=gpar(fontsize=6, col = "dodgerblue2"))
    grid.text("Nuclear reads", x=.85, y=.36, just="bottom", gp=gpar(fontsize=6, col = "firebrick2"))
    grid.text("Cytoplasm reads", x=0.85, y=.01, just="bottom", gp=gpar(fontsize=6, col = "dodgerblue2"))
    grid.text(gene_name, x=.06, y=.88, just="bottom", gp=gpar(fontsize=10))#, rot=90)
    dev.off()
  }
  if (save == TRUE ){
    print(paste0("Saving figure as /camp/lab/luscomben/home/users/ziffo/projects/motor-neuron-mislocalisation/apa/coverage/", gene_name, ".nuc_cyt_coverage_plot.png"))
    png(file = paste0("/camp/lab/luscomben/home/users/ziffo/projects/motor-neuron-mislocalisation/apa/coverage/", gene_name, ".nuc_cyt_coverage_plot.png"), height = 3.5, width = 6.5, units ="in", res = 300)
    vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
    grid.text("Nuclear reads", x=.85, y=.76, just="bottom", gp=gpar(fontsize=6, col = "firebrick2"))
    grid.text("Cytoplasm reads", x=.85, y=.41, just="bottom", gp=gpar(fontsize=6, col = "dodgerblue2"))
    grid.text("Nuclear reads", x=.85, y=.36, just="bottom", gp=gpar(fontsize=6, col = "firebrick2"))
    grid.text("Cytoplasm reads", x=0.85, y=.01, just="bottom", gp=gpar(fontsize=6, col = "dodgerblue2"))
    grid.text(gene_name, x=.06, y=.88, just="bottom", gp=gpar(fontsize=10))#, rot=90)
    dev.off()
  }}



track_viewer_nuc_cyt_tx <- function(transcript_id = "ENST00000650847", mut.nuc.bam = ac_nuc_vcp.bam, ctrl.nuc.bam = ac_nuc_ctrl.bam, mut.cyt.bam = ac_cyt_vcp.bam, ctrl.cyt.bam = ac_cyt_ctrl.bam,
                                 print2screen = TRUE, save = TRUE){
  gtf.gene <- gtf %>% filter(transcript == transcript_id) %>% group_by(transcript, ensemblID, gene_name, strand, seqnames) %>% summarise(start = min(start), end = max(end), width = max(width))
  print(gtf.gene)
  if (nrow(gtf.gene) == 0){
    stop("transcript_id not found in GTF")
  }
  gene_name <-  gtf.gene$gene_name
  thechr <-  gtf.gene$seqnames  # without chr before number
  st <- gtf.gene$start
  en <- gtf.gene$end
  strand <- gtf.gene$strand
  gr <- GRanges(thechr, IRanges(st - 100, en + 100), strand = strand) # Granges for retained intron with padding either side
  
  CTRL <- importBam(ctrl.nuc.bam, ctrl.cyt.bam, ranges=gr, pairs = TRUE) # import BAMs
  VCP <- importBam(mut.nuc.bam, mut.cyt.bam, ranges=gr, pairs = TRUE)
  CTRL$dat <- coverageGR(CTRL$dat)    # calculate coverage
  VCP$dat <- coverageGR(VCP$dat)
  
  gene <- geneModelFromTxdb(txdb, org.Hs.eg.db, gr = gr, strand = as.character(strand))   # Build Gene model
  transcript <- gene[[transcript_id]]
  # Set track arguments
  optSty <- optimizeStyle(trackList(VCP, CTRL, transcript))
  trackList <- optSty$tracks
  viewerStyle <- optSty$style  # viewerStyle <- trackViewerStyle()
  setTrackViewerStyleParam(viewerStyle, "xaxis", FALSE)
  setTrackViewerStyleParam(viewerStyle, "margin", c(.01, .10, -0.1, .01))
  
  # optimise coverage plots
  for(i in 1:2) {
    setTrackXscaleParam(trackList[[i]], "gp", list(cex=0.9))
    setTrackStyleParam(trackList[[i]], "height", .4)
    setTrackStyleParam(trackList[[i]], "color", c("firebrick2", "dodgerblue2"))
    setTrackStyleParam(trackList[[i]], "ylabgp", list(cex=0.8))
  }
  # modify the transcript track
  for(i in 3:length(trackList)) {
    setTrackStyleParam(trackList[[i]], "height", .15)
    setTrackStyleParam(trackList[[i]], "ylabgp", list(cex=0))
    setTrackStyleParam(trackList[[i]], "color", "black")
  }
  if (print2screen == TRUE ){
    print("Printing to screen.")
    vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
    grid.text("Nuclear reads", x=.85, y=.76, just="bottom", gp=gpar(fontsize=6, col = "firebrick2"))
    grid.text("Cytoplasm reads", x=.85, y=.41, just="bottom", gp=gpar(fontsize=6, col = "dodgerblue2"))
    grid.text("Nuclear reads", x=.85, y=.36, just="bottom", gp=gpar(fontsize=6, col = "firebrick2"))
    grid.text("Cytoplasm reads", x=0.85, y=.01, just="bottom", gp=gpar(fontsize=6, col = "dodgerblue2"))
    grid.text(gene_name, x=.05, y=.9, just="bottom", gp=gpar(fontsize=12))#, rot=90)
    dev.off()
  }
  if (save == TRUE ){
    print(paste0("Saving figure as /camp/lab/luscomben/home/users/ziffo/projects/motor-neuron-mislocalisation/expression/coverage/", gene_name, ".nuc_cyt_tx_coverage_plot.png"))
    png(file = paste0("/camp/lab/luscomben/home/users/ziffo/projects/motor-neuron-mislocalisation/expression/coverage/", gene_name, ".nuc_cyt_tx_coverage_plot.png"), height = 3.5, width = 6.5, units ="in", res = 300)
    vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
    grid.text("Nuclear reads", x=.85, y=.76, just="bottom", gp=gpar(fontsize=6, col = "firebrick2"))
    grid.text("Cytoplasm reads", x=.85, y=.41, just="bottom", gp=gpar(fontsize=6, col = "dodgerblue2"))
    grid.text("Nuclear reads", x=.85, y=.36, just="bottom", gp=gpar(fontsize=6, col = "firebrick2"))
    grid.text("Cytoplasm reads", x=0.85, y=.01, just="bottom", gp=gpar(fontsize=6, col = "dodgerblue2"))
    grid.text(gene_name, x=.05, y=.9, just="bottom", gp=gpar(fontsize=12))#, rot=90)
    dev.off()
  }}



track_viewer_nuc_cyt_ir <- function(Intron.GeneName.GeneID.Coords = "SAP25/ENSG00000205307/known-exon/7:100573379-100573596:-", IRFinder = irfinder.mn.nuc_vs_cyt.vcp_vs_ctrl, mut.nuc.bam = mn_nuc_vcp.bam, ctrl.nuc.bam = mn_nuc_ctrl.bam, mut.cyt.bam = mn_cyt_vcp.bam, ctrl.cyt.bam = mn_cyt_ctrl.bam, print2screen = TRUE, save = TRUE){
  if (!Intron.GeneName.GeneID.Coords %in% IRFinder$Intron.GeneName.GeneID.Coords){
    stop("Intron.GeneName.GeneID.Coords not found in IRFinder output")
  }

  intron_type = str_split_fixed(Intron.GeneName.GeneID.Coords, "/", 4)[,3] 
  gene_name <-  str_split_fixed(Intron.GeneName.GeneID.Coords, "/", 4)[,1]
  gene_id <-  str_split_fixed(Intron.GeneName.GeneID.Coords, "/", 4)[,2]
  coords <- str_split_fixed(Intron.GeneName.GeneID.Coords, "/",4)[,4]
  thechr <-  str_split_fixed(coords, ":", 3)[,1]  # without chr before number
  st <- as.numeric(str_split_fixed(coords, ":|-", 4)[,2])
  en <- as.numeric(str_split_fixed(coords, ":|-", 4)[,3])
  strand <- str_split_fixed(coords, ":|-", 4)[,4]
  intron_length = GenomicRanges::width(as.GenomicRange(Intron.GeneName.GeneID.Coords))
  xpad <- intron_length[1]*3.8
  gr <- GRanges(thechr, IRanges(st - xpad, en + xpad), strand = strand) # Granges for entire gene (not just the retained intron)
  # gr <- GRanges(thechr, IRanges(st - 100, en + 100), strand = strand) # Granges for retained intron with padding either side
  
  CTRL <- importBam(ctrl.nuc.bam, ctrl.cyt.bam, ranges=gr, pairs = TRUE) # import BAMs
  VCP <- importBam(mut.nuc.bam, mut.cyt.bam, ranges=gr, pairs = TRUE)
  CTRL$dat <- coverageGR(CTRL$dat)    # calculate coverage
  VCP$dat <- coverageGR(VCP$dat)
  
  gene <- geneModelFromTxdb(txdb, org.Hs.eg.db, gr = gr, strand = as.character(strand))   # Build Gene model
  
  # Set track arguments
  optSty <- optimizeStyle(trackList(VCP, CTRL, gene))
  trackList <- optSty$tracks
  viewerStyle <- optSty$style  # viewerStyle <- trackViewerStyle()
  setTrackViewerStyleParam(viewerStyle, "xaxis", FALSE)
  setTrackViewerStyleParam(viewerStyle, "margin", c(.01, .10, -0.1, .01))
  
  # optimise coverage plots
  for(i in 1:2) {
    setTrackXscaleParam(trackList[[i]], "gp", list(cex=0.9))
    setTrackStyleParam(trackList[[i]], "height", .36)
    setTrackStyleParam(trackList[[i]], "color", c("firebrick2", "dodgerblue2"))
    setTrackStyleParam(trackList[[i]], "ylabgp", list(cex=0.8))
  }
  
  if (length(trackList) > 8 ){
    print(paste("There are ", length(trackList) - 4, " ensembl transcripts. Plotting only the first 6.")) 
    for(i in 9:length(trackList)) {
      print(trackList[[i]]$name) 
    }
    idx = c(9:length(trackList))
    trackList <- trackList[-idx] # only plot first 6 transcript tracks - remove remaining transcript tracks
  }
  # modify the transcript tracks
  for(i in 3:length(trackList)) {
    setTrackStyleParam(trackList[[i]], "height", .03)
    setTrackStyleParam(trackList[[i]], "ylabgp", list(cex=0))
    setTrackStyleParam(trackList[[i]], "color", "black")
  }
  if (print2screen == TRUE ){
  print("Printing to screen.")
  vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
  addGuideLine(c(st, en), vp=vp, col = "grey", lwd = 3, lty = "dashed")
  addGuideLine(c(st, en), vp=vp, col = "grey", lwd = 3, lty = "dashed")
  grid.text("Nuclear reads", x=.85, y=.76, just="bottom", gp=gpar(fontsize=6, col = "firebrick2"))
  grid.text("Cytoplasm reads", x=.85, y=.41, just="bottom", gp=gpar(fontsize=6, col = "dodgerblue2"))
  grid.text("Nuclear reads", x=.85, y=.36, just="bottom", gp=gpar(fontsize=6, col = "firebrick2"))
  grid.text("Cytoplasm reads", x=0.85, y=.01, just="bottom", gp=gpar(fontsize=6, col = "dodgerblue2"))
  grid.text(gene_name, x=.06, y=.88, just="bottom", gp=gpar(fontsize=10))#, rot=90)
  }
  if (save == TRUE ){
  print(paste0("Saving figure as /camp/lab/luscomben/home/users/ziffo/projects/motor-neuron-mislocalisation/splicing/IRFinder/coverage/", gene_name, ".", thechr, ".", st, ".", en, ".nuc_cyt_IR_coverage_plot.png"))
  png(file = paste0("/camp/lab/luscomben/home/users/ziffo/projects/motor-neuron-mislocalisation/splicing/IRFinder/coverage/", gene_name, ".", thechr, ".", st, ".", en, ".nuc_cyt_IR_coverage_plot.png"), height = 3.5, width = 6.5, units ="in", res = 300)
  vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
  addGuideLine(c(st, en), vp=vp, col = "grey", lwd = 3, lty = "dashed")
  addGuideLine(c(st, en), vp=vp, col = "grey", lwd = 3, lty = "dashed")
  grid.text("Nuclear reads", x=.85, y=.76, just="bottom", gp=gpar(fontsize=6, col = "firebrick2"))
  grid.text("Cytoplasm reads", x=.85, y=.41, just="bottom", gp=gpar(fontsize=6, col = "dodgerblue2"))
  grid.text("Nuclear reads", x=.85, y=.36, just="bottom", gp=gpar(fontsize=6, col = "firebrick2"))
  grid.text("Cytoplasm reads", x=0.85, y=.01, just="bottom", gp=gpar(fontsize=6, col = "dodgerblue2"))
  grid.text(gene_name, x=.06, y=.88, just="bottom", gp=gpar(fontsize=10))#, rot=90)
  dev.off()
  }
}


track_viewer_nuc_cyt_clip <- function(gene.name = "FMR1", mut.nuc.bam = mn_nuc_vcp.bam, ctrl.nuc.bam = mn_nuc_ctrl.bam, mut.cyt.bam = mn_cyt_vcp.bam, ctrl.cyt.bam = mn_cyt_ctrl.bam,
                                      plot_transcript_number = 5, print2screen = TRUE, save = TRUE, limit.max.ctrl = "NA", limit.min.ctrl = "NA", limit.max.vcp = "NA", limit.min.vcp = "NA"){
  gtf.gene <- gtf %>% filter(gene_name == gene.name) %>% group_by(ensemblID, gene_name, strand, seqnames) %>% summarise(start = min(start), end = max(end), width = max(width))
  print(gtf.gene)
  if (nrow(gtf.gene) == 0){
    stop("ensembl ID not found in GTF")
  }
  gene_name <-  gtf.gene$gene_name
  thechr <-  gtf.gene$seqnames  # without chr before number
  st <- gtf.gene$start
  en <- gtf.gene$end
  strand <- gtf.gene$strand
  gr <- GRanges(thechr, IRanges(st - 100, en + 100), strand = strand) # Granges for retained intron with padding either side
  
  CTRL <- importBam(ctrl.nuc.bam, ctrl.cyt.bam, ranges=gr, pairs = TRUE) # import BAMs
  VCP <- importBam(mut.nuc.bam, mut.cyt.bam, ranges=gr, pairs = TRUE)

  CTRL$dat <- coverageGR(CTRL$dat)    # calculate coverage
  VCP$dat <- coverageGR(VCP$dat)
  
  gene <- geneModelFromTxdb(txdb, org.Hs.eg.db, gr = gr, strand = as.character(strand))   # Build Gene model
  
  # Import iCLIP peaks
  if(file.exists(paste0("/camp/home/ziffo/home/genomes/rbp-databases/clip-bedgraphs/HepG2/", gene.name, ".merged.bed.gz"))){
    print(paste0(gene.name, " clip bed found in HepG2"))
    CLIP <- import(paste0("/camp/home/ziffo/home/genomes/rbp-databases/clip-bedgraphs/HepG2/", gene.name, ".merged.bed.gz"), "BED")
  }else if(file.exists(paste0("/camp/home/ziffo/home/genomes/rbp-databases/clip-bedgraphs/K562/", gene.name, ".merged.bed.gz"))){
    print(paste0(gene.name, " clip bed found in K562"))
    CLIP <- import(paste0("/camp/home/ziffo/home/genomes/rbp-databases/clip-bedgraphs/K562/", gene.name, ".merged.bed.gz"), "BED")
  } else if(file.exists(paste0("/camp/home/ziffo/home/genomes/rbp-databases/clip-bedgraphs/Hela_hg19/", gene.name, ".merged.bed.gz"))){
    print(paste0(gene.name, " clip bed found in Hela_hg19"))
    CLIP <- import(paste0("/camp/home/ziffo/home/genomes/rbp-databases/clip-bedgraphs/Hela_hg19/", gene.name, ".merged.bed.gz"), "BED")
  }
  seqlevelsStyle(CLIP) <- "NCBI"
  CLIP <- new("track", dat=CLIP, type="data", format="BED")

  # Set track arguments
  optSty <- optimizeStyle(trackList(CLIP, VCP, CTRL, gene))
  trackList <- optSty$tracks
  viewerStyle <- optSty$style  # viewerStyle <- trackViewerStyle()
  setTrackViewerStyleParam(viewerStyle, "xaxis", FALSE)
  setTrackViewerStyleParam(viewerStyle, "margin", c(.01, .10, -0.1, .01))
  
  # optimise CLIP track
  setTrackXscaleParam(trackList[[1]], "gp", list(cex=0.9))
  setTrackStyleParam(trackList[[1]], "height", .1)
  setTrackStyleParam(trackList[[1]], "color", "forestgreen")
  setTrackStyleParam(trackList[[1]], "ylabgp", list(cex=0.8))
  setTrackStyleParam(trackList[[1]], "ylim", c(0, 1))  ## Adjust the limit of y-axis
  
  # optimise coverage tracks
  for(i in 2:3) {
    setTrackXscaleParam(trackList[[i]], "gp", list(cex=0.9))
    setTrackStyleParam(trackList[[i]], "height", .35)
    setTrackStyleParam(trackList[[i]], "color", c("firebrick2", "dodgerblue2"))
    setTrackStyleParam(trackList[[i]], "ylabgp", list(cex=0.8))
  }
  
  # modify the transcript tracks
  if (length(trackList) > (plot_transcript_number + 3) ){
    print(paste("There are ", length(trackList), " ensembl transcripts. Plotting only the first ", plot_transcript_number)) 
    for(i in (plot_transcript_number + 4):length(trackList)) {
      print(trackList[[i]]$name) 
    }
    idx = c((plot_transcript_number + 4):length(trackList))
    trackList <- trackList[-idx] # only plot first selected transcript tracks - remove remaining transcript tracks
  }
  for(i in 4:length(trackList)) {
    setTrackStyleParam(trackList[[i]], "height", .2/(plot_transcript_number+2))
    setTrackStyleParam(trackList[[i]], "ylabgp", list(cex=0))
    setTrackStyleParam(trackList[[i]], "color", "black")
  }
  
  if (limit.max.ctrl != "NA"){ 
    setTrackStyleParam(trackList[[2]], "ylim", c(limit.min.ctrl, limit.max.ctrl))  ## Adjust the limit of CTRL y-axis
  }
  if (limit.max.vcp != "NA"){
    setTrackStyleParam(trackList[[1]], "ylim", c(limit.min.vcp, limit.max.vcp))  ## Adjust the limit of VCP y-axis
  }
  
  if (print2screen == TRUE ){
    print("Printing to screen.")
    vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
    grid.text("Nuclear reads", x=.85, y=.80, just="bottom", gp=gpar(fontsize=6, col = "firebrick2"))
    grid.text("Cytoplasm reads", x=.85, y=.51, just="bottom", gp=gpar(fontsize=6, col = "dodgerblue2"))
    grid.text("Nuclear reads", x=.85, y=.45, just="bottom", gp=gpar(fontsize=6, col = "firebrick2"))
    grid.text("Cytoplasm reads", x=0.85, y=.16, just="bottom", gp=gpar(fontsize=6, col = "dodgerblue2"))
    grid.text(gene_name, x=.06, y=.90, just="bottom", gp=gpar(fontsize=10))#, rot=90)
    dev.off()
  }
  if (save == TRUE ){
    print(paste0("Saving figure as /camp/lab/luscomben/home/users/ziffo/projects/motor-neuron-mislocalisation/expression/coverage/", gene_name, ".nuc_cyt_CLIP_coverage_plot.png"))
    png(file = paste0("/camp/lab/luscomben/home/users/ziffo/projects/motor-neuron-mislocalisation/expression/coverage/", gene_name, ".nuc_cyt_CLIP_coverage_plot.png"), height = 3.5, width = 6.5, units ="in", res = 300)
    vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
    grid.text("Nuclear reads", x=.85, y=.80, just="bottom", gp=gpar(fontsize=6, col = "firebrick2"))
    grid.text("Cytoplasm reads", x=.85, y=.51, just="bottom", gp=gpar(fontsize=6, col = "dodgerblue2"))
    grid.text("Nuclear reads", x=.85, y=.45, just="bottom", gp=gpar(fontsize=6, col = "firebrick2"))
    grid.text("Cytoplasm reads", x=0.85, y=.16, just="bottom", gp=gpar(fontsize=6, col = "dodgerblue2"))
    grid.text(gene_name, x=.06, y=.90, just="bottom", gp=gpar(fontsize=10))#, rot=90)
    dev.off()
  }}


track_viewer <- function(gene_test = "SFPQ", txt = ac_who_vcp_vs_ctrl, mut.bam = ac_who_vcp, ctrl.bam = ac_who_ctrl, xpad = 1000){
  top_event <- txt %>% dplyr::filter(p0.05 == "significant") %>% dplyr::filter(gene_name == gene_test) %>% top_n(., 1, abs(IRratio.diff))
  IRratio_ctrl <- top_event %>% dplyr::pull(B.IRratio) # B = ctrl
  IRratio_ctrl <- format(round(IRratio_ctrl,2),nsmall = 2)
  IRratio_vcp <- top_event %>% dplyr::pull(A.IRratio) # A = vcp
  IRratio_vcp <- format(round(IRratio_vcp,2),nsmall = 2)
  # intron coords
  thechr <-  top_event %>% dplyr::pull(Chr) # without chr before number
  st <-  top_event %>% dplyr::pull(Start) # 35184593
  en <-  top_event %>% dplyr::pull(End) # 35187000
  strand <- top_event %>% dplyr::pull(Direction) 
  # whole gene coords
  # gene_test_gtf <- gtf %>% filter(gene_name %in% gene_test)
  # gene_chr <- gene_test_gtf %>% dplyr::pull(seqnames) %>% unique
  # gene_st <- gene_test_gtf %>% dplyr::pull(start) %>% min
  # gene_en <- gene_test_gtf %>% dplyr::pull(end) %>% max
  # gr <- GRanges(gene_chr, IRanges(gene_st, gene_en), strand= strand) # set width of tracks as entire gene
  gr <- GRanges(thechr, IRanges(st - xpad, en + xpad), strand= strand)
  CTRL <- importBam(ctrl.bam, ranges=gr, pairs = TRUE) # import BAM
  VCP <- importBam(mut.bam, ranges=gr, pairs = TRUE) # import BAM
  CTRL$dat <- coverageGR(CTRL$dat) # calculate coverage
  VCP$dat <- coverageGR(VCP$dat) # calculate coverage
  # Build Gene model
  # trs <- geneModelFromTxdb(txdb, org.Hs.eg.db, gr = gr) # build gene model showing transcripts
  # ids <- getGeneIDsFromTxDb(gr, txdb) # get ENSG of gene_test
  ids <- top_event %>% dplyr::pull(gene_id) # get ENSG of gene_test
  gene_track <- geneTrack(ids,txdb)[[1]] # gene all exons on single line
  
  # View the tracks
  optSty <- optimizeStyle(trackList(VCP, CTRL, gene_track))
  trackList <- optSty$tracks
  names(trackList)[3] <- gene_test # rename gene_track with gene_test name - shows as the gene label
  viewerStyle <- optSty$style
  setTrackViewerStyleParam(viewerStyle, "xaxis", FALSE)
  setTrackViewerStyleParam(viewerStyle, "margin", c(.01, .1, -0.1, .01))
  setTrackXscaleParam(trackList[[1]], "draw", TRUE)
  # setTrackXscaleParam(trackList[[1]], "gp", list(cex=0.9))
  setTrackStyleParam(trackList[[1]], "height", .4)
  setTrackStyleParam(trackList[[2]], "height", .4)
  setTrackStyleParam(trackList[[3]], "height", .1)
  setTrackStyleParam(trackList[[3]], "ylabgp", list(cex=1.5))
  setTrackStyleParam(trackList[[1]], "color", c("darkgrey", "black"))
  setTrackStyleParam(trackList[[2]], "color", c("darkgrey", "black"))
  setTrackStyleParam(trackList[[3]], "color", c("dodgerblue2", "black"))
  
  # save figure
  png(file = paste("/camp/home/ziffo/home/projects/astrocyte-ir-als/splicing/IRFinder/figures/coverage/", gene_test, "_coverage_plot.png", sep = ""), height = 3.5, width = 10, units ="in", res = 300)
  vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
  addGuideLine(c(st, en), vp=vp, col = "firebrick2", lwd = 3, lty = "dashed")
  grid.text(IRratio_ctrl, x=.545, y=.6, just="bottom")
  grid.text(IRratio_vcp, x=.545, y=.2, just="bottom")
  dev.off()
}

# track_viewer(gene_test = "SFPQ", txt = ac_who_vcp_vs_ctrl, mut.bam = ac_who_vcp, ctrl.bam = ac_who_ctrl, xpad = 7000)

track_viewer_tx <- function(gene_test = "RECK", txt = ac_who_vcp_vs_ctrl, mut.bam = ac_who_vcp, ctrl.bam = ac_who_ctrl, path = NULL){
  txt %>% dplyr::filter(p0.05 == "significant" & gene_name == gene_test) %>% print
  top_event <- txt %>% dplyr::filter(p0.05 == "significant" & gene_name == gene_test) %>% top_n(., 1, abs(IRratio.diff))
  txt_nrow <- txt  %>% dplyr::filter(p0.05_reliable == "significant" & gene_name == gene_test) %>% nrow
  if(txt_nrow > 1){
    print(paste(gene_test, " has ", txt_nrow, " significant reliable IR events. Only the top deltaIR event will be plotted.", sep = ""))
  }
  IRratio_ctrl <- top_event %>% dplyr::pull(B.IRratio) # B = ctrl
  IRratio_ctrl <- format(round(IRratio_ctrl,2),nsmall = 2)
  IRratio_vcp <- top_event %>% dplyr::pull(A.IRratio) # A = vcp
  IRratio_vcp <- format(round(IRratio_vcp,2),nsmall = 2)
  thechr <-  top_event %>% dplyr::pull(Chr) # without chr before number
  st <-  top_event %>% dplyr::pull(Start) # 35184593
  en <-  top_event %>% dplyr::pull(End) # 35187000
  strand <- top_event %>% dplyr::pull(Direction) 
  intron_length <- top_event %>% dplyr::pull(intron_length)
  xpad <- intron_length[1]*3.8
  gr <- GRanges(thechr, IRanges(st - xpad, en + xpad), strand = strand) # Granges for entire gene (not just the retained intron)
  CTRL <- importBam(ctrl.bam, ranges=gr, pairs = TRUE) # import BAM
  VCP <- importBam(mut.bam, ranges=gr, pairs = TRUE) # import BAM
  CTRL$dat <- coverageGR(CTRL$dat) # calculate coverage
  VCP$dat <- coverageGR(VCP$dat) # calculate coverage
  # Build Gene model
  trs <- geneModelFromTxdb(txdb, org.Hs.eg.db, gr = gr) # build gene model showing transcripts
  # ids <- getGeneIDsFromTxDb(gr, txdb) # get ENSG of gene_test
  ids <- top_event %>% dplyr::pull(gene_id) # get ENSG of gene_test
  # gene_track <- geneTrack(ids,txdb)[[1]] # gene all exons on single line
  
  print("Creating the coverage and gene model tracks.")
  optSty <- optimizeStyle(trackList(VCP, CTRL, trs))
  trackList <- optSty$tracks
  viewerStyle <- optSty$style
  setTrackViewerStyleParam(viewerStyle, "xaxis", FALSE)
  setTrackViewerStyleParam(viewerStyle, "margin", c(.01, .15, -0.1, .01))
  setTrackXscaleParam(trackList[[1]], "draw", TRUE)
  # setTrackXscaleParam(trackList[[1]], "gp", list(cex=0.9))
  setTrackStyleParam(trackList[[1]], "height", .36)
  setTrackStyleParam(trackList[[2]], "height", .36)
  setTrackXscaleParam(trackList[[1]], "draw", TRUE)
  setTrackXscaleParam(trackList[[1]], "gp", list(cex=0.9))
  setTrackStyleParam(trackList[[1]], "color", c("darkgrey", "black"))
  setTrackStyleParam(trackList[[2]], "color", c("darkgrey", "black"))
  
  if (length(trackList) > 8 ){
    print(paste("There are ", length(trackList), " ensembl transcripts. Plotting only the first 6.")) 
    for(i in 9:length(trackList)) {
      print(trackList[[i]]$name) 
    }
    idx = c(9:length(trackList))
    trackList <- trackList[-idx] # only plot first 6 transcript tracks - remove remaining transcript tracks
  }
  # modify the transcript tracks
  for(i in 3:length(trackList)) {
    setTrackStyleParam(trackList[[i]], "height", .03)
    setTrackStyleParam(trackList[[i]], "ylabgp", list(cex=0))
    setTrackStyleParam(trackList[[i]], "color", c("dodgerblue2", "black"))
  }
  
  vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
  addGuideLine(c(st, en), vp=vp, col = "firebrick2", lwd = 3, lty = "dashed")
  addGuideLine(c(st, en), vp=vp, col = "firebrick2", lwd = 3, lty = "dashed")
  grid.text(IRratio_ctrl, x=.57, y=.6, just="bottom")
  grid.text(IRratio_vcp, x=.57, y=.2, just="bottom")
  grid.text(gene_test, x=.09, y=.87, just="bottom", gp=gpar(cex=1.5))
  
  print("Saving figure.")
  png(file = paste("/camp/home/ziffo/home/projects/astrocyte-ir-als/splicing/IRFinder/figures/coverage/", path, gene_test, "_tx_coverage_plot.png", sep = ""), height = 2.2, width = 7, units ="in", res = 300)
  vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
  addGuideLine(c(st, en), vp=vp, col = "firebrick2", lwd = 2, lty = "dashed")
  grid.text(IRratio_ctrl, x=.57, y=.65, just="bottom")
  grid.text(IRratio_vcp, x=.57, y=.2, just="bottom")
  grid.text(gene_test, x=.09, y=.87, just="bottom", gp=gpar(cex=1.5))
  dev.off()
}

# plot top delta IR ratio event irrespective of p.diff from IRFinder
track_viewer_tx_top_IRratio_event <- function(gene_test = "RECK", txt = ac_who_vcp_vs_ctrl, mut.bam = ac_who_vcp, ctrl.bam = ac_who_ctrl, path = NULL){
  txt %>% dplyr::filter(gene_name == gene_test) %>% print
  top_event <- txt %>% dplyr::filter(gene_name == gene_test) %>% top_n(., 1, abs(IRratio.diff))
  txt_nrow <- txt  %>% dplyr::filter(p0.05_reliable == "significant" & gene_name == gene_test) %>% nrow
  cat(blue(gene_test, " has ", txt_nrow, " significant reliable IR events", sep = ""))
  txt_all_nrow <- txt  %>% dplyr::filter(gene_name == gene_test) %>% nrow
  if (txt_all_nrow == 0){
    stop("Gene has no events in IRFinder differential output")
  }
  IRratio_ctrl <- top_event %>% dplyr::pull(B.IRratio) # B = ctrl
  IRratio_ctrl <- format(round(IRratio_ctrl,2),nsmall = 2)
  IRratio_vcp <- top_event %>% dplyr::pull(A.IRratio) # A = vcp
  IRratio_vcp <- format(round(IRratio_vcp,2),nsmall = 2)
  thechr <-  top_event %>% dplyr::pull(Chr) # without chr before number
  st <-  top_event %>% dplyr::pull(Start) # 35184593
  en <-  top_event %>% dplyr::pull(End) # 35187000
  strand <- top_event %>% dplyr::pull(Direction) 
  intron_length <- top_event %>% dplyr::pull(intron_length)
  xpad <- intron_length[1]*3.8
  gr <- GRanges(thechr, IRanges(st - xpad, en + xpad), strand = strand) # Granges for entire gene (not just the retained intron)
  CTRL <- importBam(ctrl.bam, ranges=gr, pairs = TRUE) # import BAM
  VCP <- importBam(mut.bam, ranges=gr, pairs = TRUE) # import BAM
  CTRL$dat <- coverageGR(CTRL$dat) # calculate coverage
  VCP$dat <- coverageGR(VCP$dat) # calculate coverage
  # Build Gene model
  trs <- geneModelFromTxdb(txdb, org.Hs.eg.db, gr = gr) # build gene model showing transcripts
  # ids <- getGeneIDsFromTxDb(gr, txdb) # get ENSG of gene_test
  ids <- top_event %>% dplyr::pull(gene_id) # get ENSG of gene_test
  # gene_track <- geneTrack(ids,txdb)[[1]] # gene all exons on single line
  
  print("Creating the coverage and gene model tracks.")
  optSty <- optimizeStyle(trackList(VCP, CTRL, trs))
  trackList <- optSty$tracks
  viewerStyle <- optSty$style
  setTrackViewerStyleParam(viewerStyle, "xaxis", FALSE)
  setTrackViewerStyleParam(viewerStyle, "margin", c(.01, .15, -0.1, .01))
  setTrackXscaleParam(trackList[[1]], "draw", TRUE)
  # setTrackXscaleParam(trackList[[1]], "gp", list(cex=0.9))
  setTrackStyleParam(trackList[[1]], "height", .36)
  setTrackStyleParam(trackList[[2]], "height", .36)
  setTrackXscaleParam(trackList[[1]], "draw", TRUE)
  setTrackXscaleParam(trackList[[1]], "gp", list(cex=0.9))
  setTrackStyleParam(trackList[[1]], "color", c("darkgrey", "black"))
  setTrackStyleParam(trackList[[2]], "color", c("darkgrey", "black"))
  
  if (length(trackList) > 8 ){
    print(paste("There are ", length(trackList), " ensembl transcripts. Plotting only the first 6.")) 
    for(i in 9:length(trackList)) {
      print(trackList[[i]]$name) 
    }
    idx = c(9:length(trackList))
    trackList <- trackList[-idx] # only plot first 6 transcript tracks - remove remaining transcript tracks
  }
  # modify the transcript tracks
  for(i in 3:length(trackList)) {
    setTrackStyleParam(trackList[[i]], "height", .03)
    setTrackStyleParam(trackList[[i]], "ylabgp", list(cex=0))
    setTrackStyleParam(trackList[[i]], "color", c("dodgerblue2", "black"))
  }
  
  vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
  addGuideLine(c(st, en), vp=vp, col = "firebrick2", lwd = 3, lty = "dashed")
  addGuideLine(c(st, en), vp=vp, col = "firebrick2", lwd = 3, lty = "dashed")
  grid.text(IRratio_ctrl, x=.57, y=.6, just="bottom")
  grid.text(IRratio_vcp, x=.57, y=.2, just="bottom")
  grid.text(gene_test, x=.09, y=.87, just="bottom", gp=gpar(cex=1.5))
  
  print("Saving figure.")
  png(file = paste("/camp/home/ziffo/home/projects/astrocyte-ir-als/splicing/IRFinder/figures/coverage/", path, gene_test, "_topDeltaIRratio_tx_coverage_plot.png", sep = ""), height = 2.2, width = 7, units ="in", res = 300)
  vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
  addGuideLine(c(st, en), vp=vp, col = "firebrick2", lwd = 2, lty = "dashed")
  grid.text(IRratio_ctrl, x=.57, y=.65, just="bottom")
  grid.text(IRratio_vcp, x=.57, y=.2, just="bottom")
  grid.text(gene_test, x=.09, y=.87, just="bottom", gp=gpar(cex=1.5))
  dev.off()
}


# Plot multiple IR events on a single gene coverage plot
track_viewer_tx_multiple <- function(gene_test = "SFPQ", txt = ac_who_vcp_vs_ctrl, mut.bam = ac_who_vcp, ctrl.bam = ac_who_ctrl, path = NULL){
  txt_nrow <- txt  %>% dplyr::filter(p0.05_reliable == "significant" & gene_name == gene_test) %>% nrow
  if(txt_nrow > 1){
    print(paste(gene_test, " has ", txt_nrow, " significant reliable IR events.", sep = ""))
  }
  events <- txt %>% dplyr::filter(p0.05_reliable == "significant" & gene_name == gene_test) 
  
  # set GRanges for entire gene
  # whole gene coords
  gene_test_gtf <- gtf %>% filter(gene_name %in% gene_test)
  gene_chr <- gene_test_gtf %>% dplyr::pull(seqnames) %>% unique %>% as.numeric
  gene_st <- gene_test_gtf %>% dplyr::pull(start) %>% min %>% as.numeric
  gene_en <- gene_test_gtf %>% dplyr::pull(end) %>% max %>% as.numeric
  strand <- gene_test_gtf %>% dplyr::pull(strand) %>% unique %>% as.character
  gr <- GRanges(gene_chr, IRanges(gene_st, gene_en), strand= strand) # set width of tracks as entire gene
  CTRL <- importBam(ctrl.bam, ranges=gr, pairs = TRUE) # import BAM
  VCP <- importBam(mut.bam, ranges=gr, pairs = TRUE) # import BAM
  CTRL$dat <- coverageGR(CTRL$dat) # calculate coverage
  VCP$dat <- coverageGR(VCP$dat) # calculate coverage
  trs <- geneModelFromTxdb(txdb, org.Hs.eg.db, gr = gr) # build gene model showing transcripts
  
  print("Creating the coverage and gene model tracks.")
  optSty <- optimizeStyle(trackList(VCP, CTRL, trs))
  trackList <- optSty$tracks
  viewerStyle <- optSty$style
  # setTrackViewerStyleParam(viewerStyle, "xaxis", FALSE)
  # setTrackViewerStyleParam(viewerStyle, "margin", c(.01, .1, 0.03, .01))
  # setTrackXscaleParam(trackList[[1]], "draw", TRUE)
  # setTrackXscaleParam(trackList[[1]], "gp", list(cex=0.9))
  setTrackViewerStyleParam(viewerStyle, "xaxis", FALSE)
  setTrackViewerStyleParam(viewerStyle, "margin", c(.01, .15, -0.1, .01))
  setTrackStyleParam(trackList[[1]], "height", .36)
  setTrackStyleParam(trackList[[2]], "height", .36)
  setTrackXscaleParam(trackList[[1]], "draw", TRUE)
  setTrackXscaleParam(trackList[[1]], "gp", list(cex=0.9))
  setTrackStyleParam(trackList[[1]], "color", c("darkgrey", "black"))
  setTrackStyleParam(trackList[[2]], "color", c("darkgrey", "black"))
  
  # only show first 6 transcript tracks - remove remaining transcript tracks
  if (length(trackList) > 8 ){
    print("removing transcript tracks > 6:")
    for(i in 9:length(trackList)) {
      print(trackList[[i]]$name)
    }
    idx = c(9:length(trackList))
    trackList <- trackList[-idx]
    # [[9:length(trackList)]] <- NULL
  }
  # modify the transcript tracks
  for(i in 3:length(trackList)) {
    setTrackStyleParam(trackList[[i]], "height", .03)
    setTrackStyleParam(trackList[[i]], "ylabgp", list(cex=0))
    setTrackStyleParam(trackList[[i]], "color", c("dodgerblue2", "black"))
  }
  
  print("Extracting each IR event coords.")
  event <- list()
  IRratio_ctrl <- list()
  IRratio_vcp <- list()
  thechr <- list()
  st <- list()
  en <- list()
  strand <- list()
  for (i in 1:nrow(events)) {
    event[[i]] <- events[i, ] # iterate through rows of events df  
    IRratio_ctrl[[i]] <- event[[i]] %>% dplyr::pull(B.IRratio) # B = ctrl
    IRratio_ctrl[[i]] <- format(round(IRratio_ctrl[[i]],2),nsmall = 2)
    IRratio_vcp[[i]] <- event[[i]] %>% dplyr::pull(A.IRratio) # A = vcp
    IRratio_vcp[[i]] <- format(round(IRratio_vcp[[i]],2),nsmall = 2)
    thechr[[i]] <-  event[[i]] %>% dplyr::pull(Chr) # without chr before number
    st[[i]] <-  event[[i]] %>% dplyr::pull(Start) # 35184593
    en[[i]] <-  event[[i]] %>% dplyr::pull(End) # 35187000
    strand[[i]] <- event[[i]] %>% dplyr::pull(Direction) 
  }
  
  vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
  for (i in 1:nrow(events)) {
    addGuideLine(c(st[[i]], en[[i]]), vp=vp, col = "firebrick2", lwd = 1, lty = "dashed")
    # grid.text(IRratio_ctrl[[i]], x=.545, y=.6, just="bottom")
    # grid.text(IRratio_vcp[[i]], x=.545, y=.2, just="bottom")
  }
  grid.text(gene_test, x=.09, y=.87, just="bottom", gp=gpar(cex=1.5))
  
  print("Saving figure.")
  png(file = paste("/camp/home/ziffo/home/projects/astrocyte-ir-als/splicing/IRFinder/figures/coverage/", path, gene_test, "_tx_multiple_coverage_plot.png", sep = ""), height = 2.5, width = 8, units ="in", res = 300)
  vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
  for (i in 1:nrow(events)) {
    addGuideLine(c(st[[i]], en[[i]]), vp=vp, col = "firebrick2", lwd = 1, lty = "dashed")
    # grid.text(IRratio_ctrl[[i]], x=.545, y=.6, just="bottom")
    # grid.text(IRratio_vcp[[i]], x=.545, y=.2, just="bottom")
  }
  grid.text(gene_test, x=.09, y=.87, just="bottom", gp=gpar(cex=1.5))
  dev.off()
}


track_viewer_fractions <- function(gene_test = "PRPF4", txt = irfinder.ac_vcp_vs_ctrl.filt, mut.nuc.bam = ac_nuc_vcp.bam, ctrl.nuc.bam = ac_nuc_ctrl.bam, mut.cyt.bam = ac_cyt_vcp.bam, ctrl.cyt.bam = ac_cyt_ctrl.bam, path = NULL){
  txt %>% drop_na %>% arrange(IRratio.diff.nuc) %>% dplyr::filter(gene_name == gene_test) %>% print
  top_event <- txt %>% drop_na %>% arrange(IRratio.diff.nuc) %>% dplyr::filter(gene_name == gene_test) %>% top_n(., 1, abs(IRratio.diff.nuc))
  txt_nrow <- txt %>% arrange(IRratio.diff.nuc) %>% dplyr::filter(gene_name == gene_test) %>% nrow
  if(txt_nrow > 1){
    print(paste(gene_test, " has ", txt_nrow, " significant reliable nuclear IR events. Only the top deltaIR event will be plotted.", sep = ""))
  }
  IRratio_nuc_ctrl <- top_event %>% dplyr::pull(IRratio.ctrl.nuc) # B = ctrl
  IRratio_nuc_ctrl <- format(round(IRratio_nuc_ctrl,2),nsmall = 2)
  IRratio_nuc_vcp <- top_event %>% dplyr::pull(IRratio.vcp.nuc) # A = vcp
  IRratio_nuc_vcp <- format(round(IRratio_nuc_vcp,2),nsmall = 2)
  IRratio_cyt_ctrl <- top_event %>% dplyr::pull(IRratio.ctrl.cyt) # B = ctrl
  IRratio_cyt_ctrl <- format(round(IRratio_cyt_ctrl,2),nsmall = 2)
  IRratio_cyt_vcp <- top_event %>% dplyr::pull(IRratio.vcp.cyt) # A = vcp
  IRratio_cyt_vcp <- format(round(IRratio_cyt_vcp,2),nsmall = 2)
  thechr <-  top_event %>% dplyr::pull(Chr) # without chr before number
  st <-  top_event %>% dplyr::pull(Start) # 35184593
  en <-  top_event %>% dplyr::pull(End) # 35187000
  strand <- top_event %>% dplyr::pull(Direction) 
  intron_length <- top_event %>% dplyr::pull(intron_length)
  xpad <- intron_length[1]*3.8
  gr <- GRanges(thechr, IRanges(st - xpad, en + xpad), strand= strand) # Granges for retained intron with padding either side
  nuc.CTRL <- importBam(ctrl.nuc.bam, ranges=gr, pairs = TRUE) # import BAMs
  nuc.VCP <- importBam(mut.nuc.bam, ranges=gr, pairs = TRUE)
  cyt.CTRL <- importBam(ctrl.cyt.bam, ranges=gr, pairs = TRUE)
  cyt.VCP <- importBam(mut.cyt.bam, ranges=gr, pairs = TRUE)
  nuc.CTRL$dat <- coverageGR(nuc.CTRL$dat)    # calculate coverage
  nuc.VCP$dat <- coverageGR(nuc.VCP$dat)
  cyt.CTRL$dat <- coverageGR(cyt.CTRL$dat)
  cyt.VCP$dat <- coverageGR(cyt.VCP$dat)
  # Build Gene model
  trs <- geneModelFromTxdb(txdb, org.Hs.eg.db, gr = gr) # build gene model showing transcripts
  # ids <- getGeneIDsFromTxDb(gr, txdb) # get ENSG of gene_test
  ids <- top_event %>% dplyr::pull(gene_id) # get ENSG of gene_test
  # gene_track <- geneTrack(ids,txdb)[[1]] # gene all exons on single line
  
  print("Creating the coverage and gene model tracks.")
  optSty <- optimizeStyle(trackList(cyt.VCP, cyt.CTRL, nuc.VCP, nuc.CTRL, trs))
  trackList <- optSty$tracks
  viewerStyle <- optSty$style
  setTrackViewerStyleParam(viewerStyle, "xaxis", FALSE)
  setTrackViewerStyleParam(viewerStyle, "margin", c(.01, .10, -0.1, .01))
  # optimise coverage plots
  for(i in 1:4) {
    # setTrackStyleParam(trackList[[i]], "draw", TRUE)
    setTrackXscaleParam(trackList[[i]], "gp", list(cex=0.9))
    setTrackStyleParam(trackList[[i]], "height", .18)
    setTrackStyleParam(trackList[[i]], "color", c("darkgrey", "black"))
    setTrackStyleParam(trackList[[i]], "ylabgp", list(cex=0.8))
  }
  # optimise the transcript tracks
  if (length(trackList) > 10 ){
    print(paste("There are ", length(trackList), " ensembl transcripts. Plotting only the first 6.")) 
    for(i in 11:length(trackList)) {
      print(trackList[[i]]$name) 
    }
    idx = c(11:length(trackList))
    trackList <- trackList[-idx] # only plot first 6 transcript tracks - remove remaining transcript tracks
  }
  for(i in 5:length(trackList)) {
    setTrackStyleParam(trackList[[i]], "height", .03)
    setTrackStyleParam(trackList[[i]], "ylabgp", list(cex=0))
    setTrackStyleParam(trackList[[i]], "color", c("dodgerblue2", "black"))
  }
  
  vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
  addGuideLine(c(st, en), vp=vp, col = "firebrick2", lwd = 3, lty = "dashed")
  grid.text(IRratio_cyt_vcp, x=.545, y=.15, just="bottom")
  grid.text(IRratio_cyt_ctrl, x=.545, y=.35, just="bottom")
  grid.text(IRratio_nuc_vcp, x=.545, y=.55, just="bottom")
  grid.text(IRratio_nuc_ctrl, x=.545, y=.75, just="bottom")
  grid.text(gene_test, x=.05, y=.87, just="bottom", gp=gpar(cex=1.0))
  
  print("Saving figure.")
  png(file = paste("/camp/lab/luscomben/home/users/ziffo/projects/motor-neuron-mislocalisation/splicing/IRFinder/figures/coverage/", path, gene_test, "_tx_coverage_plot.png", sep = ""), height = 4, width = 7, units ="in", res = 300)
  vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
  addGuideLine(c(st, en), vp=vp, col = "firebrick2", lwd = 3, lty = "dashed")
  grid.text(IRratio_cyt_vcp, x=.545, y=.15, just="bottom")
  grid.text(IRratio_cyt_ctrl, x=.545, y=.35, just="bottom")
  grid.text(IRratio_nuc_vcp, x=.545, y=.55, just="bottom")
  grid.text(IRratio_nuc_ctrl, x=.545, y=.75, just="bottom")
  grid.text(gene_test, x=.06, y=.87, just="bottom", gp=gpar(cex=1.0))
  dev.off()
}

# track_viewer_fractions(gene_test = "RND3", txt = ac_vcp_vs_ctrl.filt, mut.nuc.bam = ac_nuc_vcp, ctrl.nuc.bam = ac_nuc_ctrl, mut.cyt.bam = ac_cyt_vcp, ctrl.cyt.bam = ac_cyt_ctrl, xpad = 700)


