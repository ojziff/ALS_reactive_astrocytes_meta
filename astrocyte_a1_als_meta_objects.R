# load("/camp/home/ziffo/home/projects/astrocyte-a1-als-meta/expression/deseq2/astrocyte_a1_als_meta.RData")
mass_spec.res <- readRDS("/camp/lab/luscomben/home/users/ziffo/projects/astrocyte-a1-als-meta/expression/mass-spec/mass_spec.res_vcp_vs_ctrl.rds")
source("/camp/lab/luscomben/home/users/ziffo/projects/astrocyte-a1-als-meta/scripts/ALS_reactive_astrocytes_meta/astrocyte_a1_als_meta_functions.R")


# Metadata ----------------------------------------------------------------
dir <- "/camp/home/ziffo/home/projects/astrocyte-a1-als-meta/alignment/salmon"

### hiPSC datasets

# Patani VCP astrocytes alone
patani.metadata <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/astrocyte-fractionation-bulk-rnaseq/sample-details/sample_table.csv") %>% filter(fraction == "whole", cellline != "ctrl2") %>% distinct(fraction_cellline, .keep_all = TRUE) %>% # 
  mutate(replicate = c(1,2,1,2), sample = cellline, file_salmon = file.path(dir, paste0("ac_who_", vcp, "_R", replicate), "quant.sf"), dataset = "patani", 
         group = gsub("vcp", "als", vcp), group = factor(group, levels = c("ctrl", "als")), region = "ipsc", celltype_source = "This_study", rpt = 1:4, name = paste("astrocyte_ipsc", rpt, sep = "_"), irfinder = paste0("ac_who_",sample)) %>% 
  select(group, replicate, sample, dataset, file_salmon, celltype_source, region, name, irfinder)

# Tyzack SOD1 astrocytes alone
tyzack.metadata <- read_excel("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/astrocyte-sod1-ipsc-tyzack-2017/sample-details/metadata.xlsx") %>%
  mutate(replicate = c(1,2,3,4,1,2,3), sample = paste("ac_tyz_", condition, "_R", replicate, sep = ""), file_salmon = file.path(dir, paste0("ac_tyz_", condition, "_R", replicate), "quant.sf"), dataset = "tyzack", 
         group = gsub("sod1", "als", condition), group = factor(group, levels = c("ctrl", "als")), rpt = 1:7, name = paste0("ac_tyz_astrocyte_ipsc_R", rpt), region = "ipsc", celltype_source = "Tyzack_ipsc", irfinder = Experiment) %>% 
  select(group, replicate, sample, dataset, file_salmon, celltype_source, region, name, irfinder)

# Birger C9orf72 alone
birger.metadata <- read_excel("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/astrocyte-c9orf72-ipsc-birger-2019/sample-details/metadata.xlsx") %>%
  mutate(replicate = c(1,2,1,2), condition = gsub("C9orf72", "c9orf72", condition), condition = gsub("control", "ctrl", condition), file_salmon = file.path(dir, paste0("ac_bir_", condition, "_R", replicate), "quant.sf"),
         group =  gsub("c9orf72", "als", condition), group = factor(group, levels = c("ctrl", "als")), dataset = "birger", rpt = c(1,2,3,4), region = "ipsc", celltype_source = "Birger_ipsc", name = paste0("ac_bir_astrocyte_ipsc_R", rpt), 
         irfinder = paste0("ac_bir_",sample), irfinder = gsub("C9orf72_L", "c9orf", irfinder), irfinder = gsub("control_L","ctrl",irfinder)) %>% 
  select(group, replicate, sample, dataset, file_salmon, celltype_source, region, name,irfinder)

# Neyrinck FUS astrocytes alone - both day 23 & 48 since both valid astrocytes in original paper
neyrinck.metadata <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/astrocyte-fus-ipsc-neyrinck-2021/sample-details/samplesheet.csv") %>% filter(celltype == "ipsc", overexpressed == "sox9oe", day == 48) %>% # remove primary astrocytes & no TF oe
  mutate(file_salmon = file.path(dir, paste0(group, "_R", replicate), "quant.sf"), sample = paste("ac_ney_", group, "_R", replicate, sep = ""),
         condition = group, group = gsub("wt", "ctrl", mutation), group =  gsub("fus", "als", group), group = factor(group, levels = c("ctrl", "als")),
         dataset = "neyrinck", rpt = c(1:12), region = case_when(grepl("ipsc", condition) ~ "ipsc", TRUE ~ "cortex"), celltype_source = paste0("Neyrinck_", region), name = paste0("ac_ney_astrocyte_", region, "_R", rpt), irfinder = run_accession) %>% 
  select(group, replicate, sample, dataset, file_salmon, celltype_source, region, name,irfinder)

# ALS hiPSC meta-analysis: datasets VCP, SOD1 and C9orf72
als_datasets.metadata <- bind_rows(patani.metadata, birger.metadata, tyzack.metadata, neyrinck.metadata) %>% mutate(group = factor(group, levels = c("ctrl", "als")), dataset = as.factor(dataset))

# Barbar cytokine-stimulated A1 astrocytes
barbar.metadata <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/reactive-astrocyte-barbar-2020/sample-details/samplesheet.csv") %>% 
  mutate(file_salmon = file.path(dir, paste(group, "_R", replicate, sep = ""), "quant.sf"), reactivity = c("a0","a1","a0","a1","a0","a1"), group = factor(reactivity, levels = c("a0", "a1")), rpt = c(1:6), region = "ipsc", 
         celltype_source = "Barbar_ipsc", name = paste("CD49_astrocyte_ipsc", rpt, sep = "_"), dataset = "barbar", irfinder = sample, cellline = paste0(dataset,"_",replicate)) %>% 
  select(group, replicate, sample, dataset, file_salmon, celltype_source, region, name, cellline, irfinder)

# Leng cytokine-stimulated A1 astrocytes - 5 protocols - remove Fernandopulle as odd clustering
leng.metadata <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/reactive-astrocyte-ipsc-leng-2021/sample-details/samplesheet.csv") %>% 
  filter(protocol != "Fernandopulle") %>% # remove protocol Fernandopulle
  mutate(sample = paste0(group, "_R", replicate), file_salmon = file.path(dir, sample, "quant.sf"), rpt = c(rep(1:6, each = 4)), region = "ipsc", celltype_source = paste0(protocol,"_",region), name = paste0(protocol,"_ipsc_",rpt), group = factor(treatment, levels = c("a0", "a1")), irfinder = run_accession) %>% 
  select(group, replicate, sample, dataset = protocol, file_salmon, celltype_source, region, name, cellline, irfinder)

# a1 meta-analysis
a1_datasets.metadata <- bind_rows(barbar.metadata, leng.metadata) %>% mutate(group = factor(group, levels = c("a0", "a1")), cellline = as.factor(cellline))

# Patani & Tyzack & Birger & Zhang & Barbar & Bradley together - all samples for PCA
zhang.metadata <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/purified-cortex-human-zhang-2016/sample-details/samplesheet.csv") %>% 
  filter(!group %in% c("astrocyte_mouse_cortex", "astrocyte_adult_tumour", "astrocyte_adult_tumour", "astrocyte_adult_hippocampus", "endothelial_adult_temporal_lobe", "whole_adult_temporal_lobe")) %>% 
  mutate(file_salmon = file.path(dir, paste(group, "_R", replicate, sep = ""), "quant.sf"), celltype_source = paste(celltype,source, sep = "_"), region = gsub("temporal_lobe", "cortex", region), name = sample, dataset = "zhang") %>% 
  select(group, replicate, sample, dataset, file_salmon, region, celltype_source, name) # celltype, 
bradley.metadata <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/regional-astrocytes-ipsc-bradley-2019/sample-details/samplesheet.csv") %>%
  mutate(file_salmon = file.path(dir, paste(group, "_R", replicate, sep = ""), "quant.sf"), sample = paste(group, "_R", replicate, sep = ""), region = "ipsc", celltype_source = rep(c("Bradley_FB_ipsc", "Bradley_FB_ipsc", "Bradley_SC_ipsc", "Bradley_SC_ipsc"), each = 6), celltype= "astrocyte", group = "astrocyte_ipsc", name = sample, dataset = "bradley") %>%  select(group, replicate, sample, dataset, file_salmon, region, celltype_source, name) %>% distinct(sample, .keep_all = TRUE)

zhang_bradley_patani_tyzack_birger_neyrinck.metadata <- bind_rows(zhang.metadata, patani.metadata, tyzack.metadata, neyrinck.metadata, birger.metadata, bradley.metadata) %>% 
  mutate(celltype_source = factor(celltype_source, levels = c("This_study", "Tyzack_ipsc", "Birger_ipsc", "Neyrinck_ipsc", "Bradley_SC_ipsc","Bradley_FB_ipsc", "astrocyte_fetal", "astrocyte_adult", "oligodendrocyte_adult", "neuron_adult", "myeloid_adult", "endothelial_adult", "whole_adult"))) # "Neyrinck_ipsc.SOX9", "Neyrinck_fetal",

# CTRLs from this study & Zhang & Bradley together for characterisation
ctrls.zhang_bradley_patani.metadata <- bind_rows(zhang.metadata, filter(patani.metadata, group == "ctrl"), filter(bradley.metadata, grepl("ventral_spinalcord", name))) %>%
  mutate(celltype_source = factor(celltype_source, levels = c("This_study", "Bradley_VSC_ipsc", "astrocyte_fetal", "astrocyte_adult", "oligodendrocyte_adult", "neuron_adult", "myeloid_adult")))

### Mouse datasets

# Guttenplan mouse cytokine-stimulated A1 astrocytes
guttenplan.A1.metadata <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/reactive-astrocyte-mouse-guttenplan-2020/sample-details/samplesheet.csv") %>% 
  filter(celltype == "astrocyte", mutation == "wt") %>%
  mutate(sample = paste(group, "_R", replicate, sep = ""), file_salmon = file.path(dir, sample, "quant.sf"), reactivity = c("A1","A1","A1","A0","A0","A0"), group = factor(reactivity, levels = c("A0", "A1")), dataset = "guttenplan", irfinder = run_accession) %>% 
  select(group, replicate, sample, dataset, file_salmon, irfinder)

# Guttenplan mouse SOD1 mutant astrocytes
guttenplan.sod1.metadata <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/reactive-astrocyte-mouse-guttenplan-2020/sample-details/samplesheet.csv") %>% 
  filter(celltype == "astrocyte", treatment == "untx") %>%
  mutate(sample = paste(group, "_R", replicate, sep = ""), file_salmon = file.path(dir, sample, "quant.sf"), group = factor(mutation, levels = c("wt", "sod1")), condition = case_when(group == "sod1" ~ "als", group == "wt" ~ "ctrl"), dataset = "guttenplan", irfinder = run_accession) %>% 
  select(group, replicate, sample, file_salmon, condition, dataset, irfinder)

# Cleveland bacTRAP
cleveland.metadata <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/astrocyte-sod1-mouse-bacTRAP-cleveland-2015/sample-details/samplesheet.csv") %>% filter(grepl("ac", group)) %>%
  mutate(sample = paste(group, "_R", replicate, sep = ""), file_salmon = file.path(dir, paste(group, "_R", replicate, sep = ""), "quant.sf"), condition = gsub("ac_sun_", "", group), group = factor(condition, levels = c("ctrl", "sod1")),
         condition = case_when(group == "sod1" ~ "als", group == "ctrl" ~ "ctrl"), dataset = "cleveland") %>% 
  select(group, replicate, sample, file_salmon, condition, dataset)

# Peng TDP43 alone
peng.metadata <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/astrocyte-tdp43-mouse-peng-2020/sample-details/samplesheet.csv") %>% filter(group != "mou_tdp43_het") %>%
  mutate(sample = paste(group, "_R", replicate, sep = ""), file_salmon = file.path(dir, paste(group, "_R", replicate, sep = ""), "quant.sf"), condition = gsub("mou_tdp43_", "", group), group = factor(condition, levels = c("ctrl", "ko")), condition = case_when(group == "ko" ~ "als", group == "ctrl" ~ "ctrl"), dataset = "peng", irfinder = experiment_accession) %>% 
  select(group, replicate, sample, file_salmon, condition, dataset, irfinder)

# Jiang Membralin alone
jiang.metadata <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/astrocyte-membralin-mouse-jiang-2019/sample-details/samplesheet.csv") %>% 
  mutate(sample = paste(group, "_R", replicate, sep = ""), file_salmon = file.path(dir, paste(group, "_R", replicate, sep = ""), "quant.sf"), condition = gsub("mou_mem_", "", group), group = factor(condition, levels = c("ctrl", "ko")), condition = case_when(group == "ko" ~ "als", group == "ctrl" ~ "ctrl"), dataset = "jiang", irfinder = experiment_accession) %>% 
  select(group, replicate, sample, file_salmon, condition, dataset, irfinder)

# Mouse model meta-analysis: datasets SOD1 TRAP, TDP-43 ko, membralin ko
mouse_als_datasets.metadata <- bind_rows(guttenplan.sod1, peng, jiang) %>% mutate(condition = factor(condition, levels = c("ctrl", "als")), dataset = as.factor(dataset))

# Spinal cord injury model of protective scar forming astrocytes
anderson.metadata <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/astrocyte-sci-mouse-anderson-2017/sample-details/samplesheet.csv") %>% 
  mutate(sample = paste(group, "_R", replicate, sep = ""), file_salmon = file.path(dir, paste(group, "_R", replicate, sep = ""), "quant.sf"), group = factor(group, levels = c("ctrl", "sci")), 
         condition = case_when(group == "sci" ~ "als", group == "ctrl" ~ "ctrl"), dataset = "anderson", irfinder = run_accession) %>% 
  select(group, replicate, sample, file_salmon, condition, dataset, irfinder)

# Zamanian MCAO microarray https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35338 Click **Analyse with GEO2R**
# Define groups and samples: group names (A2 day 1; A0 day 1 sham)  - timepoint at 24h defined by Liddelow https://www.nature.com/articles/nature21029
# Assign samples to group (drag over samples & click group name)
# Options: Select comparisons of interest (A2 day 1 vs A0 sham day 1) >>> Click Reanalyse
# Top differentially expressed genes: Select columns (add Gene symbol, Gene ID) > Set > Download full table
zamanian.A2.res = read_tsv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/reactive-astrocyte-mouse-zamanian-2012/GSE35338.top.table.A2.day1_vs_A0.day1.tsv") %>%
  left_join(mouse2human, by = c("Gene.symbol"="mouse")) %>% drop_na(human) %>% group_by(human) %>% summarise(log2FoldChange = mean(logFC, na.rm = TRUE), padj = mean(adj.P.Val, na.rm = TRUE), pvalue = mean(P.Value, na.rm = TRUE)) %>% rename(gene_name = human) %>% ungroup

# Ferraiuolo sporadic microarray **Analyse with GEO2R**
ferraiuolo.res = read_tsv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/sporadic-human-ferraiuolo-2016/GSE87385.top.table.astrocyte.als_vs_ctrl.tsv") %>%
  mutate(gene_name = str_split_fixed(.$Gene.symbol, "///", 2)[,1]) %>% drop_na(gene_name) %>% select(-Gene.symbol) %>% group_by(gene_name) %>% 
  summarise(log2FoldChange = mean(logFC, na.rm = TRUE), padj = mean(adj.P.Val, na.rm = TRUE), pvalue = mean(P.Value, na.rm = TRUE)) %>% arrange(pvalue)# %>% left_join(gene2ens) # 3101 gene_name unrecognised ensembl IDs, this is duplicating gene_names


# DESeq2 ------------------------------------------------------------------

### hiPSC
patani = DESeq.analysis(metadata = patani.metadata, design = ~ group, contrast = "group_als_vs_ctrl", species = "human", transcript.level = FALSE)
tyzack = DESeq.analysis(metadata = tyzack.metadata, design = ~ group, contrast = "group_als_vs_ctrl", species = "human", transcript.level = FALSE)
birger = DESeq.analysis(metadata = birger.metadata, design = ~ group, contrast = "group_als_vs_ctrl", species = "human", transcript.level = FALSE)
neyrinck = DESeq.analysis(metadata = neyrinck.metadata, design = ~ group, contrast = "group_als_vs_ctrl", species = "human", transcript.level = FALSE)
als_datasets = DESeq.analysis(metadata = als_datasets.metadata, design = ~ dataset + group, contrast = "group_als_vs_ctrl", species = "human", transcript.level = FALSE)
a1_datasets = DESeq.analysis(metadata = a1_datasets.metadata, design = ~ cellline + group, contrast = "group_a1_vs_a0", species = "human", transcript.level = FALSE) # cellline contains dataset & replicate info


### ERROR
zhang_bradley_patani_tyzack_birger_neyrinck = DESeq.analysis(metadata = zhang_bradley_patani_tyzack_birger_neyrinck.metadata, design = ~ 1, contrast = "NA", species = "human", transcript.level = FALSE)
ctrls.zhang_bradley_patani = DESeq.analysis(metadata = ctrls.zhang_bradley_patani.metadata, design = ~ 1, contrast = "NA", species = "human", transcript.level = FALSE)

### mouse
guttenplan.A1 = DESeq.analysis(metadata = guttenplan.A1.metadata, design = ~ replicate + group, contrast = "group_A1_vs_A0", species = "mouse", transcript.level = FALSE) # remove replicate (cellline effects)
guttenplan.sod1 = DESeq.analysis(metadata = guttenplan.sod1.metadata, design = ~ group, contrast = "group_sod1_vs_wt", species = "mouse", transcript.level = FALSE)
cleveland = DESeq.analysis(metadata = cleveland.metadata, design = ~ group, contrast = "group_sod1_vs_ctrl", species = "mouse", transcript.level = FALSE)
peng = DESeq.analysis(metadata = peng.metadata, design = ~ group, contrast = "group_ko_vs_ctrl", species = "mouse", transcript.level = FALSE)
jiang = DESeq.analysis(metadata = jiang.metadata, design = ~ group, contrast = "group_ko_vs_ctrl", species = "mouse", transcript.level = FALSE)
mouse_als_datasets = DESeq.analysis(metadata = mouse_als_datasets.metadata, design = ~ dataset + condition, contrast = "condition_als_vs_ctrl", species = "mouse", transcript.level = FALSE)
anderson = DESeq.analysis(metadata = anderson.metadata, design = ~ group, contrast = "group_sci_vs_ctrl", species = "mouse", transcript.level = FALSE)




# ########## OLD
# anderson.res <- DESeq2::results(anderson.dds, name = "group_sci_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(mus.musculus.gene2ens) %>% mutate(gene_name.mouse = gene_name, gene_name = toupper(gene_name)) %>% arrange(pvalue)
# # anderson.res <- DESeq2::results(anderson.dds, name = "group_sci_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(mus.musculus.gene2ens) %>% left_join(mouse2human, by = c("gene_name"="mouse")) %>% rename(gene_name.mouse = gene_name, gene_name = human) %>% arrange(pvalue)


# patani.files <- patani$file_salmon
# names(patani.files) <- patani$sample
# rownames(patani) <- patani$sample
# patani.txi <- tximport(patani.files, type="salmon", tx2gene=tx2gene)
# patani.dds <- DESeqDataSetFromTximport(patani.txi, colData = patani, design = ~ group)
# patani.dds <- DESeq(patani.dds) # run DESeq2
# patani.vsd <- vst(patani.dds, blind=FALSE)
# patani.res <- DESeq2::results(patani.dds, name = "group_als_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(pvalue)
# 
# tyzack.files <- tyzack$file_salmon
# names(tyzack.files) <- tyzack$sample
# rownames(tyzack) <- tyzack$sample
# tyzack.txi <- tximport(tyzack.files, type="salmon", tx2gene=tx2gene)
# tyzack.dds <- DESeqDataSetFromTximport(tyzack.txi, colData = tyzack, design = ~ group)
# tyzack.dds <- DESeq(tyzack.dds) # run DESeq2
# tyzack.vsd <- vst(tyzack.dds, blind=FALSE)
# tyzack.res <- DESeq2::results(tyzack.dds, name = "group_als_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(pvalue)
# 
# birger.files <- birger$file_salmon
# names(birger.files) <- birger$sample
# rownames(birger) <- birger$sample
# birger.txi <- tximport(birger.files, type="salmon", tx2gene=tx2gene)
# birger.dds <- DESeqDataSetFromTximport(birger.txi, colData = birger, design = ~ group)
# birger.dds <- DESeq(birger.dds) # run DESeq2
# birger.vsd <- vst(birger.dds, blind=FALSE)
# birger.res <- DESeq2::results(birger.dds, name = "group_als_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(pvalue)
# 
# neyrinck.files <- neyrinck$file_salmon
# names(neyrinck.files) <- neyrinck$sample
# rownames(neyrinck) <- neyrinck$sample
# neyrinck.txi <- tximport(neyrinck.files, type="salmon", tx2gene=tx2gene)
# neyrinck.dds <- DESeqDataSetFromTximport(neyrinck.txi, colData = neyrinck, design = ~ group)
# neyrinck.dds <- DESeq(neyrinck.dds) # run DESeq2
# neyrinck.vsd <- vst(neyrinck.dds, blind=FALSE)
# neyrinck.res <- DESeq2::results(neyrinck.dds, name = "group_als_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(pvalue)
# 
# als_datasets.files <- als_datasets$file_salmon
# names(als_datasets.files) <- als_datasets$sample
# rownames(als_datasets) <- als_datasets$sample
# als_datasets.txi <- tximport(als_datasets.files, type="salmon", tx2gene=tx2gene)
# als_datasets.dds <- DESeqDataSetFromTximport(als_datasets.txi, colData = als_datasets, design = ~ dataset + group)
# als_datasets.dds <- DESeq(als_datasets.dds)
# als_datasets.vsd <- vst(als_datasets.dds, blind=FALSE)
# als_datasets.res <- DESeq2::results(als_datasets.dds, name = "group_als_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(pvalue)
# # saveRDS(als_datasets.dds, "/camp/home/ziffo/home/projects/astrocyte-a1-als-meta/expression/deseq2/als_datasets.dds.rds")
# 
# barbar.files <- barbar$file_salmon
# names(barbar.files) <- barbar$sample
# rownames(barbar) <- barbar$sample
# barbar.txi <- tximport(barbar.files, type="salmon", tx2gene=tx2gene)
# barbar.dds <- DESeqDataSetFromTximport(barbar.txi, colData = barbar, design = ~ group)
# barbar.dds <- DESeq(barbar.dds)
# barbar.res <- DESeq2::results(barbar.dds, name = "group_a1_vs_a0") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(pvalue)
# 
# 
# krencik = leng %>% filter(dataset == "Krencik")
# krencik.files <- krencik$file_salmon
# names(krencik.files) <- krencik$sample
# rownames(krencik) <- krencik$sample
# krencik.txi <- tximport(krencik.files, type="salmon", tx2gene=tx2gene)
# krencik.dds <- DESeqDataSetFromTximport(krencik.txi, colData = krencik, design = ~ group)
# krencik.dds <- DESeq(krencik.dds)
# krencik.res <- DESeq2::results(krencik.dds, name = "group_a1_vs_a0") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(pvalue)
# 
# leng.filt = leng %>% filter(dataset == "Leng")
# leng.files <- leng.filt$file_salmon
# names(leng.files) <- leng.filt$sample
# rownames(leng.filt) <- leng.filt$sample
# leng.txi <- tximport(leng.files, type="salmon", tx2gene=tx2gene)
# leng.dds <- DESeqDataSetFromTximport(leng.txi, colData = leng.filt, design = ~ group)
# leng.dds <- DESeq(leng.dds)
# leng.res <- DESeq2::results(leng.dds, name = "group_a1_vs_a0") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(pvalue)
# 
# TCW = leng %>% filter(dataset == "TCW")
# TCW.files <- TCW$file_salmon
# names(TCW.files) <- TCW$sample
# rownames(TCW) <- TCW$sample
# TCW.txi <- tximport(TCW.files, type="salmon", tx2gene=tx2gene)
# TCW.dds <- DESeqDataSetFromTximport(TCW.txi, colData = TCW, design = ~ group)
# TCW.dds <- DESeq(TCW.dds)
# TCW.res <- DESeq2::results(TCW.dds, name = "group_a1_vs_a0") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(pvalue)
# 
# li = leng %>% filter(dataset == "Li")
# li.files <- li$file_salmon
# names(li.files) <- li$sample
# rownames(li) <- li$sample
# li.txi <- tximport(li.files, type="salmon", tx2gene=tx2gene)
# li.dds <- DESeqDataSetFromTximport(li.txi, colData = li, design = ~ group)
# li.dds <- DESeq(li.dds)
# li.res <- DESeq2::results(li.dds, name = "group_a1_vs_a0") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(pvalue)
# 
# # leng.files <- leng$file_salmon
# # names(leng.files) <- leng$sample
# # rownames(leng) <- leng$sample
# # leng.txi <- tximport(leng.files, type="salmon", tx2gene=tx2gene)
# # leng.dds <- DESeqDataSetFromTximport(leng.txi, colData = leng, design = ~ dataset + group)
# # leng.dds <- DESeq(leng.dds)
# # leng.vsd <- vst(leng.dds, blind=FALSE)
# # leng.res <- DESeq2::results(leng.dds, name = "group_a1_vs_a0") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(pvalue)
# 
# a1_datasets.files <- a1_datasets$file_salmon
# names(a1_datasets.files) <- a1_datasets$sample
# rownames(a1_datasets) <- a1_datasets$sample
# a1_datasets.txi <- tximport(a1_datasets.files, type="salmon", tx2gene=tx2gene)
# a1_datasets.dds <- DESeqDataSetFromTximport(a1_datasets.txi, colData = a1_datasets, design = ~ dataset + group)
# a1_datasets.dds <- DESeq(a1_datasets.dds)
# a1_datasets.vsd <- vst(a1_datasets.dds, blind=FALSE)
# a1_datasets.res <- DESeq2::results(a1_datasets.dds, name = "group_a1_vs_a0") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(pvalue)
# # saveRDS(a1_datasets.dds, "/camp/home/ziffo/home/projects/astrocyte-a1-als-meta/expression/deseq2/a1_datasets.dds.rds")
# 
# zhang_bradley_patani_tyzack_birger_neyrinck.files <- zhang_bradley_patani_tyzack_birger_neyrinck.sampleTable$file_salmon
# names(zhang_bradley_patani_tyzack_birger_neyrinck.files) <- zhang_bradley_patani_tyzack_birger_neyrinck.sampleTable$name
# rownames(zhang_bradley_patani_tyzack_birger_neyrinck.sampleTable) <- zhang_bradley_patani_tyzack_birger_neyrinck.sampleTable$name
# zhang_bradley_patani_tyzack_birger_neyrinck.txi <- tximport(zhang_bradley_patani_tyzack_birger_neyrinck.files, type="salmon", tx2gene=tx2gene)
# zhang_bradley_patani_tyzack_birger_neyrinck.dds <- DESeqDataSetFromTximport(zhang_bradley_patani_tyzack_birger_neyrinck.txi, colData = zhang_bradley_patani_tyzack_birger_neyrinck.sampleTable, design = ~ group)
# zhang_bradley_patani_tyzack_birger_neyrinck.vsd <- vst(zhang_bradley_patani_tyzack_birger_neyrinck.dds, blind=FALSE)
# 
# ctrls.zhang_bradley_patani.files <- ctrls.zhang_bradley_patani.sampleTable$file_salmon
# names(ctrls.zhang_bradley_patani.files) <- ctrls.zhang_bradley_patani.sampleTable$name
# rownames(ctrls.zhang_bradley_patani.sampleTable) <- ctrls.zhang_bradley_patani.sampleTable$name
# ctrls.zhang_bradley_patani.txi <- tximport(ctrls.zhang_bradley_patani.files, type="salmon", tx2gene=tx2gene)
# ctrls.zhang_bradley_patani.dds <- DESeqDataSetFromTximport(ctrls.zhang_bradley_patani.txi, colData = ctrls.zhang_bradley_patani.sampleTable, design = ~ 1)
# ctrls.zhang_bradley_patani.vsd <- vst(ctrls.zhang_bradley_patani.dds, blind=FALSE)
# 
# guttenplan.A1.files <- guttenplan.A1$file_salmon
# names(guttenplan.A1.files) <- guttenplan.A1$sample
# rownames(guttenplan.A1) <- guttenplan.A1$sample
# guttenplan.A1.txi <- tximport(guttenplan.A1.files, type="salmon", tx2gene=mus.musculus.tx2gene)
# guttenplan.A1.dds <- DESeqDataSetFromTximport(guttenplan.A1.txi, colData = guttenplan.A1, design = ~ group) ### add CELLLINE to model
# guttenplan.A1.dds <- DESeq(guttenplan.A1.dds)
# guttenplan.A1.vsd <- vst(guttenplan.A1.dds, blind=FALSE)
# guttenplan.A1.res <- DESeq2::results(guttenplan.A1.dds, name = "group_A1_vs_A0") %>% as_tibble(rownames = "gene_id") %>% left_join(mus.musculus.gene2ens) %>% mutate(gene_name.mouse = gene_name, gene_name = toupper(gene_name)) %>% arrange(pvalue)
# 
# guttenplan.sod1.files <- guttenplan.sod1$file_salmon
# names(guttenplan.sod1.files) <- guttenplan.sod1$sample
# rownames(guttenplan.sod1) <- guttenplan.sod1$sample
# guttenplan.sod1.txi <- tximport(guttenplan.sod1.files, type="salmon", tx2gene=mus.musculus.tx2gene)
# guttenplan.sod1.dds <- DESeqDataSetFromTximport(guttenplan.sod1.txi, colData = guttenplan.sod1, design = ~ group)
# guttenplan.sod1.dds <- DESeq(guttenplan.sod1.dds)
# guttenplan.sod1.vsd <- vst(guttenplan.sod1.dds, blind=FALSE)
# guttenplan.sod1.res <- DESeq2::results(guttenplan.sod1.dds, name = "group_sod1_vs_wt") %>% as_tibble(rownames = "gene_id") %>% left_join(mus.musculus.gene2ens) %>% mutate(gene_name.mouse = gene_name, gene_name = toupper(gene_name)) %>% arrange(pvalue)
# 
# 
# cleveland.files <- cleveland$file_salmon
# names(cleveland.files) <- cleveland$sample
# rownames(cleveland) <- cleveland$sample
# cleveland.txi <- tximport(cleveland.files, type="salmon", tx2gene = mus.musculus.tx2gene)
# cleveland.dds <- DESeqDataSetFromTximport(cleveland.txi, colData = cleveland, design = ~ group)
# cleveland.dds <- DESeq(cleveland.dds)
# cleveland.vsd <- vst(cleveland.dds, blind=FALSE)
# resultsNames(cleveland.dds)
# cleveland.res <- DESeq2::results(cleveland.dds, name = "group_sod1_vs_ctrl")  %>% as_tibble(rownames = "gene_id") %>% left_join(mus.musculus.gene2ens) %>% mutate(gene_name.mouse = gene_name, gene_name = toupper(gene_name)) %>% arrange(pvalue)
# 
# peng.files <- peng$file_salmon
# names(peng.files) <- peng$sample
# rownames(peng) <- peng$sample
# peng.txi <- tximport(peng.files, type="salmon", tx2gene = mus.musculus.tx2gene)
# peng.dds <- DESeqDataSetFromTximport(peng.txi, colData = peng, design = ~ group)
# peng.dds <- DESeq(peng.dds) # run DESeq2
# peng.vsd <- vst(peng.dds, blind=FALSE) # VST transformation)
# resultsNames(peng.dds)
# peng.res <- DESeq2::results(peng.dds, name = "group_ko_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(mus.musculus.gene2ens) %>% mutate(gene_name.mouse = gene_name, gene_name = toupper(gene_name)) %>% arrange(pvalue) # join gene gene_nameS to res tibble & order by padj
# 
# jiang.files <- jiang$file_salmon
# names(jiang.files) <- jiang$sample
# rownames(jiang) <- jiang$sample
# jiang.txi <- tximport(jiang.files, type="salmon", tx2gene = mus.musculus.tx2gene)
# jiang.dds <- DESeqDataSetFromTximport(jiang.txi, colData = jiang, design = ~ group)
# jiang.dds <- DESeq(jiang.dds) # run DESeq2
# jiang.vsd <- vst(jiang.dds, blind=FALSE) # VST transformation)
# resultsNames(jiang.dds)
# jiang.res <- DESeq2::results(jiang.dds, name = "group_ko_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(mus.musculus.gene2ens) %>% mutate(gene_name.mouse = gene_name, gene_name = toupper(gene_name)) %>% arrange(pvalue) 
# 
# mouse_als_datasets.files <- mouse_als_datasets$file_salmon
# names(mouse_als_datasets.files) <- mouse_als_datasets$sample
# rownames(mouse_als_datasets) <- mouse_als_datasets$sample
# mouse_als_datasets.txi <- tximport(mouse_als_datasets.files, type="salmon", tx2gene=mus.musculus.tx2gene)
# mouse_als_datasets.dds <- DESeqDataSetFromTximport(mouse_als_datasets.txi, colData = mouse_als_datasets, design = ~ dataset + condition)
# mouse_als_datasets.dds <- DESeq(mouse_als_datasets.dds)
# mouse_als_datasets.vsd <- vst(mouse_als_datasets.dds, blind=FALSE)
# mouse_als_datasets.res <- DESeq2::results(mouse_als_datasets.dds, name = "condition_als_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(mus.musculus.gene2ens) %>% mutate(gene_name.mouse = gene_name, gene_name = toupper(gene_name)) %>% arrange(pvalue)
# # saveRDS(mouse_als_datasets.dds, "/camp/home/ziffo/home/projects/astrocyte-a1-als-meta/expression/deseq2/mouse_als_datasets.dds.rds")
# 
# anderson.files <- anderson$file_salmon
# names(anderson.files) <- anderson$sample
# rownames(anderson) <- anderson$sample
# anderson.txi <- tximport(anderson.files, type="salmon", tx2gene = mus.musculus.tx2gene)
# anderson.dds <- DESeqDataSetFromTximport(anderson.txi, colData = anderson, design = ~ group)
# anderson.dds <- DESeq(anderson.dds) # run DESeq2
# anderson.vsd <- vst(anderson.dds, blind=FALSE) # VST transformation)
# resultsNames(anderson.dds)
# 

# IRFinder 
ir.als_datasets = IRFinder.analysis(metadata = als_datasets, sample.names = "sample", condition = "group", ctrl = "ctrl", mut = "als", batch = "dataset", file.var = "irfinder", ge.res = als_datasets.res, irfinder.dir = "/camp/home/ziffo/home/projects/astrocyte-a1-als-meta/splicing/irfinder")
ir.a1_datasets = IRFinder.analysis(metadata = a1_datasets, sample.names = "sample", condition = "group", ctrl = "a0", mut = "a1", batch = "dataset", file.var = "irfinder", ge.res = a1_datasets.res, irfinder.dir = "/camp/home/ziffo/home/projects/astrocyte-a1-als-meta/splicing/irfinder")
ir.mouse_als_datasets = IRFinder.analysis(metadata = mouse_als_datasets, sample.names = "sample", condition = "condition", ctrl = "ctrl", mut = "als", batch = "dataset", file.var = "irfinder", ge.res = mouse_als_datasets.res, animal = "mouse", irfinder.dir = "/camp/home/ziffo/home/projects/astrocyte-a1-als-meta/splicing/irfinder")
ir.guttenplan.A1 = IRFinder.analysis(metadata = guttenplan.A1, sample.names = "sample", condition = "group", ctrl = "A0", mut = "A1", batch = "NA", file.var = "irfinder", ge.res = guttenplan.A1.res, animal = "mouse", irfinder.dir = "/camp/home/ziffo/home/projects/astrocyte-a1-als-meta/splicing/irfinder")

print("saving")
save.image("/camp/home/ziffo/home/projects/astrocyte-a1-als-meta/expression/deseq2/astrocyte_a1_als_meta.RData")
print("complete")

# ###
# dir <- "/camp/home/ziffo/home/projects/astrocyte-a1-als-meta/alignment/salmon"
# 
# # Patani VCP astrocytes alone
# patani <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/astrocyte-fractionation-bulk-rnaseq/sample-details/sample_table.csv") %>% filter(fraction == "whole", cellline != "ctrl2") %>% distinct(fraction_cellline, .keep_all = TRUE) %>% #
#   mutate(replicate = c(1,2,1,2), sample = cellline, file_salmon = file.path(dir, paste0("ac_who_", vcp, "_R", replicate), "quant.sf"), dataset = "patani",
#          group = gsub("vcp", "als", vcp), group = factor(group, levels = c("ctrl", "als")), region = "ipsc", celltype_source = "This_study", rpt = 1:4, name = paste("astrocyte_ipsc", rpt, sep = "_"), irfinder = paste0("ac_who_",sample)) %>%
#   select(group, replicate, sample, dataset, file_salmon, celltype_source, region, name, irfinder)
# patani.files <- patani$file_salmon
# names(patani.files) <- patani$sample
# rownames(patani) <- patani$sample
# patani.txi <- tximport(patani.files, type="salmon", tx2gene=tx2gene)
# patani.dds <- DESeqDataSetFromTximport(patani.txi, colData = patani, design = ~ group)
# patani.dds <- DESeq(patani.dds) # run DESeq2
# patani.vsd <- vst(patani.dds, blind=FALSE)
# patani.res <- DESeq2::results(patani.dds, name = "group_als_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(pvalue)
# 
# # Tyzack SOD1 astrocytes alone
# tyzack <- read_excel("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/astrocyte-sod1-ipsc-tyzack-2017/sample-details/metadata.xlsx") %>%
#   mutate(replicate = c(1,2,3,4,1,2,3), sample = paste("ac_tyz_", condition, "_R", replicate, sep = ""), file_salmon = file.path(dir, paste0("ac_tyz_", condition, "_R", replicate), "quant.sf"), dataset = "tyzack",
#          group = gsub("sod1", "als", condition), group = factor(group, levels = c("ctrl", "als")), rpt = 1:7, name = paste0("ac_tyz_astrocyte_ipsc_R", rpt), region = "ipsc", celltype_source = "Tyzack_ipsc", irfinder = Experiment) %>%
#   select(group, replicate, sample, dataset, file_salmon, celltype_source, region, name, irfinder)
# tyzack.files <- tyzack$file_salmon
# names(tyzack.files) <- tyzack$sample
# rownames(tyzack) <- tyzack$sample
# tyzack.txi <- tximport(tyzack.files, type="salmon", tx2gene=tx2gene)
# tyzack.dds <- DESeqDataSetFromTximport(tyzack.txi, colData = tyzack, design = ~ group)
# tyzack.dds <- DESeq(tyzack.dds) # run DESeq2
# tyzack.vsd <- vst(tyzack.dds, blind=FALSE)
# tyzack.res <- DESeq2::results(tyzack.dds, name = "group_als_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(pvalue)
# 
# # Birger C9orf72 alone
# birger <- read_excel("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/astrocyte-c9orf72-ipsc-birger-2019/sample-details/metadata.xlsx") %>%
#   mutate(replicate = c(1,2,1,2), condition = gsub("C9orf72", "c9orf72", condition), condition = gsub("control", "ctrl", condition), file_salmon = file.path(dir, paste0("ac_bir_", condition, "_R", replicate), "quant.sf"),
#          group =  gsub("c9orf72", "als", condition), group = factor(group, levels = c("ctrl", "als")), dataset = "birger", rpt = c(1,2,3,4), region = "ipsc", celltype_source = "Birger_ipsc", name = paste0("ac_bir_astrocyte_ipsc_R", rpt),
#          irfinder = paste0("ac_bir_",sample), irfinder = gsub("C9orf72_L", "c9orf", irfinder), irfinder = gsub("control_L","ctrl",irfinder)) %>%
#   select(group, replicate, sample, dataset, file_salmon, celltype_source, region, name,irfinder)
# birger.files <- birger$file_salmon
# names(birger.files) <- birger$sample
# rownames(birger) <- birger$sample
# birger.txi <- tximport(birger.files, type="salmon", tx2gene=tx2gene)
# birger.dds <- DESeqDataSetFromTximport(birger.txi, colData = birger, design = ~ group)
# birger.dds <- DESeq(birger.dds) # run DESeq2
# birger.vsd <- vst(birger.dds, blind=FALSE)
# birger.res <- DESeq2::results(birger.dds, name = "group_als_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(pvalue)
# 
# # Neyrinck FUS astrocytes alone - both day 23 & 48 since both valid astrocytes in original paper
# neyrinck <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/astrocyte-fus-ipsc-neyrinck-2021/sample-details/samplesheet.csv") %>% filter(celltype == "ipsc", overexpressed == "sox9oe", day == 48) %>% # remove primary astrocytes & no TF oe
#   mutate(file_salmon = file.path(dir, paste0(group, "_R", replicate), "quant.sf"), sample = paste("ac_ney_", group, "_R", replicate, sep = ""),
#          condition = group, group = gsub("wt", "ctrl", mutation), group =  gsub("fus", "als", group), group = factor(group, levels = c("ctrl", "als")),
#          dataset = "neyrinck", rpt = c(1:12),
#          region = case_when(grepl("ipsc", condition) ~ "ipsc", TRUE ~ "cortex"), celltype_source = paste0("Neyrinck_", region), name = paste0("ac_ney_astrocyte_", region, "_R", rpt), irfinder = run_accession) %>%
#   select(group, replicate, sample, dataset, file_salmon, celltype_source, region, name,irfinder)
# neyrinck.files <- neyrinck$file_salmon
# names(neyrinck.files) <- neyrinck$sample
# rownames(neyrinck) <- neyrinck$sample
# neyrinck.txi <- tximport(neyrinck.files, type="salmon", tx2gene=tx2gene)
# neyrinck.dds <- DESeqDataSetFromTximport(neyrinck.txi, colData = neyrinck, design = ~ group)
# neyrinck.dds <- DESeq(neyrinck.dds) # run DESeq2
# neyrinck.vsd <- vst(neyrinck.dds, blind=FALSE)
# neyrinck.res <- DESeq2::results(neyrinck.dds, name = "group_als_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(pvalue)
# 
# # # include primary astrocytes from Neyrinck for clustering
# # neyrinck.all <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/astrocyte-fus-ipsc-neyrinck-2021/sample-details/samplesheet.csv") %>%
# #   mutate(file_salmon = file.path(dir, paste0(group, "_R", replicate), "quant.sf"), sample = paste("ac_ney_", group, "_R", replicate, sep = ""),
# #          condition = group, group = gsub("wt", "ctrl", mutation), group =  gsub("fus", "als", group), group = factor(group, levels = c("ctrl", "als")),
# #          dataset = "neyrinck", rpt = c(1:6,1:30), region = case_when(grepl("ipsc", condition) ~ "ipsc", TRUE ~ "cortex"), source = case_when(grepl("fetal", condition) ~ "fetal", grepl("noTFoe", condition) ~ "ipsc",  TRUE ~ "ipsc.SOX9"),
# #          celltype_source = paste0("Neyrinck_", source), name = paste0("ac_ney_astrocyte_", region, "_R", rpt)) %>%
# #  select(group, replicate, sample, dataset, file_salmon, celltype_source, region, name)
# 
# # ALS hiPSC meta-analysis: datasets VCP, SOD1 and C9orf72
# als_datasets <- bind_rows(patani, birger, tyzack, neyrinck) %>% mutate(group = factor(group, levels = c("ctrl", "als")))
# als_datasets.files <- als_datasets$file_salmon
# names(als_datasets.files) <- als_datasets$sample
# rownames(als_datasets) <- als_datasets$sample
# als_datasets.txi <- tximport(als_datasets.files, type="salmon", tx2gene=tx2gene)
# als_datasets.dds <- DESeqDataSetFromTximport(als_datasets.txi, colData = als_datasets, design = ~ dataset + group)
# als_datasets.dds <- DESeq(als_datasets.dds)
# als_datasets.vsd <- vst(als_datasets.dds, blind=FALSE)
# als_datasets.res <- DESeq2::results(als_datasets.dds, name = "group_als_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(pvalue)
# # saveRDS(als_datasets.dds, "/camp/home/ziffo/home/projects/astrocyte-a1-als-meta/expression/deseq2/als_datasets.dds.rds")
# 
# # Barbar cytokine-stimulated A1 astrocytes
# barbar <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/reactive-astrocyte-barbar-2020/sample-details/samplesheet.csv") %>%
#   mutate(file_salmon = file.path(dir, paste(group, "_R", replicate, sep = ""), "quant.sf"), reactivity = c("a0","a1","a0","a1","a0","a1"), group = factor(reactivity, levels = c("a0", "a1")), rpt = c(1,2,3,4,5,6), region = "ipsc", celltype_source = "Barbar_ipsc", name = paste("CD49_astrocyte_ipsc", rpt, sep = "_"), dataset = "barbar", irfinder = sample) %>% select(group, replicate, sample, dataset, file_salmon, celltype_source, region, name, irfinder)
# barbar.files <- barbar$file_salmon
# names(barbar.files) <- barbar$sample
# rownames(barbar) <- barbar$sample
# barbar.txi <- tximport(barbar.files, type="salmon", tx2gene=tx2gene)
# barbar.dds <- DESeqDataSetFromTximport(barbar.txi, colData = barbar, design = ~ group)
# barbar.dds <- DESeq(barbar.dds)
# barbar.res <- DESeq2::results(barbar.dds, name = "group_a1_vs_a0") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(pvalue)
# 
# # Leng a1 astrocytes
# leng <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/reactive-astrocyte-ipsc-leng-2021/sample-details/samplesheet.csv") %>% filter(protocol != "Fernandopulle") %>% # remove protocol Fernandopulle
#   mutate(sample = paste0(group, "_R", replicate), file_salmon = file.path(dir, sample, "quant.sf"), rpt = c(rep(1:6, each = 4)), region = "ipsc", celltype_source = paste0(protocol,"_",region), name = paste0(protocol,"_ipsc_",rpt), group = factor(treatment, levels = c("a0", "a1")), irfinder = Run) %>% select(group, replicate, sample, dataset = protocol, file_salmon, celltype_source, region, name, irfinder)
# 
# krencik = leng %>% filter(dataset == "Krencik")
# krencik.files <- krencik$file_salmon
# names(krencik.files) <- krencik$sample
# rownames(krencik) <- krencik$sample
# krencik.txi <- tximport(krencik.files, type="salmon", tx2gene=tx2gene)
# krencik.dds <- DESeqDataSetFromTximport(krencik.txi, colData = krencik, design = ~ group)
# krencik.dds <- DESeq(krencik.dds)
# krencik.res <- DESeq2::results(krencik.dds, name = "group_a1_vs_a0") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(pvalue)
# 
# leng.filt = leng %>% filter(dataset == "Leng")
# leng.files <- leng.filt$file_salmon
# names(leng.files) <- leng.filt$sample
# rownames(leng.filt) <- leng.filt$sample
# leng.txi <- tximport(leng.files, type="salmon", tx2gene=tx2gene)
# leng.dds <- DESeqDataSetFromTximport(leng.txi, colData = leng.filt, design = ~ group)
# leng.dds <- DESeq(leng.dds)
# leng.res <- DESeq2::results(leng.dds, name = "group_a1_vs_a0") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(pvalue)
# 
# TCW = leng %>% filter(dataset == "TCW")
# TCW.files <- TCW$file_salmon
# names(TCW.files) <- TCW$sample
# rownames(TCW) <- TCW$sample
# TCW.txi <- tximport(TCW.files, type="salmon", tx2gene=tx2gene)
# TCW.dds <- DESeqDataSetFromTximport(TCW.txi, colData = TCW, design = ~ group)
# TCW.dds <- DESeq(TCW.dds)
# TCW.res <- DESeq2::results(TCW.dds, name = "group_a1_vs_a0") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(pvalue)
# 
# li = leng %>% filter(dataset == "Li")
# li.files <- li$file_salmon
# names(li.files) <- li$sample
# rownames(li) <- li$sample
# li.txi <- tximport(li.files, type="salmon", tx2gene=tx2gene)
# li.dds <- DESeqDataSetFromTximport(li.txi, colData = li, design = ~ group)
# li.dds <- DESeq(li.dds)
# li.res <- DESeq2::results(li.dds, name = "group_a1_vs_a0") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(pvalue)
# 
# # leng.files <- leng$file_salmon
# # names(leng.files) <- leng$sample
# # rownames(leng) <- leng$sample
# # leng.txi <- tximport(leng.files, type="salmon", tx2gene=tx2gene)
# # leng.dds <- DESeqDataSetFromTximport(leng.txi, colData = leng, design = ~ dataset + group)
# # leng.dds <- DESeq(leng.dds)
# # leng.vsd <- vst(leng.dds, blind=FALSE)
# # leng.res <- DESeq2::results(leng.dds, name = "group_a1_vs_a0") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(pvalue)
# 
# # a1 meta-analysis
# a1_datasets <- bind_rows(barbar, leng) %>% mutate(group = factor(group, levels = c("a0", "a1")))
# a1_datasets.files <- a1_datasets$file_salmon
# names(a1_datasets.files) <- a1_datasets$sample
# rownames(a1_datasets) <- a1_datasets$sample
# a1_datasets.txi <- tximport(a1_datasets.files, type="salmon", tx2gene=tx2gene)
# a1_datasets.dds <- DESeqDataSetFromTximport(a1_datasets.txi, colData = a1_datasets, design = ~ dataset + group)
# a1_datasets.dds <- DESeq(a1_datasets.dds)
# a1_datasets.vsd <- vst(a1_datasets.dds, blind=FALSE)
# a1_datasets.res <- DESeq2::results(a1_datasets.dds, name = "group_a1_vs_a0") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(pvalue)
# # saveRDS(a1_datasets.dds, "/camp/home/ziffo/home/projects/astrocyte-a1-als-meta/expression/deseq2/a1_datasets.dds.rds")
# 
# # Patani & Tyzack & Birger & Zhang & Barbar & Bradley together - all samples for PCA
# zhang <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/purified-cortex-human-zhang-2016/sample-details/samplesheet.csv") %>%
#   filter(!group %in% c("astrocyte_mouse_cortex", "astrocyte_adult_tumour", "astrocyte_adult_tumour", "astrocyte_adult_hippocampus", "endothelial_adult_temporal_lobe", "whole_adult_temporal_lobe")) %>%
#   mutate(file_salmon = file.path(dir, paste(group, "_R", replicate, sep = ""), "quant.sf"), celltype_source = paste(celltype,source, sep = "_"), region = gsub("temporal_lobe", "cortex", region), name = sample, dataset = "zhang") %>%
#   select(group, replicate, sample, dataset, file_salmon, region, celltype_source, name) # celltype,
# bradley <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/regional-astrocytes-ipsc-bradley-2019/sample-details/samplesheet.csv") %>%
#   mutate(file_salmon = file.path(dir, paste(group, "_R", replicate, sep = ""), "quant.sf"), sample = paste(group, "_R", replicate, sep = ""), region = "ipsc", celltype_source = rep(c("Bradley_FB_ipsc", "Bradley_FB_ipsc", "Bradley_SC_ipsc", "Bradley_SC_ipsc"), each = 6), celltype= "astrocyte", group = "astrocyte_ipsc", name = sample, dataset = "bradley") %>%  select(group, replicate, sample, dataset, file_salmon, region, celltype_source, name) %>% distinct(sample, .keep_all = TRUE)
# 
# zhang_bradley_patani_tyzack_birger_neyrinck.sampleTable <- bind_rows(zhang, patani, tyzack, neyrinck, birger, bradley) %>% mutate(celltype_source = factor(celltype_source, levels = c("This_study", "Tyzack_ipsc", "Birger_ipsc", "Neyrinck_ipsc", "Bradley_SC_ipsc","Bradley_FB_ipsc", "astrocyte_fetal", "astrocyte_adult", "oligodendrocyte_adult", "neuron_adult", "myeloid_adult", "endothelial_adult", "whole_adult"))) # "Neyrinck_ipsc.SOX9", "Neyrinck_fetal",
# zhang_bradley_patani_tyzack_birger_neyrinck.files <- zhang_bradley_patani_tyzack_birger_neyrinck.sampleTable$file_salmon
# names(zhang_bradley_patani_tyzack_birger_neyrinck.files) <- zhang_bradley_patani_tyzack_birger_neyrinck.sampleTable$name
# rownames(zhang_bradley_patani_tyzack_birger_neyrinck.sampleTable) <- zhang_bradley_patani_tyzack_birger_neyrinck.sampleTable$name
# zhang_bradley_patani_tyzack_birger_neyrinck.txi <- tximport(zhang_bradley_patani_tyzack_birger_neyrinck.files, type="salmon", tx2gene=tx2gene)
# zhang_bradley_patani_tyzack_birger_neyrinck.dds <- DESeqDataSetFromTximport(zhang_bradley_patani_tyzack_birger_neyrinck.txi, colData = zhang_bradley_patani_tyzack_birger_neyrinck.sampleTable, design = ~ group)
# zhang_bradley_patani_tyzack_birger_neyrinck.vsd <- vst(zhang_bradley_patani_tyzack_birger_neyrinck.dds, blind=FALSE)
# 
# # CTRLs from this study & Zhang & Bradley together for characterisation
# ctrls.zhang_bradley_patani.sampleTable <- bind_rows(zhang, filter(patani, group == "ctrl"), filter(bradley, grepl("ventral_spinalcord", name))) %>%
#   mutate(celltype_source = factor(celltype_source, levels = c("This_study", "Bradley_VSC_ipsc", "astrocyte_fetal", "astrocyte_adult", "oligodendrocyte_adult", "neuron_adult", "myeloid_adult")))
# ctrls.zhang_bradley_patani.files <- ctrls.zhang_bradley_patani.sampleTable$file_salmon
# names(ctrls.zhang_bradley_patani.files) <- ctrls.zhang_bradley_patani.sampleTable$name
# rownames(ctrls.zhang_bradley_patani.sampleTable) <- ctrls.zhang_bradley_patani.sampleTable$name
# ctrls.zhang_bradley_patani.txi <- tximport(ctrls.zhang_bradley_patani.files, type="salmon", tx2gene=tx2gene)
# ctrls.zhang_bradley_patani.dds <- DESeqDataSetFromTximport(ctrls.zhang_bradley_patani.txi, colData = ctrls.zhang_bradley_patani.sampleTable, design = ~ 1)
# ctrls.zhang_bradley_patani.vsd <- vst(ctrls.zhang_bradley_patani.dds, blind=FALSE)
# 
# # Guttenplan mouse cytokine-stimulated A1 astrocytes
# guttenplan.A1 <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/reactive-astrocyte-mouse-guttenplan-2020/sample-details/samplesheet.csv") %>% filter(celltype == "astrocyte", mutation == "wt") %>%
#   mutate(sample = paste(group, "_R", replicate, sep = ""), file_salmon = file.path(dir, sample, "quant.sf"), reactivity = c("A1","A1","A1","A0","A0","A0"), group = factor(reactivity, levels = c("A0", "A1")), dataset = "guttenplan", irfinder = run_accession) %>%
#   select(group, replicate, sample, dataset, file_salmon, irfinder)
# guttenplan.A1.files <- guttenplan.A1$file_salmon
# names(guttenplan.A1.files) <- guttenplan.A1$sample
# rownames(guttenplan.A1) <- guttenplan.A1$sample
# guttenplan.A1.txi <- tximport(guttenplan.A1.files, type="salmon", tx2gene=mus.musculus.tx2gene)
# guttenplan.A1.dds <- DESeqDataSetFromTximport(guttenplan.A1.txi, colData = guttenplan.A1, design = ~ group)
# guttenplan.A1.dds <- DESeq(guttenplan.A1.dds)
# guttenplan.A1.vsd <- vst(guttenplan.A1.dds, blind=FALSE)
# guttenplan.A1.res <- DESeq2::results(guttenplan.A1.dds, name = "group_A1_vs_A0") %>% as_tibble(rownames = "gene_id") %>% left_join(mus.musculus.gene2ens) %>% mutate(gene_name.mouse = gene_name, gene_name = toupper(gene_name)) %>% arrange(pvalue)
# 
# # Guttenplan mouse SOD1 mutant astrocytes
# guttenplan.sod1 <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/reactive-astrocyte-mouse-guttenplan-2020/sample-details/samplesheet.csv") %>% filter(celltype == "astrocyte", treatment == "untx") %>%
#   mutate(sample = paste(group, "_R", replicate, sep = ""), file_salmon = file.path(dir, sample, "quant.sf"), group = factor(mutation, levels = c("wt", "sod1")), condition = case_when(group == "sod1" ~ "als", group == "wt" ~ "ctrl"), dataset = "guttenplan", irfinder = run_accession) %>%
#   select(group, replicate, sample, file_salmon, condition, dataset, irfinder)
# guttenplan.sod1.files <- guttenplan.sod1$file_salmon
# names(guttenplan.sod1.files) <- guttenplan.sod1$sample
# rownames(guttenplan.sod1) <- guttenplan.sod1$sample
# guttenplan.sod1.txi <- tximport(guttenplan.sod1.files, type="salmon", tx2gene=mus.musculus.tx2gene)
# guttenplan.sod1.dds <- DESeqDataSetFromTximport(guttenplan.sod1.txi, colData = guttenplan.sod1, design = ~ group)
# guttenplan.sod1.dds <- DESeq(guttenplan.sod1.dds)
# guttenplan.sod1.vsd <- vst(guttenplan.sod1.dds, blind=FALSE)
# guttenplan.sod1.res <- DESeq2::results(guttenplan.sod1.dds, name = "group_sod1_vs_wt") %>% as_tibble(rownames = "gene_id") %>% left_join(mus.musculus.gene2ens) %>% mutate(gene_name.mouse = gene_name, gene_name = toupper(gene_name)) %>% arrange(pvalue)
# 
# 
# # Cleveland bacTRAP alone
# cleveland <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/astrocyte-sod1-mouse-bacTRAP-cleveland-2015/sample-details/samplesheet.csv") %>% filter(grepl("ac", group)) %>%
#   mutate(sample = paste(group, "_R", replicate, sep = ""), file_salmon = file.path(dir, paste(group, "_R", replicate, sep = ""), "quant.sf"), condition = gsub("ac_sun_", "", group), group = factor(condition, levels = c("ctrl", "sod1")),
#          condition = case_when(group == "sod1" ~ "als", group == "ctrl" ~ "ctrl"), dataset = "cleveland") %>% select(group, replicate, sample, file_salmon, condition, dataset)
# cleveland.files <- cleveland$file_salmon
# names(cleveland.files) <- cleveland$sample
# rownames(cleveland) <- cleveland$sample
# cleveland.txi <- tximport(cleveland.files, type="salmon", tx2gene = mus.musculus.tx2gene)
# cleveland.dds <- DESeqDataSetFromTximport(cleveland.txi, colData = cleveland, design = ~ group)
# cleveland.dds <- DESeq(cleveland.dds)
# cleveland.vsd <- vst(cleveland.dds, blind=FALSE)
# resultsNames(cleveland.dds)
# cleveland.res <- DESeq2::results(cleveland.dds, name = "group_sod1_vs_ctrl")  %>% as_tibble(rownames = "gene_id") %>% left_join(mus.musculus.gene2ens) %>% mutate(gene_name.mouse = gene_name, gene_name = toupper(gene_name)) %>% arrange(pvalue)
# 
# # Peng TDP43 alone
# peng <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/astrocyte-tdp43-mouse-peng-2020/sample-details/samplesheet.csv") %>% filter(group != "mou_tdp43_het") %>%
#   mutate(sample = paste(group, "_R", replicate, sep = ""), file_salmon = file.path(dir, paste(group, "_R", replicate, sep = ""), "quant.sf"), condition = gsub("mou_tdp43_", "", group), group = factor(condition, levels = c("ctrl", "ko")), condition = case_when(group == "ko" ~ "als", group == "ctrl" ~ "ctrl"), dataset = "peng", irfinder = experiment_accession) %>% select(group, replicate, sample, file_salmon, condition, dataset, irfinder)
# peng.files <- peng$file_salmon
# names(peng.files) <- peng$sample
# rownames(peng) <- peng$sample
# peng.txi <- tximport(peng.files, type="salmon", tx2gene = mus.musculus.tx2gene)
# peng.dds <- DESeqDataSetFromTximport(peng.txi, colData = peng, design = ~ group)
# peng.dds <- DESeq(peng.dds) # run DESeq2
# peng.vsd <- vst(peng.dds, blind=FALSE) # VST transformation)
# resultsNames(peng.dds)
# peng.res <- DESeq2::results(peng.dds, name = "group_ko_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(mus.musculus.gene2ens) %>% mutate(gene_name.mouse = gene_name, gene_name = toupper(gene_name)) %>% arrange(pvalue) # join gene gene_nameS to res tibble & order by padj
# 
# # Jiang Membralin alone
# jiang <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/astrocyte-membralin-mouse-jiang-2019/sample-details/samplesheet.csv") %>%
#   mutate(sample = paste(group, "_R", replicate, sep = ""), file_salmon = file.path(dir, paste(group, "_R", replicate, sep = ""), "quant.sf"), condition = gsub("mou_mem_", "", group), group = factor(condition, levels = c("ctrl", "ko")), condition = case_when(group == "ko" ~ "als", group == "ctrl" ~ "ctrl"), dataset = "jiang", irfinder = experiment_accession) %>% select(group, replicate, sample, file_salmon, condition, dataset, irfinder)
# jiang.files <- jiang$file_salmon
# names(jiang.files) <- jiang$sample
# rownames(jiang) <- jiang$sample
# jiang.txi <- tximport(jiang.files, type="salmon", tx2gene = mus.musculus.tx2gene)
# jiang.dds <- DESeqDataSetFromTximport(jiang.txi, colData = jiang, design = ~ group)
# jiang.dds <- DESeq(jiang.dds) # run DESeq2
# jiang.vsd <- vst(jiang.dds, blind=FALSE) # VST transformation)
# resultsNames(jiang.dds)
# jiang.res <- DESeq2::results(jiang.dds, name = "group_ko_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(mus.musculus.gene2ens) %>% mutate(gene_name.mouse = gene_name, gene_name = toupper(gene_name)) %>% arrange(pvalue)
# 
# # Mouse model meta-analysis: datasets SOD1 TRAP, TDP-43 ko, membralin ko
# mouse_als_datasets <- bind_rows(guttenplan.sod1, peng, jiang) %>% mutate(condition = factor(condition, levels = c("ctrl", "als")))
# mouse_als_datasets.files <- mouse_als_datasets$file_salmon
# names(mouse_als_datasets.files) <- mouse_als_datasets$sample
# rownames(mouse_als_datasets) <- mouse_als_datasets$sample
# mouse_als_datasets.txi <- tximport(mouse_als_datasets.files, type="salmon", tx2gene=mus.musculus.tx2gene)
# mouse_als_datasets.dds <- DESeqDataSetFromTximport(mouse_als_datasets.txi, colData = mouse_als_datasets, design = ~ dataset + condition)
# mouse_als_datasets.dds <- DESeq(mouse_als_datasets.dds)
# mouse_als_datasets.vsd <- vst(mouse_als_datasets.dds, blind=FALSE)
# mouse_als_datasets.res <- DESeq2::results(mouse_als_datasets.dds, name = "condition_als_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(mus.musculus.gene2ens) %>% mutate(gene_name.mouse = gene_name, gene_name = toupper(gene_name)) %>% arrange(pvalue)
# # saveRDS(mouse_als_datasets.dds, "/camp/home/ziffo/home/projects/astrocyte-a1-als-meta/expression/deseq2/mouse_als_datasets.dds.rds")
# 
# # Spinal cord injury model of protective scar forming astrocytes
# anderson <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/astrocyte-sci-mouse-anderson-2017/sample-details/samplesheet.csv") %>%
#   mutate(sample = paste(group, "_R", replicate, sep = ""), file_salmon = file.path(dir, paste(group, "_R", replicate, sep = ""), "quant.sf"), group = factor(group, levels = c("ctrl", "sci")),
#          condition = case_when(group == "sci" ~ "als", group == "ctrl" ~ "ctrl"), dataset = "anderson", irfinder = run_accession) %>% select(group, replicate, sample, file_salmon, condition, dataset, irfinder)
# anderson.files <- anderson$file_salmon
# names(anderson.files) <- anderson$sample
# rownames(anderson) <- anderson$sample
# anderson.txi <- tximport(anderson.files, type="salmon", tx2gene = mus.musculus.tx2gene)
# anderson.dds <- DESeqDataSetFromTximport(anderson.txi, colData = anderson, design = ~ group)
# anderson.dds <- DESeq(anderson.dds) # run DESeq2
# anderson.vsd <- vst(anderson.dds, blind=FALSE) # VST transformation)
# resultsNames(anderson.dds)
# anderson.res <- DESeq2::results(anderson.dds, name = "group_sci_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(mus.musculus.gene2ens) %>% mutate(gene_name.mouse = gene_name, gene_name = toupper(gene_name)) %>% arrange(pvalue)
# # anderson.res <- DESeq2::results(anderson.dds, name = "group_sci_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(mus.musculus.gene2ens) %>% left_join(mouse2human, by = c("gene_name"="mouse")) %>% rename(gene_name.mouse = gene_name, gene_name = human) %>% arrange(pvalue)
# 
# # Zamanian 2012 A1 & A2 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35338Click **Analyse with GEO2R**
# # Define groups:
# # - group names (A2 day 1; A0 day 1 sham)  - timepoint at 24h defined by Liddelow https://www.nature.com/articles/nature21029
# # - assign samples to group (drag over samples & click group name)
# # - logFC only available when only 2 groups defined
# # Options: Select comparisons of interest (A2 day 1 vs A0 sham day 1) >>> Click Reanalyse
# # Top differentially expressed genes:
# # - Select columns (add Gene symbol, Gene ID) > Set
# # - Download full table
# # Read tsv into R
# zamanian.A1.res = read_tsv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/reactive-astrocyte-mouse-zamanian-2012/GSE35338.top.table.A1_vs_A0.tsv") %>%
#   left_join(mouse2human, by = c("Gene.symbol"="mouse")) %>% drop_na(human) %>% group_by(human) %>% summarise(log2FoldChange = mean(logFC, na.rm = TRUE), padj = mean(adj.P.Val, na.rm = TRUE), pvalue = mean(P.Value, na.rm = TRUE)) %>% rename(gene_name = human) %>% ungroup
# zamanian.A2.res = read_tsv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/reactive-astrocyte-mouse-zamanian-2012/GSE35338.top.table.A2.day1_vs_A0.day1.tsv") %>%
#   left_join(mouse2human, by = c("Gene.symbol"="mouse")) %>% drop_na(human) %>% group_by(human) %>% summarise(log2FoldChange = mean(logFC, na.rm = TRUE), padj = mean(adj.P.Val, na.rm = TRUE), pvalue = mean(P.Value, na.rm = TRUE)) %>% rename(gene_name = human) %>% ungroup
# # zamanian.res = read_tsv(
# #   "/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/reactive-astrocyte-mouse-zamanian-2012/GSE35338.top.table.A2.day1_vs_A0.day1.tsv"
# #   # "/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/reactive-astrocyte-mouse-zamanian-2012/GSE35338.top.table.A2.day3_vs_A0.day3.tsv"
# #   # "/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/reactive-astrocyte-mouse-zamanian-2012/GSE35338.top.table.A2.day7_vs_A0.day7.tsv"
# #   # "/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/reactive-astrocyte-mouse-zamanian-2012/GSE35338.top.table.A2.alldays_vs_A0.alldays.tsv"
# #   # "/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/reactive-astrocyte-mouse-zamanian-2012/GSE35338.top.table.A2.alldays_vs_A0.alldays.saline.tsv"
# #   # "/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/reactive-astrocyte-mouse-zamanian-2012/GSE35338.top.table.A2.day1_vs_A0.saline.tsv"
# #   # "/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/reactive-astrocyte-mouse-zamanian-2012/GSE35338.top.table.A1_vs_A0.tsv"
# #   ) %>% left_join(mouse2human, by = c("Gene.symbol"="mouse")) %>% drop_na(human) %>% group_by(human) %>% summarise(log2FoldChange = mean(logFC, na.rm = TRUE), padj = mean(adj.P.Val, na.rm = TRUE), pvalue = mean(P.Value, na.rm = TRUE)) %>% rename(gene_name = human)
# 
# # Ferraiuolo sporadic microarray **Analyse with GEO2R**
# ferraiuolo.res = read_tsv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/sporadic-human-ferraiuolo-2016/GSE87385.top.table.astrocyte.als_vs_ctrl.tsv") %>%
#   mutate(gene_name = str_split_fixed(.$Gene.symbol, "///", 2)[,1]) %>% drop_na(gene_name) %>% select(-Gene.symbol) %>% group_by(gene_name) %>%
#   summarise(log2FoldChange = mean(logFC, na.rm = TRUE), padj = mean(adj.P.Val, na.rm = TRUE), pvalue = mean(P.Value, na.rm = TRUE)) %>% arrange(pvalue)# %>%
# #%>%  left_join(gene2ens) # 3101 gene_name unrecognised ensembl IDs, this is duplicating gene_names
# 
# # IRFinder
# # paths.als_datasets = file.path("/camp/home/ziffo/home/projects/astrocyte-a1-als-meta/splicing/irfinder", als_datasets$irfinder, "IRFinder-IR-dir.txt") # File names vector
# # experiment.als_datasets = als_datasets %>% select(SampleNames = sample, Condition = group, dataset) %>% remove_rownames() %>% mutate(Condition = factor(Condition, levels = c("ctrl", "als")))
# # metaList.als_datasets = DESeqDataSetFromIRFinder(filePaths=paths.als_datasets, designMatrix=experiment.als_datasets, designFormula = ~ dataset + Condition * IRFinder )
# # irfinder.dds.als_datasets = DESeq(metaList.als_datasets$DESeq2Object)
# # resultsNames(irfinder.dds.als_datasets)
# # # [1] "Intercept"                               "condition_vcp_vs_ctrl"
# # # [3] "fraction_nuclear_vs_cytoplasmic"         "IRFinder_IR_vs_Splice"
# # # [5] "conditionvcp.fractionnuclear"            "conditionvcp.IRFinderIR"
# # # [7] "fractionnuclear.IRFinderIR"              "conditionvcp.fractionnuclear.IRFinderIR"
# # irfinder.als_datasets.res = DESeq2::results(irfinder.dds.als_datasets, name="Conditionals.IRFinderIR")
# # irfinder.als_datasets = cleanIRFinder.multidesign(IRFinder = irfinder.als_datasets.res)
# # ir_ge.als_datasets <- summariseIRjoinGE(cleanIRFinder = irfinder.als_datasets, GE = als_datasets.res)
# 
# ir.als_datasets = IRFinder.analysis(metadata = als_datasets, sample.names = "sample", condition = "group", ctrl = "ctrl", mut = "als", batch = "dataset", file.var = "irfinder", ge.res = als_datasets.res, irfinder.dir = "/camp/home/ziffo/home/projects/astrocyte-a1-als-meta/splicing/irfinder")
# ir.a1_datasets = IRFinder.analysis(metadata = a1_datasets, sample.names = "sample", condition = "group", ctrl = "a0", mut = "a1", batch = "dataset", file.var = "irfinder", ge.res = a1_datasets.res, irfinder.dir = "/camp/home/ziffo/home/projects/astrocyte-a1-als-meta/splicing/irfinder")
# ir.mouse_als_datasets = IRFinder.analysis(metadata = mouse_als_datasets, sample.names = "sample", condition = "condition", ctrl = "ctrl", mut = "als", batch = "dataset", file.var = "irfinder", ge.res = mouse_als_datasets.res, animal = "mouse", irfinder.dir = "/camp/home/ziffo/home/projects/astrocyte-a1-als-meta/splicing/irfinder")
# ir.guttenplan.A1 = IRFinder.analysis(metadata = guttenplan.A1, sample.names = "sample", condition = "group", ctrl = "A0", mut = "A1", batch = "NA", file.var = "irfinder", ge.res = guttenplan.A1.res, animal = "mouse", irfinder.dir = "/camp/home/ziffo/home/projects/astrocyte-a1-als-meta/splicing/irfinder")
# 
# print("saving")
# save.image("/camp/home/ziffo/home/projects/astrocyte-a1-als-meta/expression/deseq2/astrocyte_a1_als_meta.RData")
# print("complete")