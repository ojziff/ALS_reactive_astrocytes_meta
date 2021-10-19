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
         celltype_source = "Barbar_ipsc", name = paste("CD49_astrocyte_ipsc", rpt, sep = "_"), dataset = "barbar", irfinder = sample, cellline = paste0(dataset,"_",replicate), group = factor(group, levels = c("a0", "a1")), cellline = as.factor(cellline)) %>% 
  select(group, replicate, sample, dataset, file_salmon, celltype_source, region, name, cellline, irfinder)

# Leng cytokine-stimulated A1 astrocytes - 5 protocols - remove Fernandopulle as odd clustering
leng.protocols.metadata <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/reactive-astrocyte-ipsc-leng-2021/sample-details/samplesheet.csv") %>% 
  filter(protocol != "Fernandopulle") %>% # remove protocol Fernandopulle
  mutate(sample = paste0(group, "_R", replicate), file_salmon = file.path(dir, sample, "quant.sf"), rpt = c(rep(1:6, each = 4)), region = "ipsc", celltype_source = paste0(protocol,"_",region), name = paste0(protocol,"_ipsc_",rpt), group = factor(treatment, levels = c("a0", "a1")), irfinder = run_accession) %>% 
  select(group, replicate, sample, dataset = protocol, file_salmon, celltype_source, region, name, cellline, irfinder)

leng.metadata  = leng.protocols.metadata %>% filter(dataset == "Leng") %>% mutate(group = factor(group, levels = c("a0", "a1")), cellline = as.factor(cellline))
li.metadata = leng.protocols.metadata %>% filter(dataset == "Li") %>% mutate(group = factor(group, levels = c("a0", "a1")), cellline = as.factor(cellline))
krencik.metadata = leng.protocols.metadata %>% filter(dataset == "Krencik") %>% mutate(group = factor(group, levels = c("a0", "a1")), cellline = as.factor(cellline))
TCW.metadata = leng.protocols.metadata %>% filter(dataset == "TCW") %>% mutate(group = factor(group, levels = c("a0", "a1")), cellline = as.factor(cellline))

# a1 meta-analysis
a1_datasets.metadata <- bind_rows(barbar.metadata, leng.metadata, li.metadata, krencik.metadata, TCW.metadata) %>% mutate(group = factor(group, levels = c("a0", "a1")), cellline = as.factor(cellline))

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
mouse_als_datasets.metadata <- bind_rows(guttenplan.sod1.metadata, peng.metadata, jiang.metadata) %>% mutate(condition = factor(condition, levels = c("ctrl", "als")), dataset = as.factor(dataset))

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

barbar = DESeq.analysis(metadata = barbar.metadata, design = ~ cellline + group, contrast = "group_a1_vs_a0", species = "human", transcript.level = FALSE)
leng = DESeq.analysis(metadata = leng.metadata, design = ~ cellline + group, contrast = "group_a1_vs_a0", species = "human", transcript.level = FALSE)
li = DESeq.analysis(metadata = li.metadata, design = ~ cellline + group, contrast = "group_a1_vs_a0", species = "human", transcript.level = FALSE)
krencik = DESeq.analysis(metadata = krencik.metadata, design = ~ cellline + group, contrast = "group_a1_vs_a0", species = "human", transcript.level = FALSE)
TCW = DESeq.analysis(metadata = TCW.metadata, design = ~ cellline + group, contrast = "group_a1_vs_a0", species = "human", transcript.level = FALSE)
a1_datasets = DESeq.analysis(metadata = a1_datasets.metadata, design = ~ cellline + group, contrast = "group_a1_vs_a0", species = "human", transcript.level = FALSE) # cellline contains dataset & replicate info

zhang_bradley_patani_tyzack_birger_neyrinck = DESeq.analysis(metadata = zhang_bradley_patani_tyzack_birger_neyrinck.metadata, design = ~ 1, contrast = "NA", species = "human", transcript.level = FALSE)
ctrls.zhang_bradley_patani = DESeq.analysis(metadata = ctrls.zhang_bradley_patani.metadata, design = ~ 1, contrast = "NA", species = "human", transcript.level = FALSE)

print("saving")
save.image("/camp/home/ziffo/home/projects/astrocyte-a1-als-meta/expression/deseq2/astrocyte_a1_als_meta.RData")
print("running mouse")

### mouse
guttenplan.A1 = DESeq.analysis(metadata = guttenplan.A1.metadata, design = ~ replicate + group, contrast = "group_A1_vs_A0", species = "mouse", transcript.level = FALSE) # remove replicate (cellline effects)
guttenplan.sod1 = DESeq.analysis(metadata = guttenplan.sod1.metadata, design = ~ group, contrast = "group_sod1_vs_wt", species = "mouse", transcript.level = FALSE)
cleveland = DESeq.analysis(metadata = cleveland.metadata, design = ~ group, contrast = "group_sod1_vs_ctrl", species = "mouse", transcript.level = FALSE)
peng = DESeq.analysis(metadata = peng.metadata, design = ~ group, contrast = "group_ko_vs_ctrl", species = "mouse", transcript.level = FALSE)
jiang = DESeq.analysis(metadata = jiang.metadata, design = ~ group, contrast = "group_ko_vs_ctrl", species = "mouse", transcript.level = FALSE)
mouse_als_datasets = DESeq.analysis(metadata = mouse_als_datasets.metadata, design = ~ dataset + condition, contrast = "condition_als_vs_ctrl", species = "mouse", transcript.level = FALSE)
anderson = DESeq.analysis(metadata = anderson.metadata, design = ~ group, contrast = "group_sci_vs_ctrl", species = "mouse", transcript.level = FALSE)

# anderson.res <- DESeq2::results(anderson.dds, name = "group_sci_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(mus.musculus.gene2ens) %>% mutate(gene_name.mouse = gene_name, gene_name = toupper(gene_name)) %>% arrange(pvalue)
# # anderson.res <- DESeq2::results(anderson.dds, name = "group_sci_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(mus.musculus.gene2ens) %>% left_join(mouse2human, by = c("gene_name"="mouse")) %>% rename(gene_name.mouse = gene_name, gene_name = human) %>% arrange(pvalue)

print("saving")
save.image("/camp/home/ziffo/home/projects/astrocyte-a1-als-meta/expression/deseq2/astrocyte_a1_als_meta.RData")
print("running IRFinder")

# IRFinder ----------------------------------------------------------------
ir.als_datasets = IRFinder.analysis(metadata = a1_datasets.metadata, sample.names = "sample", condition = "group", ctrl = "ctrl", mut = "als", batch = "dataset", file.var = "irfinder", ge.res = als_datasets$res, irfinder.dir = "/camp/home/ziffo/home/projects/astrocyte-a1-als-meta/splicing/irfinder")
ir.a1_datasets = IRFinder.analysis(metadata = a1_datasets.metadata, sample.names = "sample", condition = "group", ctrl = "a0", mut = "a1", batch = "dataset", file.var = "irfinder", ge.res = a1_datasets$res, irfinder.dir = "/camp/home/ziffo/home/projects/astrocyte-a1-als-meta/splicing/irfinder")
ir.mouse_als_datasets = IRFinder.analysis(metadata = mouse_als_datasets.metadata, sample.names = "sample", condition = "condition", ctrl = "ctrl", mut = "als", batch = "dataset", file.var = "irfinder", ge.res = mouse_als_datasets$res, animal = "mouse", irfinder.dir = "/camp/home/ziffo/home/projects/astrocyte-a1-als-meta/splicing/irfinder")
ir.guttenplan.A1 = IRFinder.analysis(metadata = guttenplan.A1.metadata, sample.names = "sample", condition = "group", ctrl = "A0", mut = "A1", batch = "NA", file.var = "irfinder", ge.res = guttenplan.A1$res, animal = "mouse", irfinder.dir = "/camp/home/ziffo/home/projects/astrocyte-a1-als-meta/splicing/irfinder")

print("saving")
save.image("/camp/home/ziffo/home/projects/astrocyte-a1-als-meta/expression/deseq2/astrocyte_a1_als_meta.RData")
print("complete")
