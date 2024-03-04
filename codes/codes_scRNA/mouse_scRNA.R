library(Seurat)
library(plyr)
library(scater)
library(stringr)
library(future)
library(ggplot2)
library(cowplot)
library(reshape2)
library(data.table)
library(tidyverse)
library(future.apply)
library(RColorBrewer)
library(SingleCellExperiment)
source("public_function_mouse.R")
source("citeseq_function.R")


# pre-process =====================================================================

path = "/home/ug0302/CITEseq3/03_scRNA/data"
data <- read10x_all(path = path, merge = T)
data$group <- substring(data$sample, 7, 7)

name = "figure/01_basic/QC_plot.pdf"
pdf(file = name, width = 10, height = 7)
qc_plot(data = data, species = "mouse")
dev.off()

datafilt <- qc_filter(data = data,
                      nFeature_dn = 300,
                      nFeature_up = 7000,
                      nCount_dn = 1000,
                      nCount_up = 40000,
                      mito_up = 15,
                      ribo_dn = 1,
                      ribo_up = 50,
                      gene_filter = 3,
                      species = "mouse")

load("anno/gtf_v22_mouse.Rdata")
common <- intersect(rownames(datafilt), gene_info$SYMBOL)
datafilt <- datafilt[common,]
saveRDS(datafilt, "data/datafilt_raw.rds")



# annotate =====================================================================

source("/home/ug0302/CITEseq3/public_data/autoan_function_mus.R")
scGate_DB <- readRDS("/home/ug0302/CITEseq3/public_data/scGate_DB.rds")

# gene_anno

non_epi <- c("Krt5-", "Krt14-", "Krt6a-", "Dsp-", "Krt17-", "Lgals7-")

# Immune annotation

dataimmu <- auto_immune(datafilt, scGate_DB = scGate_DB,
                        non_epi = non_epi, min_cell = 100, ncore = 30)

# validation ====================

# umap

filename = "figure/01_basic/umap_celltype_sig2.pdf"
dimplot_hca(dataimmu, filename = filename,
            psize = 1, width = 9, height = 7.5)

# marker

Tcell = c("Cd3d", "Cd3e")
CD8T = c("Cd8a", "Cd8b1")
CD4T = c("Cd4")
Treg = c("Foxp3", "Il2ra", "Il7r")
NKcell = c("Gzma", "Klrd1", "Klrk1")
Bcell = c("Cd19", "Cd79a", "Ms4a1", "Igkc")
Macro = c("Apoe", "Arg1", "Lyz2", "C1qa", "C1qb", "C1qc")
DC = c("Cst3", "H2-Aa", "H2-Ab1", "H2-Eb1")
Neu = c("Csf3r", "Cxcr2", "Clec4d", "S100a8", "S100a9")

gene_list <- list(Tcell = Tcell,
                  CD8T = CD8T,
                  CD4T = CD4T,
                  Treg = Treg,
                  NKcell = NKcell,
                  Bcell = Bcell,
                  Macro = Macro,
                  DC = DC,
                  Neu = Neu)

name = "figure/01_basic/marker_dotplot.pdf"
dotplot_marker(dataimmu,
               group.by = "celltype_sig2",
               marker = gene_list,
               species = NULL,
               output = name,
               height = 6)

# Visualization1

plot <- dimplot_new(dataimmu,
                    reduction = "umap",
                    pt.size = 0.2, label = T,
                    group.by = c("seurat_clusters"))

name = "figure/01_basic/umap_clusters.png"
ggsave(name, plot, dpi = 300, width = 6, height = 5)

# Visualization2

plot <- dimplot_new(dataimmu,
                    reduction = "umap",
                    pt.size = 0.5, label = T,
                    group.by = c("sample"))

name = "figure/01_basic/umap_sample.png"
ggsave(name, plot, dpi = 300, width = 6, height = 5)

# Visualization3

plot <- dimplot_new(dataimmu,
                    reduction = "umap",
                    pt.size = 0.5, label = T,
                    group.by = c("group"))

name = "figure/01_basic/umap_group.png"
ggsave(name, plot, dpi = 300, width = 6, height = 5)

# marker visualization

name = "figure/01_basic/markers_umap.png"
featureplot_marker(dataimmu,
                   reduction = "umap",
                   pt.size = 0.5,
                   outlier.rm = F,
                   species = "mouse",
                   output = name)



# prop diff ===================================================================

# umap density

plot <- prop_density(datafilt = dataimmu,
                     group = "group",
                     coord = "umap")

name = "figure/01_basic/prop_density.pdf"
ggsave(name, plot, width = 10, height = 5)

# barplot

plot <- prop_back2back(datafilt = dataimmu,
                       group = "group",
                       cluster = "celltype_sig2",
                       order = TRUE)

name = "figure/01_basic/prop_back2back.pdf"
ggsave(name, plot, width = 7, height = 5)



# CD8 signature ==================================================================

# signature score ----------

dataimmu <- readRDS("data/dataimmu_flt.rds")
subfilt <- dataimmu[,dataimmu$celltype_sig2 == "CD8T"]

path = "/home/ug0302/CITEseq3/03_scRNA/signature/sig_collection1.csv"
scores <- seurat_score(subfilt, source = path,
                       geneset = NULL,
                       min.sz = 3)

subfilt <- AddMetaData(subfilt, scores)

info <- subfilt@meta.data
input <- data.frame(type = info$group,
                    value = info$sig1_Progenitor_Ex)

plot <- ecdf_plot(input, ylim = c(0,1))
name = "figure/03_CD8/ecdf_CD8_sig1_Progenitor_Ex.pdf"
ggsave(name, plot, width = 6, height = 5)



# Neutrophil signature ===============================================================

dataimmu <- readRDS("data/dataimmu.rds")
subfilt <- dataimmu[,dataimmu$celltype_sig2 == "Neutrophil"]

path = "/home/ug0302/CITEseq3/03_scRNA/signature/neutrophil_sig.csv"
scores <- seurat_score(subfilt, source = path,
                       geneset = NULL,
                       min.sz = 3)

subfilt <- AddMetaData(subfilt, scores)

info <- subfilt@meta.data
input <- data.frame(type = info$group,
                    value = info$mNeu_10_Ifit1)

plot <- ecdf_plot(input, ylim = c(0,1))
name = "figure/02_neutrophil/ecdf_neu_mNeu_10_Ifit1.pdf"
ggsave(name, plot, width = 6, height = 5)

# JSD ----------

# L-R data

name = "/home/ug0302/CITEseq2/data/LRanalysis/NicheNet_LR_pairs.csv"
lrinfo <- read.table(name, sep = ",", header = T, check.names = F)

load("/home/ug0302/CITEseq2/data/LRanalysis/Science.Rdata")
membrane = protein_anno[["Science_protein_class"]][["Predicted membrane proteins"]]

load("anno/gtf_v22_mouse.Rdata")
membrane_mus <- gene_info$SYMBOL[gene_info$human_gene_symbol %in% membrane]
receptor_mus <- gene_info$SYMBOL[gene_info$human_gene_symbol %in% unique(lrinfo$to)]

receptor <- intersect(receptor_mus, rownames(dataimmu))
receptor <- intersect(receptor, membrane_mus)

data_r <- GetAssayData(dataimmu, slot = "data", assay = "RNA")[receptor,]
data_r <- as.matrix(data_r)

cluster <- data.frame(id = colnames(dataimmu), 
                      type = dataimmu$celltype_sig2)

jsd_r <- jsd_score(data = data_r,
                   cluster = cluster,
                   select = "Neutrophil")

exp_r <- meanexpression(dataimmu[receptor,], slot = "data",
                        group = "celltype_sig2", method = "mean")

exp_r <- data.frame(id = rownames(exp_r), exp = exp_r[,"Neutrophil"])
input_r <- merge(exp_r, jsd_r, by = "id")

plot <- scatter_highlight(input_r, topn = 5, color = "#0070b2")
name = "figure/02_neutrophil/receptor_jsd.pdf"
ggsave(name, plot, width = 6, height = 5.5)


