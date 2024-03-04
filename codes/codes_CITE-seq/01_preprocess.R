
library(Seurat)
library(plyr)
library(scater)
library(future)
library(ggplot2)
library(cowplot)
library(reshape2)
library(data.table)
library(tidyverse)
library(future.apply)
library(RColorBrewer)
library(SingleCellExperiment)

load("gtf_v22.Rdata")
source("citeseq_function.R")
dsb.norm = TRUE


# readdata ========================================================

filename <- list.files("data")
filename <- setdiff(filename, "raw_mat")

allexp = list()
allpro = list()

for (i in filename) {
  
  name <- paste0("data/", i)
  data <- Read10X(name)
  
  dataexp <- data[["Gene Expression"]]
  dataexp <- CreateSeuratObject(dataexp)
  dataexp$sample <- i
  
  datapro <- data[["Antibody Capture"]]
  datapro <- CreateSeuratObject(datapro)
  datapro$sample <- i
  
  dataexp <- RenameCells(dataexp, add.cell.id = i)
  datapro <- RenameCells(datapro, add.cell.id = i)
  
  allexp <- c(allexp, list(dataexp))
  allpro <- c(allpro, list(datapro))
}

names(allexp) <- filename
names(allpro) <- filename



# DSB normalization =========================================================

filename <- list.files("data")
filename <- setdiff(filename, "raw_mat")

allexp = list()
allpro = list()

for (i in filename) {
  
  name <- paste0("data/", i)
  data <- Read10X(name)
  
  dataexp <- data[["Gene Expression"]]
  dataexp <- CreateSeuratObject(dataexp)
  dataexp$sample <- i
  
  countpro <- data[["Antibody Capture"]]
  name1 <- paste0("data/raw_mat/", i)
  name2 <- paste0("data/", i)
  
  raw = Seurat::Read10X(name1)
  cells = Seurat::Read10X(name2)
  
  datapro <- dsb_norm(raw, cells, isotype = isotype)
  countpro <- countpro[rownames(datapro), colnames(datapro)]
  
  datapro <- CreateSeuratObject(datapro)
  datapro@assays[["RNA"]]@counts <- countpro
  datapro$sample <- i
  
  dataexp <- RenameCells(dataexp, add.cell.id = i)
  datapro <- RenameCells(datapro, add.cell.id = i)
  
  allexp <- c(allexp, list(dataexp))
  allpro <- c(allpro, list(datapro))
}

names(allexp) <- filename
names(allpro) <- filename

# QC =====================================================================

for (i in names(allexp)) {
  subdata <- allexp[[i]]
  
  name = paste0("figure/01_basic/QC_plot_", i, ".pdf")
  pdf(file = name, width = 10, height = 7)
  qc_plot(data = subdata)
  dev.off()
}

allexp <- lapply(allexp, function(i){
  qc_filter(data = i,
            nFeature_dn = 1000,
            nFeature_up = 7000,
            nCount_dn = 4000,
            nCount_up = 40000,
            mito_up = 15,
            ribo_dn = 5,
            ribo_up = 50,
            gene_filter = 3,
            species = "human")
})

plan(multisession, workers = length(allexp))
options(future.globals.maxSize = 20000*1024^2)

allexp <- future_lapply(allexp, function(i){
  doublet_rm(i, prop_doublet = 0.05)
}, future.seed = TRUE)
plan(sequential)

allexp <- lapply(allexp, function(i){
  common <- intersect(rownames(i), gene_info$SYMBOL)
  i[common,]
})


# ADT QC --------------------

allpro <- lapply(allpro, function(i){
  
  counts <- GetAssayData(i, slot = "counts", assay = "RNA")
  data <- GetAssayData(i, slot = "data", assay = "RNA")
  
  rownames(counts) <- do.call(rbind, strsplit(rownames(counts), '_'))[,2]
  rownames(data) <- do.call(rbind, strsplit(rownames(data), '-'))[,2]
  
  object <- CreateSeuratObject(data)
  object@assays[["RNA"]]@counts <- counts
  object
})

# allpro <- lapply(allpro, function(i){
#   pro <- i
#   select <- c("IgG1", "IgG2a", "IgG2b", "IgG2a.1")
#   data <- exp_dicho(pro, select)
#   pro[,rownames(data)[rowSums(data) < 2]]
# })

for (j in names(allpro)) {
  
  all_plots <- protein_distribution(data = allpro[[j]], 
                                    dsb.norm = dsb.norm, thres = 7)
  
  all_plots <- plot_grid(plotlist = all_plots, ncol = 5)
  name = paste0("figure/01_basic/dense_protein_", j, ".png")
  ggsave(name, all_plots, dpi = 200,
         width = 30, height = 20, limitsize = FALSE)
}

allprop <- do.call(cbind, lapply(names(allpro), function(i){
  subpro <- as.matrix(GetAssayData(allpro[[i]], slot = "data", assay = "RNA"))
  subpro <- ifelse(subpro > 6.5, 1, 0)
  prop <- data.frame(rowSums(subpro) / ncol(subpro))
  names(prop) <- i
  prop
}))

name = "figure/01_basic/heatmap_protein_prop.pdf"
pdf(file = name, width = 7, height = 8)
heatmap_text(input = allprop, col_rot = 45, 
             color = "white2blue", cutoff = 0,
             order_name = F, cluster_row = T, cluster_col = F)
dev.off()

select_al <- c("CD13", "CD24", "CD44", "CD47", "CD49f", "CD54",
               "CD73", "CD90", "CD117", "CD133", "CD326", "CD338",
               "CD183", "CD184", 
               "CD274", "CD155",
               "IgG1", "IgG2a", "IgG2b")

allpro <- lapply(allpro, function(i){
  i[select_al,]})


# Integration --------------------

name = names(allexp)
dataexp = list()
datapro = list()

for (i in name) {
  
  subexp <- allexp[[i]]; subpro <- allpro[[i]]
  common <- intersect(colnames(subexp), colnames(subpro))
  subexp <- subexp[,common]; subpro <- subpro[,common]
  
  dataexp <- c(dataexp, list(subexp))
  datapro <- c(datapro, list(subpro))
}

names(dataexp) <- name; allexp <- dataexp
names(datapro) <- name; allpro <- datapro


# savedata --------------------

save(allexp, allpro, file = "savedata/data1.Rdata")



# normalization =====================================================================

allexp <- lapply(allexp, function(i){
  subexp <- autocluster_exp(i, nfeatures = 2000, ndim = 20,
                            neigh = 20, dist = 1, res = 0.5)
  
  subexp$cluster_exp <- subexp$seurat_clusters
  subexp
})

all_plots <- lapply(names(allexp), function(i){
  dimplot_new(allexp[[i]], pt.size = 0.75, label = T,
              group.by = c("seurat_clusters")) + ggtitle(i)
})

name = "figure/01_basic/umap_cluster_exp_sep.png"
export_dimplot(all_plots = all_plots, ncol = 2,
               dpi = 300, output = name)

select <- c("CD13", "CD24", "CD44", "CD47", "CD49f", "CD54",
            "CD73", "CD90", "CD117", "CD133", "CD326", "CD338", "CD184")

allpro <- lapply(allpro, function(i){
  subpro <- autocluster_pro(i, ndim = 5, features = select,
                            neigh = 20, dist = 1, res = 0.5,
                            dsb.norm = dsb.norm)
  
  subpro$cluster_pro <- subpro$seurat_clusters
  subpro
})

all_plots <- lapply(names(allpro), function(i){
  dimplot_new(allpro[[i]], pt.size = 0.75, label = T,
              group.by = c("seurat_clusters")) + ggtitle(i)
})

name = "figure/01_basic/umap_cluster_pro_sep.png"
export_dimplot(all_plots = all_plots, ncol = 2,
               dpi = 300, output = name)


# merge data =============================================================

megexp <- allexp[[1]]

for (i in 2:length(allexp)) {
  data <- allexp[[i]]
  megexp <- merge(megexp, data)
}

megexp <- autocluster_exp(megexp, nfeatures = 2000, 
                          ndim = 20, neigh = 20, dist = 1, res = 0.5)

plot <- dimplot_new(megexp, pt.size = 0.75, label = T,
                    group.by = c("sample"))

name = "figure/01_basic/umap_cluster_exp_meg.png"
ggsave(name, plot, dpi = 300, width = 8.5, height = 7)

megpro <- allpro[[1]]

for (i in 2:length(allpro)) {
  data <- allpro[[i]]
  megpro <- merge(megpro, data)
}

megpro <- autocluster_pro(megpro, ndim = 5,
                          features = rownames(megpro),
                          neigh = 20, dist = 1, res = 0.5,
                          dsb.norm = dsb.norm)

plot <- dimplot_new(megpro, pt.size = 0.75, label = T,
                    group.by = c("orig.ident"))

name = "figure/01_basic/umap_cluster_pro_meg.png"
ggsave(name, plot, dpi = 300, width = 8.5, height = 7)

input <- data.frame(table(megexpmc$sample))

plot <- ggplot(input, aes(x = Var1, y = Freq))+
  geom_bar(stat = 'identity', width = 0.6) + 
  theme_classic() + coord_flip()

name = "figure/01_basic/barplot_allexpmc.pdf"
ggsave(name, plot, width = 5, height = 7)


# scanpy --------------------

library(sceasy)
library(reticulate)
use_condaenv('scanpy')

name = "/home/ug0302/CITEseq/savedata/megexp.h5ad"
convertFormat(megexp, outFile = name,
              from = "seurat", to = "anndata", assay = "RNA")


# legoplot --------------------

input <- GetAssayData(megpro, slot = "data", assay = "RNA")
input <- SetAssayData(megpro, slot = "counts",
                      new.data = input, assay = "RNA")

input <- AverageExpression(input, group.by = "orig.ident",
                           slot = "counts", return.seurat = FALSE)[[1]]

# input <- t(input)
input <- t(allprop)
write.table(input, "figure/01_basic/input_legoplot.txt",
            sep = "\t", quote = F, row.names = T, col.names = NA)

save(allexp, allpro, file = "savedata/data1.Rdata")
save(megexp, megpro, file = "savedata/data1_meg.Rdata")



# variation evaluation ===========================================================

input <- GetAssayData(megpro, slot = "data", assay = "RNA")
input <- SetAssayData(megpro, slot = "counts",
                      new.data = input, assay = "RNA")

input <- AverageExpression(input, group.by = "orig.ident",
                           slot = "counts", return.seurat = FALSE)[[1]]
input <- t(input)
# input <- input[!(rownames(input) %in% c("HC031", "PLCPR5")),]

cv_scores <- do.call(rbind, lapply(1:ncol(input), function(i){
  score <- sd(input[,i]) / mean(input[,i])
  data.frame(id = colnames(input)[i], score = score)
}))


# ITH --------------------

library(nlme)

input <- as.matrix(GetAssayData(input, slot = "data", assay = "RNA"))
input <- input[,sample(colnames(input), 10000)]

id <- do.call(rbind, strsplit(colnames(input), '_'))[,1]
input <- data.frame(id = id, scale(t(input)))

var_result = data.frame()
for(i in colnames(input[,2:ncol(input)])){
  
  a <- input[i]
  a <- unlist(a)
  lme <- lme(a~1, random = ~1|id, data = input)
  var <- data.frame(as.character(VarCorr(lme)))[,1] %>% as.numeric()
  
  var_sum <- data.frame(id = i,
                        between_group = var[1],
                        within_group = var[2],
                        std1 = var[3],
                        std2 = var[4],
                        all = var[1]+var[2],
                        WT = var[2]/(var[1]+var[2]))
  
  var_result <- rbind(var_result, var_sum)
}

input <- as.matrix(GetAssayData(megpro, slot = "data", assay = "RNA"))
exp <- rowMeans(input)
exp <- data.frame(id = names(exp), expression = unname(exp))
var <- data.frame(id = var_result$id, variance = var_result$between_group)
input <- merge(exp, var, by = "id")

plot <- scatter_label(input = input,
                      size.p = 3.7, size.l = 5)

name = "figure/01_basic/variance_protein.pdf"
ggsave(name, plot, width = 7.2, height = 7)

