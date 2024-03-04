library(Seurat)
library(plyr)
library(future)
library(gtools)
library(ggplot2)
library(cowplot)
library(data.table)
library(tidyverse)
library(future.apply)
library(RColorBrewer)
source("citeseq_function.R")

select_dm <- c("CD13", "CD24", "CD44", "CD47", "CD49f", "CD54",
               "CD73", "CD90", "CD117", "CD133", "CD326", "CD338", "CD184")

select_vi <- c("CD13", "CD24", "CD44", "CD49f", "CD54",
               "CD73", "CD90", "CD117", "CD133", "CD326", "CD338",
               "CD183", "CD184", 
               "CD47", "CD274", "CD155")

select_al <- c("CD13", "CD24", "CD44", "CD49f", "CD54",
               "CD73", "CD90", "CD117", "CD133", "CD326", "CD338",
               "CD183", "CD184", 
               "CD47", "CD274", "CD155",
               "IgG1", "IgG2a", "IgG2b")

dsb.norm = TRUE


# metacell =================================================================

plan(multisession, workers = length(allexp))
options(future.globals.maxSize = 20000*1024^2)

allexpmc <- future_lapply(names(allexp), function(i){
  metacell_multiomics(allexp[[i]], resolution = 50, integ.method = "tvi",
                      cluster.name = paste0(i, "_metacell"))
  
}, future.seed = TRUE)
plan(sequential)
names(allexpmc) <- names(allexp)

# metacell - pro ----------

allpromc <- lapply(names(allexpmc), function(i){
  
  metacell <- allexpmc[[i]]
  rawcell <- allpro[[i]]
  
  subpromc <- metacell_object_pro(metacell = metacell, rawcell = rawcell)
  subpromc <- autocluster_pro(subpromc, ndim = 5, features = select_dm,
                              neigh = 20, dist = 0.5, res = 0.5,
                              dsb.norm = dsb.norm)
})

names(allpromc) <- names(allexpmc)
allexpmc3 <- allexpmc
allpromc3 <- allpromc

# dropout prop ============================================================

dropout1 <- do.call(rbind, lapply(names(allexp), function(i){
  mat <- GetAssayData(allexp[[i]], slot = "data", assay = "RNA")
  prop <- data.frame(type = "raw", sample = i,
                     value = sum(mat == 0)/(ncol(mat) * nrow(mat)))
}))
  
dropout2 <- do.call(rbind, lapply(names(allexpmc), function(i){
  mat <- GetAssayData(allexpmc[[i]], slot = "data", assay = "RNA")
  prop <- data.frame(type = "meta", sample = i,
                     value = sum(mat == 0)/(ncol(mat) * nrow(mat)))
}))

input <- rbind(dropout1, dropout2)

name = "figure/04_metacell/boxplot_dropout.pdf"
common_barplot(input = input, output = name, width = 7, height = 5)


# integration =============================================================

megexpmc <- allexpmc[[1]]
megexpmc$sample <- names(allexpmc)[1]

for (i in 2:length(allexpmc)) {
  data <- allexpmc[[i]]
  
  data$sample <- names(allexpmc)[i]
  megexpmc <- merge(megexpmc, data)
}

megexpmc <- autocluster_exp(megexpmc, nfeatures = 2000,
                            ndim = 20, neigh = 20, dist = 0.5, res = 0.5)

plot <- dimplot_new(megexpmc, pt.size = 0.75, label = T,
                    group.by = c("sample"))

name = "figure/04_metacell/umap_cluster_exp_meg.png"
ggsave(name, plot, dpi = 300, width = 8.5, height = 7)

megpromc <- allpromc[[1]]
megpromc$sample <- names(allpromc)[1]

for (i in 2:length(allpromc)) {
  data <- allpromc[[i]]
  
  data$sample <- names(allpromc)[i]
  megpromc <- merge(megpromc, data)
}

megpromc <- autocluster_pro(megpromc, ndim = 5,
                            features = rownames(megpromc),
                            neigh = 20, dist = 1, res = 0.5)

plot <- dimplot_new(megpromc, pt.size = 0.75, label = T,
                    group.by = c("sample"))

name = "figure/04_metacell/umap_cluster_pro_meg.png"
ggsave(name, plot, dpi = 300, width = 8.5, height = 7)

save(megexpmc, megpromc, file = "savedata/data2_meg.Rdata")

