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
library(SingleCellExperiment)

source("citeseq_function.R")


# TotalVI  =================================================================

select_dm <- c("CD13", "CD24", "CD44", "CD47", "CD49f", "CD54",
               "CD73", "CD133", "CD326", "CD338")

alltotal = list()
for (i in names(allexp)) {
  
  expdata <- allexp[[i]]; prodata <- allpro[[i]]
  prodata <- prodata[select_dm,]
  
  total <- totalvi_latent(expdata = expdata, prodata = prodata)
  alltotal <- c(alltotal, list(total))
}

names(alltotal) <- names(allexp)
save(alltotal, file = "savedata/totalvi_results.Rdata")


# integration =====================================================================

allexp <- lapply(names(allexp), function(i){
  
  subexp <- allexp[[i]]; total <- alltotal[[i]]
  subexp[["total"]] <- CreateDimReducObject(embeddings = total, key = "total")
  
  subexp <- RunUMAP(subexp, dims = 1:20,
                    n.neighbors = 20, min.dist = 1, 
                    reduction = "total", reduction.name = "umap_tvi")
  
  subexp <- FindNeighbors(subexp, dims = 1:20, reduction = "total")
  subexp <- FindClusters(subexp, resolution = 0.3, n.iter = 50)
  
  subexp$cluster_tvi <- subexp$seurat_clusters
  subexp
})

names(allexp) <- names(allpro)

# clustering --------------------

all_plots <- lapply(names(allexp), function(i){
  dimplot_new(allexp[[i]], pt.size = 0.75, reduction = "umap_tvi",
              label = T, group.by = c("cluster_tvi")) + ggtitle(i)
})

name = "figure/02_integ/cluster_tvi_sep.png"
export_dimplot(all_plots = all_plots, ncol = 2,
               dpi = 300, output = name)


