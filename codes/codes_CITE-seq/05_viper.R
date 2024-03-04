library(Seurat)
library(plyr)
library(future)
library(tibble)
library(reshape2)
library(tidyverse)
library(ggalluvial)
library(ggbeeswarm)
library(data.table)
library(future.apply)
library(RColorBrewer)
library(cowplot)
library(ggpubr)
library(ggplot2)
library(cowplot)
library(ComplexHeatmap)
source("citeseq_function.R")

select_pro <- c("CD13", "CD24", "CD44", "CD47", "CD49f", "CD54",
                "CD73", "CD90", "CD95", "CD117", "CD133", "CD155",
                "CD183", "CD184", "CD274", "CD326", "CD338")

select_gene <- c("ANPEP", "CD24", "CD44", "THY1", "PROM1", "EPCAM",
                 "NT5E", "CD47", "ICAM1", "ABCG2", "ITGA6", "KIT",
                 "CD274", "PVR", "CXCR3", "FAS", "CXCR4")


# gene - protein ----------

anno_gene <- c("ANPEP" = "CD13",
               "CD24" = "CD24",
               "CD44" = "CD44",
               "THY1" = "CD90",
               "PROM1" = "CD133",
               "EPCAM" = "CD326",
               "NT5E" = "CD73",
               "CD47" = "CD47",
               "ICAM1" = "CD54",
               "ABCG2" = "CD338",
               "ITGA6" = "CD49f",
               "KIT" = "CD117",
               "CD274" = "CD274",
               "PVR" = "CD155",
               "CXCR3" = "CD183",
               "FAS" = "CD95",
               "CXCR4" = "CD184")

anno_gene <- data.frame(gene = names(anno_gene),
                        protein = unname(anno_gene))



# ARACNe input ===========================================================

allinput <- lapply(allexp, function(i){
  metacell_aracne(i, resolution = 35,
                  filter.gene = TRUE, cluster.name = "metacell")
})

# savedata

for (i in names(allinput)) {
  
  data <- allinput[[i]]
  
  if (ncol(data) >= 270) {
    data <- data[,sample(colnames(data), size = 270, replace = F)]
  } else {data <- data}
  
  dir.create(paste0("/home/ug0302/CITEseq/savedata/aracne/", i))
  name <- paste0("/home/ug0302/CITEseq/savedata/aracne/", i, "/", i, ".txt")
  
  write.table(data, name, sep = "\t", quote = F,
              row.names = T, col.names = NA)
}


# process network files ============================================================

filename <-  list.files("savedata/aracne/01_net_raw")

for (i in filename) {
  
  name <- paste0("savedata/aracne/01_net_raw/", i)
  network <- read.table(name, sep = "\t", header = T)
  
  network <- network[,c(1:3)]
  output_name <- paste0("savedata/aracne/02_net_norm/", i)
  
  write.table(network, output_name, sep = "\t",
              quote = F, row.names = F, col.names = F)
}


# generate regulon list =============================================================

filename <-  list.files("savedata/aracne/02_net_norm")
regulon_list <- c()

for (i in filename) {
  
  name <- unlist(strsplit(i, "_"))[1]
  exp_name <- paste0("savedata/aracne/", name, "/", name, ".txt")
  net_name <- paste0("savedata/aracne/02_net_norm/", i)
  
  regulon <- aracne2regulon(afile = net_name, exp_name, format = "3col")
  regulon <- list(regulon); names(regulon) <- name
  regulon_list <- c(regulon_list, regulon)
}

saveRDS(regulon_list, "savedata/aracne/regulon_list.rds")



# viper =================================================================

library(viper)

source("viper_function.R")

plan(multisession, workers = length(allexp))
options(future.globals.maxSize = 20000*1024^2)

allvip <- future_lapply(names(allexp), function(i){
  
  data <- GetAssayData(allexp[[i]], slot = "data", assay = "RNA") %>% as.matrix()
  data <- RankTransform(data)
  
  # Meta VIPER
  
  # result <- viper(data, regulon_list_sub,
  #                 minsize = 25, cores = 10, method = 'none')
  
  result <- viper(data, regulon_list[[i]],
                  minsize = 20, cores = 10, method = 'none')
  
  CreateSeuratObject(result)
}, future.seed = TRUE)

plan(sequential)
names(allvip) <- names(allexp)


# metacell - viper --------------------

allvipmc <- lapply(names(allexpmc), function(i){
  
  subexpmc <- allexpmc[[i]]
  subvip <- allvip[[i]]
  
  info <- subexpmc@meta.data
  info <- data.frame(metacell = rownames(info),
                     cell_id = info$cell_id)
  
  info <- do.call(rbind, lapply(info$metacell, function(i){
    cell_id <- unlist(strsplit(info$cell_id[info$metacell == i], ','))
    data.frame(id = cell_id, metacell = i)
  }))
  
  info <- column_to_rownames(info, var = "id")
  subvip <- AddMetaData(subvip, info)
  
  subvipmc <- AverageExpression(subvip, group.by = "metacell",
                                slot = "counts", return.seurat = FALSE)[[1]]
  CreateSeuratObject(subvipmc)
})

names(allvipmc) <- names(allexpmc)

# savedata --------------------

name = "/home/ug0302/CITEseq/savedata/viper_results.Rdata"
save(allvip, allvipmc, file = name)


# gene protein viper ================================================

nulldata <- lapply(names(allexp), function(i){
  
  select <- intersect(select_gene, rownames(allexp[[i]]))
  all_plots <- featureplot_new(data = allexp[[i]],
                               pt.size = 0.75, reduction = "wnn",
                               color = "parula", features = select)
  
  name = paste0("figure/03_viper/wnn_feature_expression_", i, ".png")
  export_featureplot(all_plots = all_plots,
                     ncol = 5, dpi = 300, output = name)
})


allpro <- lapply(names(allpro), function(i){
  
  subexp <- allexp[[i]]
  subpro <- allpro[[i]]
  
  wnn <- Embeddings(subexp, "wnn")
  subpro[["wnn"]] <- CreateDimReducObject(wnn, key = "wnnUMAP")
  subpro
})

names(allpro) <- names(allexp)


nulldata <- lapply(names(allpro), function(i){
  
  select <- intersect(select_pro, rownames(allpro[[i]]))
  all_plots <- featureplot_new(data = allpro[[i]],
                               pt.size = 0.75, reduction = "wnn",
                               color = "parula", features = select)
  
  name = paste0("figure/03_viper/wnn_feature_proteins_", i, ".png")
  export_featureplot(all_plots = all_plots,
                     ncol = 5, dpi = 300, output = name)
})


# viper visualization --------------------

allvip <- lapply(names(allvip), function(i){
  
  subexp <- allexp[[i]]
  subpro <- allvip[[i]]
  
  wnn <- Embeddings(subexp, "wnn")
  subpro[["wnn"]] <- CreateDimReducObject(wnn, key = "wnnUMAP")
  subpro
})

names(allvip) <- names(allexp)

nulldata <- lapply(names(allvip), function(i){
  
  select <- intersect(select_gene, rownames(allvip[[i]]))
  all_plots <- featureplot_new(data = allvip[[i]],
                               pt.size = 0.75, reduction = "wnn",
                               color = "parula", features = select)
  
  name = paste0("figure/03_viper/wnn_feature_viper_", i, ".png")
  export_featureplot(all_plots = all_plots,
                     ncol = 5, dpi = 300, output = name)
})


# protein-viper consistency ====================================================

allcor <- lapply(names(allvip), function(i){
  
  subpro <- allpro[[i]]
  subvip <- allvip[[i]]
  
  subpro <- t(as.matrix(GetAssayData(subpro, slot = "data", assay = "RNA")))
  subvip <- t(as.matrix(GetAssayData(subvip, slot = "data", assay = "RNA")))
  
  common <- (anno_gene$protein %in% colnames(subpro)) &
            (anno_gene$gene %in% colnames(subvip))
  
  subpro <- subpro[,anno_gene$protein[common]]
  subvip <- subvip[,anno_gene$gene[common]]
  cor_betweenAB(subpro, subvip)
})

names(allcor) <- names(allpro)


# visualization --------------------

for (i in names(allcor)) {
  input <- allcor[[i]]
  
  name = paste0("figure/03_viper/heatmap_cor_", i, ".pdf")
  pdf(file = name, width = 8, height = 7)
  print(heatmap_text(input = input, col_rot = 45, 
                     color = "blue2red", cutoff = 0.05,
                     order_name = F, cluster_row = F, cluster_col = F))
  dev.off()
}
