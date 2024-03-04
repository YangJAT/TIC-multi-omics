library(NMF)
library(Seurat)
library(plyr)
library(dplyr)
library(future)
library(gtools)
library(ggplot2)
library(cowplot)
library(data.table)
library(tidyverse)
library(future.apply)
library(RColorBrewer)
source("citeseq_function.R")


# NMF ======================================================================

# generate input data --------------------

nmf_input <- lapply(allexpmc, function(i){
  
  data <- NormalizeData(i)
  data <- FindVariableFeatures(data, nfeatures = 7000)
  data <- ScaleData(data, do.center = F)
  
  as.matrix(GetAssayData(data,
            slot = "scale.data", assay = "RNA"))
})


# NMF program --------------------

allnmf_integ <- consensus_nmf(nmf_input,
                              topnumber = 150,
                              knumber = 5:10,
                              thres_freq = 2,
                              thres_cor = 0.8, 
                              ncore = length(nmf_input))

# savedata --------------------

name = "/home/ug0302/CITEseq/savedata/nmf_results.Rdata"
save(allnmf_integ, file = name)



# NMF - umap =============================================================

allnmf <- do.call(c, lapply(allnmf_integ, function(i){list(i$cell)}))
names(allnmf) <- names(allnmf_integ)


# program featureplot --------------------

nulldata <- lapply(names(allexpmc), function(i){
  
  subexp <- allexpmc[[i]]
  subnmf <- allnmf[[i]]
  
  subexp <- AddMetaData(subexp, subnmf)
  feature <- colnames(subnmf)
  
  all_plots <- featureplot_new(data = subexp, pt.size = 2,
                               reduction = "wnn", color = "blue2red",
                               features = feature)
  
  name = paste0("figure/05_nmf/wnn_nmf_", i, ".png")
  export_featureplot(all_plots = all_plots,
                     ncol = 5, dpi = 300, output = name)
})

all_plots <- lapply(names(allexpmc), function(i){
  
  subexp <- allexpmc[[i]]
  subnmf <- allnmf[[i]]
  
  nmfclass <- nmf_multi_cluster(subnmf, thres = 0.75)
  subexp <- AddMetaData(subexp, nmfclass)
  
  dimplot_new(subexp, pt.size = 1.5, reduction = "wnn",
              label = F, group.by = c("nmfclass")) + ggtitle(i)
})

name = "figure/05_nmf/wnn_nmfclass_sep.png"
export_dimplot(all_plots = all_plots,
               ncol = 2, dpi = 300, output = name)


# NMF - signature relationship =======================================================

allnmf <- do.call(c, lapply(allnmf_integ, function(i){list(i$cell)}))
names(allnmf) <- names(allnmf_integ)


# signature score --------------------

path = "/home/ug0302/CITEseq/public_data/tumor_sig5.csv"

allsig <- lapply(allexpmc, function(i){
  seurat_score(i, source = path)
})


# calculate --------------------

allcor <- lapply(names(allnmf), function(i){
  
  subsig <- allsig[[i]]
  subnmf <- allnmf[[i]]
  subnmf <- subnmf[rownames(subsig),]
  cor_betweenAB(subnmf, subsig)
})

names(allcor) <- names(allnmf)


# Visualization --------------------

for (i in names(allcor)) {
  input <- allcor[[i]]
  
  name = paste0("figure/05_nmf/heatmap_cor_nmf_sig_", i, ".pdf")
  pdf(file = name, width = 8, height = 7)
  print(heatmap_text(input = input, col_rot = 90, 
                     color = "blue2red", cutoff = 0.05,
                     order_name = F, cluster_row = T, cluster_col = F))
  dev.off()
}



# NMF - protein relationship =========================================================

allnmf <- do.call(c, lapply(allnmf_integ, function(i){list(i$cell)}))
names(allnmf) <- names(allnmf_integ)

allcor <- do.call(cbind, lapply(names(allnmf), function(i){
  
  subpro <- allpromc[[i]]
  subpro <- t(GetAssayData(subpro, slot = "data"))[,select_pro]
  
  subnmf <- allnmf[[i]]
  subnmf <- subnmf[rownames(subpro),]
  cor_betweenAB(subnmf, subpro)
}))

# Visualization ----------

pdf(file = "figure/05_nmf2/heatmap_allcor_nmf_pro.pdf",
    width = 7, height = 20)

heatmap_cor(t(allcor), color = "blue2red",
            rev.col = FALSE, show.name = TRUE)
dev.off()


allnmf <- do.call(c, lapply(allnmf_integ, function(i){list(i$cell)}))
names(allnmf) <- names(allnmf_integ)

allcor <- lapply(names(allnmf), function(i){
  
  subpro <- allpromc[[i]]
  subpro <- t(GetAssayData(subpro, slot = "data"))
  
  subnmf <- allnmf[[i]]
  subnmf <- subnmf[rownames(subpro),]
  cor_betweenAB(subnmf, subpro)
})

names(allcor) <- names(allnmf)

for (i in names(allcor)) {
  input <- allcor[[i]]
  
  name = paste0("figure/05_nmf/heatmap_cor_nmf_pro_", i, ".pdf")
  pdf(file = name, width = 5, height = 7)
  print(heatmap_text(input = input, col_rot = 90, 
                     color = "blue2red", cutoff = 0.05,
                     order_name = F, cluster_row = T, cluster_col = F))
  dev.off()
}

