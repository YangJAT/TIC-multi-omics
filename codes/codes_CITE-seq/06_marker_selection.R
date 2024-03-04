library(NMF)
library(Seurat)
library(plyr)
library(dplyr)
library(scater)
library(future)
library(gtools)
library(stringr)
library(ggplot2)
library(cowplot)
library(data.table)
library(tidyverse)
library(future.apply)
library(RColorBrewer)
source("citeseq_function.R")


# 1 cytotrace ==================================================================

plan(multisession, workers = length(allexp))
options(future.globals.maxSize = 20000*1024^2)

cytotrace <- future_lapply(names(allexp), function(i){
  
  subexp <- allexp[[i]]
  subexpmc <- allexpmc[[i]]
  cytotrace_metacell(metacell = subexpmc, rawcell = subexp)
  
}, future.seed = TRUE)
plan(sequential)
names(cytotrace) <- names(allexp)


# protein - cytotrace relationship --------------------

allcor <- do.call(rbind, lapply(names(allpromc), function(i){
  
  info <- cytotrace[[i]]
  subpro <- allpromc[[i]]
  
  subpro <- t(as.matrix(GetAssayData(subpro, slot = "data")))
  subpro <- subpro[rownames(info),]
  
  cor <- cor_between_aB(info, subpro)
  cor$sample <- i
  cor
}))

allcor <- data.frame(type = allcor$id,
                     value = allcor$cor,
                     sample = allcor$sample)
# visualization

name = "figure/06_optimal/protein_cytotrace.pdf"
jitter_boxplot2(input = allcor,
                start.point = 0,
                name.coord = "Value",
                order = TRUE, decreasing = F, color = col,
                output = name, width = 6, height = 7.5)


# Cell entropy ===============================================================

plan(multisession, workers = length(allexp))
options(future.globals.maxSize = 20000*1024^2)

entropy <- future_lapply(names(allexp), function(i){
  
  subexp <- allexp[[i]]
  subexpmc <- allexpmc[[i]]
  entropy_metacell(metacell = subexpmc, rawcell = subexp)
  
}, future.seed = TRUE)
plan(sequential)
names(entropy) <- names(allexp)


# protein - entropy relationship --------------------

allcor <- do.call(rbind, lapply(names(entropy), function(i){
  
  info <- entropy[[i]]
  subpro <- allpromc[[i]]
  
  subpro <- t(as.matrix(GetAssayData(subpro, slot = "data")))
  subpro <- subpro[rownames(info),]
  
  cor <- cor_between_aB(info, subpro)
  cor$sample <- i
  cor
}))

allcor <- data.frame(type = allcor$id,
                     value = allcor$cor,
                     sample = allcor$sample)
# visualization

name = "figure/06_optimal/protein_entropy.pdf"
jitter_boxplot2(input = allcor,
                start.point = 0,
                name.coord = "Value",
                order = TRUE, decreasing = F, color = col,
                output = name, width = 6, height = 7.5)


# protein - TF (viper)  =========================================

select <- read.table("public_data/CSC_TF.txt", sep = "\t", header = T)
select <- select[,1]

allcor <- lapply(names(allpromc), function(i){
  
  subvip <- allvipmc[[i]]
  subvip <- t(as.matrix(GetAssayData(subvip, slot = "data")))
  subvip <- subvip[,colnames(subvip) %in% select]
  
  subpro <- allpromc[[i]]
  subpro <- t(GetAssayData(subpro, slot = "data"))
  
  limma_betweenAB(subvip, subpro,
                  value = "es",
                  format = "long",
                  thres = 0.75)
  
  # cor_betweenAB_exact(subvip, subpro,
  #                     format = c("long"),
  #                     method = c("spearman"))
})

names(allcor) <- names(allpromc)

# annotate marker

marker = c("SALL4_CD49f", "SOX9_CD49f", "SOX9_CD44", "NOTCH3_CD24",
           "MYC_CD44", "HMGA2_CD49f", "SOX2_CD49f", "NOTCH3_CD155",
           "SALL4_CD54", "STAT3_CD155", "STAT3_CD49f", "MYC_CD49f")
select_gene <- input[input$id %in% marker,]

# visualization

plot <- ggplot(data = input, aes(logfc, -log10(pvalue))) + 
  geom_point(alpha = 1, size = 6, aes(shape = type1, fill = type2, color = type2)) + 
  scale_fill_manual(values = c(col)) + 
  scale_color_manual(values = c(col)) + 
  scale_shape_manual(values=c(0, 1, 2, 5, 21, 22, 23, 24)) + 
  
  col <- c("#B8E3EA", "#5CB3DA", "#0070B2", "#FBDD7E", "#F7AE24", "#FF7149", 
           "#F2D7EE", "#A37CB7", "#A231A1", "#ECB2C8", "#E93B8C", "#B91372", 
           "#FF9F99", "#F15E4C", "#DA1735", "#CDE391", "#8BBE53", "#679436", 
           "#98D4C6", "#00A385", "#067D69", "#B2DBBF", "#028090", "#114B5F", 
           "#FBD0C0", "#CD6981", "#A23E48", "#CCDAD1", "#9CAEA9", "#788585")
  
  geom_hline(yintercept = -log10(0.05), color = "black", linetype = "dashed", lwd = 0.75) + 
  
  geom_text_repel(data = select_gene, aes(label=id), size = 2) + 
  
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black", size = 15),
        axis.title = element_text(colour = "black", size = 15),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 1)) + 
  
  scale_y_continuous(limits = c(10, 320),  
                     breaks = c(10, 20, 40, 80, 160, 320),
                     trans = "log10") + 
  # scale_y_continuous(limits = c(0, 10),  
  #                    breaks = c(0:10)) + 
  
  scale_x_continuous(limits = c(-0.7, 0.7),  
                     breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6)) + 
  labs(x = 'Log Fold Change',y= 'Log (P-value)',title = '')

name = "figure/06_optimal/volcano_cor_vip_1.pdf"
ggsave(name, plot, width = 10, height = 6.5)

name = "figure/06_optimal/volcano_cor_vip_input.txt"
write.table(input, name, sep = "\t", quote = F,
            row.names = F, col.names = T)



# functional analysis for each marker =======================================================

select_al <- c("CD13", "CD24", "CD44", "CD47", "CD49f", "CD54",
               "CD73", "CD90", "CD117", "CD133", "CD326", "CD338",
               "CD183", "CD184", 
               "CD274", "CD155",
               "IgG1", "IgG2a", "IgG2b")

# DEA --------------------

plan(multisession, workers = length(allpromc))
options(future.globals.maxSize = 3000*1024^2)

alldiff <- future_lapply(names(allpromc), function(i){
  
  subexp <- allexpmc[[i]]
  subpro <- allpromc[[i]]
  
  subexp <- t(as.matrix(GetAssayData(subexp, slot = "data")))
  subpro <- t(as.matrix(GetAssayData(subpro, slot = "data")))[,select_al]
  
  limma_betweenAB(mat.value = subexp,
                  mat.type = subpro,
                  thres = 0.9,
                  value = "es")
  
}, future.seed = TRUE)
plan(sequential)
names(alldiff) <- names(allpromc)

# integrate --------------------

input <- do.call(rbind, lapply(names(alldiff), function(i){
  data <- data.frame(alldiff[[i]])
  data <- rownames_to_column(data, var = "id")
  reshape2::melt(data, id.vars = c("id"))
}))

input <- acast(input, id ~ variable,
               value.var = "value",
               fun.aggregate = mean_rmout)

save(diff_sig, file = "figure/06_optimal/diff_sig.Rdata")


# KnockTF --------------------

load("/home/ug0302/CITEseq/public_data/KnockTF.Rdata")
name = "/home/ug0302/CITEseq/public_data/CSC_marker_49.txt"
info <- fread(file = name, header = T, sep = '\t', fill = T, data.table = F)

select <- c("MYC_DataSet_01_237", "MYCN_DataSet_01_160",
            "KLF4_DataSet_01_341", "STAT3_DataSet_01_180",
            "SOX9_DataSet_01_308", "SOX2_DataSet_01_349",
            "SALL4_DataSet_01_277", "NANOG_DataSet_01_94",
            "POU5F1_DataSet_01_05")

pbsig <- pertubsig[str_detect(names(pertubsig),
                              paste(select, collapse = '|'))]

pbscores <- pertub_score(mat.fc = input,
                         pb.sig = pbsig,
                         nperm = 1000)

pbscores_flt <- pbscores[pbscores$pvalue < 0.05,]

pbinput <- acast(pbscores, id ~ type,
                 value.var = "zscore", fun.aggregate = mean)

name = "figure/06_optimal/heatmap_KnockTF.pdf"
pdf(file = name, width = 9, height = 5)

heatmap_text.perb(input = pbinput, col_rot = 90, 
                  cutoff = 2, order_name = F,
                  cluster_row = T, cluster_col = T,
                  color = "parula",
                  outlier.rm = T)
dev.off()

