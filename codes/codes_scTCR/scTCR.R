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
library(RColorBrewer)
source("public_function_mouse.R")

path = "/home/ug0302/CITEseq3/04_scTCR/data/CD49f_H1"
tcr1 <- tcr_basic(path)
rownames(tcr1) <- paste0("CD49f_H1_", rownames(tcr1))

path = "/home/ug0302/CITEseq3/04_scTCR/data/CD49f_L1"
tcr2 <- tcr_basic(path)
rownames(tcr2) <- paste0("CD49f_L1_", rownames(tcr2))

tcr <- rbind(tcr1, tcr2)
datafilt <- AddMetaData(datafilt, tcr)

CD49H_CD8 <- datafilt[,datafilt$group == "H" & datafilt$celltype_sig2 == "CD8T"]
CD49H_CD8 <- na.omit(data.frame(clonotype_id = CD49H_CD8$clonotype_id))

CD49L_CD8 <- datafilt[,datafilt$group == "L" & datafilt$celltype_sig2 == "CD8T"]
CD49L_CD8 <- na.omit(data.frame(clonotype_id = CD49L_CD8$clonotype_id))

CD49H_CD4 <- datafilt[,datafilt$group == "H" & datafilt$celltype_sig2 == "CD4T"]
CD49H_CD4 <- na.omit(data.frame(clonotype_id = CD49H_CD4$clonotype_id))

CD49L_CD4 <- datafilt[,datafilt$group == "L" & datafilt$celltype_sig2 == "CD4T"]
CD49L_CD4 <- na.omit(data.frame(clonotype_id = CD49L_CD4$clonotype_id))

clone1 <- clone_size(CD49H_CD8)
clone2 <- clone_size(CD49L_CD8)
clone3 <- clone_size(CD49H_CD4)
clone4 <- clone_size(CD49L_CD4)

stat <- data.frame(table(clone1$clone_type))
plot <- circ_prop(stat)
ggsave("figure/CD49H_CD8.pdf", plot, width = 5, height = 5)

stat <- data.frame(table(clone2$clone_type))
plot <- circ_prop(stat)
ggsave("figure/CD49L_CD8.pdf", plot, width = 5, height = 5)

stat <- data.frame(table(clone3$clone_type))
plot <- circ_prop(stat)
ggsave("figure/CD49H_CD4.pdf", plot, width = 5, height = 5)

stat <- data.frame(table(clone4$clone_type))
plot <- circ_prop(stat)
ggsave("figure/CD49L_CD4.pdf", plot, width = 5, height = 5)
