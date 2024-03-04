library(Seurat)
library(reshape2)
library(cowplot)

setwd("/home/ug0302/CITEseq2")
source("/home/ug0302/CITEseq/public_data/citeseq_function.R")


# data input =====================================================================

name <- "/home/ug0302/CITEseq2/data/HCC_from_tumor_Nature.rds"
datafilt <- readRDS(name)

name = "/home/ug0302/CITEseq2/data/LRanalysis/NicheNet_LR_pairs.csv"
lrinfo <- read.table(name, sep = ",", header = T, check.names = F)

all_genesets <- readRDS("/home/ug0302/CITEseq/public_data/all_genesets.rds")
secreted = all_genesets[["msigdb_all"]][["NABA_SECRETED_FACTORS"]]

ligand <- intersect(unique(lrinfo$from), rownames(datafilt))
ligand <- intersect(ligand, secreted)

load("/home/ug0302/CITEseq2/data/LRanalysis/Science.Rdata")
membrane = protein_anno[["Science_protein_class"]][["Predicted membrane proteins"]]

receptor <- intersect(unique(lrinfo$to), rownames(datafilt))
receptor <- intersect(receptor, membrane)


# JSD analysis ======================================================================

data_l <- GetAssayData(datafilt, slot = "data", assay = "RNA")[ligand,]
data_l <- as.matrix(data_l)

data_r <- GetAssayData(datafilt, slot = "data", assay = "RNA")[receptor,]
data_r <- as.matrix(data_r)

cluster <- data.frame(id = colnames(datafilt), 
                      type = datafilt$celltype_sig2)

jsd_l <- jsd_score(data = data_l,
                   cluster = cluster,
                   select = "Neutrophil")

jsd_r <- jsd_score(data = data_r,
                   cluster = cluster,
                   select = "Neutrophil")

exp_l <- meanexpression(datafilt[ligand,], slot = "data",
                        group = "celltype_sig2", method = "mean")

exp_l <- data.frame(id = rownames(exp_l), exp = exp_l[,"Neutrophil"])
input_l <- merge(exp_l, jsd_l, by = "id")

plot <- scatter_highlight(input_l, topn = 5, color = "#0070b2")
name = "figure/LRanalysis/ligands.pdf"
ggsave(name, plot, width = 6, height = 5.5)

exp_r <- meanexpression(datafilt[receptor,], slot = "data",
                        group = "celltype_sig2", method = "mean")

exp_r <- data.frame(id = rownames(exp_r), exp = exp_r[,"Neutrophil"])
input_r <- merge(exp_r, jsd_r, by = "id")

plot <- scatter_highlight(input_r, topn = 5, color = "#0070b2")
name = "figure/LRanalysis/receptor.pdf"
ggsave(name, plot, width = 6, height = 5.5)



# ligand - cell type abundance ========================================================

select = c("CXCL6", "CXCL8", "CXCL1", "CXCL2",
           "CXCL3", "CXCL5", "CXCL7")

subfilt <- datafilt[,datafilt$celltype_sig2 == "tumor"]
data <- meanexpression(subfilt[select,], slot = "data", group = "sample")
data <- t(as.matrix(data))

cellprop <- data.frame(id = colnames(datafilt),
                       sample = datafilt$sample,
                       type = datafilt$celltype_sig2)

cellprop <- cellprop[!(cellprop$type %in% c("tumor", "Endothelial", "Fibroblast")),]

cellprop <- data.frame(table(cellprop$sample, cellprop$type))
cellprop <- acast(cellprop, Var1 ~ Var2, value.var = "Freq",
                  fun.aggregate = mean, na.rm = TRUE)
cellprop <- proportions(cellprop, 1)

cor_neu <- do.call(rbind, lapply(colnames(data), function(i){
  input1 <- data.frame(id = rownames(data), value = data[,i])
  input2 <- data.frame(id = rownames(cellprop), prop = cellprop[,"Neutrophil"])
  input <- merge(input1, input2, by = "id")
  cor <- cor.test(input[,2], input[,3], method = "pearson")
  data.frame(id = i, cor = cor$estimate, pvalue = cor$p.value)
}))

input <- data.frame(id = cor_neu$id,
                    value = cor_neu$cor)

plot <- loading_plot2(input = input,
                      color = "#560047",
                      topn = NULL,
                      marker = select,
                      decreasing = T)

name = "figure/LRanalysis/cor_ligand_content.pdf"
ggsave(name, plot, width = 6.5, height = 5)

name = "figure/LRanalysis/cor_ligand_content.txt"
write.table(cor_neu, name, sep = "\t", quote = F,
            row.names = F, col.names = T)


# ligands in neutrophils - ITGA6 in tumor cells ==========================

name = "/home/ug0302/CITEseq2/data/LRanalysis/NicheNet_LR_pairs.csv"
lrinfo <- read.table(name, sep = ",", header = T, check.names = F)

all_genesets <- readRDS("/home/ug0302/CITEseq/public_data/all_genesets.rds")
secreted = all_genesets[["msigdb_all"]][["NABA_SECRETED_FACTORS"]]

ligand <- intersect(unique(lrinfo$from), rownames(datafilt))
ligand <- intersect(ligand, secreted)

subfilt <- datafilt[,datafilt$celltype_sig2 == "tumor"]
data <- meanexpression(subfilt, slot = "data", group = "sample")
data <- t(as.matrix(data))

select = "ITGA6"
value1 <- data.frame(id = rownames(data),
                     value1 = data[,select])

subfilt <- datafilt[,datafilt$celltype_sig2 == "Neutrophil"]
data <- meanexpression(subfilt, slot = "data", group = "sample")
data <- t(as.matrix(data))

cor_neu <- do.call(rbind, lapply(ligand, function(i){
  value2 <- data.frame(id = rownames(data),
                       value2 = data[,i])
  
  input <- merge(value1, value2, by = "id")
  cor <- cor.test(input[,2], input[,3], method = "spearman")
  data.frame(id = i, cor = cor$estimate, pvalue = cor$p.value)
}))

cor_neu <- na.omit(cor_neu)
cor_neu$FDR <- p.adjust(cor_neu$pvalue, method = "fdr")

select = c("CCL4", "CXCL8")

input <- data.frame(id = cor_neu$id,
                    value = cor_neu$cor)

plot <- loading_plot2(input = input,
                      color = "#560047",
                      topn = 5,
                      marker = select,
                      decreasing = T)

name = "figure/LRanalysis/loadingplot_ligand_CD49f.pdf"
ggsave(name, plot, width = 5, height = 7)



# CCL4 in different cell types - ITGA6 in tumor cells ======================

subfilt <- datafilt[,datafilt$celltype_sig2 == "tumor"]
data <- meanexpression(subfilt, slot = "data", group = "sample")
data <- t(as.matrix(data))

select = "ITGA6"
value1 <- data.frame(id = rownames(data),
                     ITGA6 = data[,select])

select <- unique(datafilt$celltype_sig2)
select <- select[!(select == "tumor")]

allplots <- lapply(select, function(i){
  
  subfilt <- datafilt[,datafilt$celltype_sig2 == i]
  data <- meanexpression(subfilt, slot = "data", group = "sample")
  data <- t(as.matrix(data))
  
  value2 <- data.frame(id = rownames(data),
                       CCL4 = data[,"CCL4"])
  
  input <- merge(value1, value2, by = "id")
  input <- input[,-1]
  plot <- corplot(input, method = "spearman")
  plot + ggtitle(i)
})

allplots <- plot_grid(plotlist = allplots, ncol = 4)
name = "figure/LRanalysis/cor_ITGA6_CCL4.pdf"
ggsave(name, allplots, width = 15, height = 12)


# CCL4 in neutrophils between CD49f-hi and CD49f-lo ----------

datafilt <- seurat_sample_group(datafilt = datafilt,
                                col.sample = "sample",
                                col.type = "celltype_sig2",
                                col.type.cell = "tumor",
                                gene.select = "ITGA6",
                                gene.metrics = "exp",
                                gene.thres = 0.75)

datafilt_flt <- datafilt[,datafilt$group %in% c("H", "L")]

type = "Neutrophil"
gene = "CCL4"

input <- datafilt_flt[,datafilt_flt$celltype_sig2 == type]
plot <- exp_distribution_ecdf(datafilt = input,
                              select = gene,
                              group = "group",
                              method = "wilcox.test",
                              ncol = 1)

name = paste0("figure/LRanalysis/Distri_exp_", type, "_", gene, ".pdf")
ggsave(name, plot, width = 5.5, height = 5)


