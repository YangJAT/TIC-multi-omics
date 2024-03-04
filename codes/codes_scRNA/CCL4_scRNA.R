library(Seurat)
library(plyr)
library(dplyr)
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

path = "/home/ug0302/CITEseq3/07_scCCL4/data"
data <- read10x_all(path = path, merge = T, gene.column = 1)

name = "figure/01_basic/QC_plot.pdf"
pdf(file = name, width = 10, height = 7)
qc_plot(data = data, species = "human")
dev.off()

datafilt <- qc_filter(data = data,
                      nFeature_dn = 200,
                      nFeature_up = 10000,
                      nCount_dn = 500,
                      nCount_up = 50000,
                      mito_up = 20,
                      ribo_dn = 1,
                      ribo_up = 50,
                      gene_filter = 3,
                      species = "human")

datafilt <- autocluster(datafilt, nfeatures = 2000,
                        ndim = 20, neigh = 20,
                        dist = 0.5, res = 3)



# CytoTRACE ================================================================

library(CytoTRACE)
library(tibble)

input <- GetAssayData(datafilt, slot = "counts",
                      assay = "RNA") %>% as.matrix()

score <- CytoTRACE(input, enableFast = TRUE, ncores = 50)
score <- data.frame(id = names(score[["CytoTRACE"]]),
                    cytotrace = unname(score[["CytoTRACE"]]))

score <- column_to_rownames(score, var = "id")
datafilt <- AddMetaData(datafilt, score)

saveRDS(datafilt, "data/datafilt.rds")
datafilt <- readRDS("data/datafilt.rds")



# Visualization ===================================================================

# umap

plot <- prop_density(datafilt = datamodi,
                     group = "sample",
                     coord = "umap")

name = "figure/02_select/prop_density.pdf"
ggsave(name, plot, width = 10, height = 5)

# cytotrace

info <- datamodi@meta.data
input <- data.frame(type = info$sample,
                    value = info$cytotrace)

ecdf_plot(input, ylim = c(0,1))

# ITGA6

exp_distribution_ecdf(datamodi,
                      assay = "RNA",
                      select = "ITGA6",
                      group = "sample",
                      method = "wilcox.test",
                      ncol = 1)
# PVR

exp_distribution_ecdf(datamodi,
                      assay = "RNA",
                      select = "PVR",
                      group = "sample",
                      method = "wilcox.test",
                      ncol = 1)


# DEG =====================================================================

data <- as.matrix(GetAssayData(datamodi, slot = "data", assay = "RNA"))
info <- data.frame(id = colnames(datamodi), type = datamodi$sample)

diff <- difflimma(data = data, info = info,
                  group1 = "tag3", group2 = "tag2")

diff <- diff[!(diff$id %in% c("BHtag2", "BHtag3", "BHtag4")),]

input <- data.frame(id = diff$id, value = diff$t)
path = "/home/ug0302/CITEseq/public_data/stem_sig.csv"
result <- gsea_analysis(input = input,
                        source = path, geneset = NULL)

plot <- gsea_barcode(input = input,
                     source = path, geneset = NULL,
                     select = "Sensitizer")

ggsave("figure/06_immune/gsea_plot_Sensitizer.pdf", plot, width = 9, height = 3)


# CCL4 signature ====================

bpdata <- readRDS("bpdata_imbrave_150.rds")
data <- t(get.exp(bpdata, state.or.type = "type", cell.name = "tumor"))
data <- log2(data + 1)

name = "ICB_HCC_cohorts_clinical.rds"
info <- readRDS(name)
info <- info[["imbrave_150"]]
info <- info[info$sample_time %in% "pre",]

info <- data.frame(id = info$id,
                   type = info$response2)
info <- na.omit(info)

# signature score

ccl4_sig <- list(up = top_n(diff, 100, logfc)$id)

scores <- sig_scores(data, method = "zscore",
                     source = ccl4_sig, geneset = NULL,
                     min.sz = 2,
                     max.sz = 1000)
scores <- data.frame(t(scores))
scores <- rownames_to_column(scores, var = "id")

input <- merge(scores, info, by = "id")
input <- input[,-1]
names(input)[1] <- "value"

plot <- common_dotbox(input, rotate = 45,
                      order = T, decreasing = T,
                      method = "wilcox.test",
                      pt.size = 2)

name = "figure/04_sig/boxplot_all.pdf"
ggsave(name, plot, width = 3.9, height = 5)


# metacell =================================================================

datameta_list <- lapply(unique(datamodi$sample), function(i){
  subfilt <- datamodi[,datamodi$sample == i]
  metacell_object(subfilt, resolution = 70,
                  cluster.name = paste0(i, "_metacell"))
})

datameta <- datameta_list[[1]]
for (i in 2:length(datameta_list)) {
  data <- datameta_list[[i]]
  datameta <- merge(datameta, data)
}

datameta <- autocluster(datameta, nfeatures = 2000,
                        ndim = 20, neigh = 20,
                        dist = 0.5, res = 0.3)

datameta$sample <- datameta$orig.ident

datameta1 <- datameta[,datameta$sample == "tag2"]
datameta2 <- datameta[,datameta$sample == "tag3"]

coord1 <- data.frame(Embeddings(datameta1, "umap"))
coord2 <- data.frame(Embeddings(datameta2, "umap"))

plot <- ggplot(coord1, aes(coord1[,1], coord1[,2])) +
  
  geom_density2d(size = 0.35, alpha = 0.1,
                 colour = "black", bins = 10, h = 3.5) + 
  geom_point(size = 2) + 
  
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black", size = 15),
        axis.title = element_text(colour = "black", size = 15),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 1))

name = "figure/03_metacell/umap_density_control.pdf"
ggsave(name, plot, dpi = 300, width = 5.5, height = 5)


plot <- dimplot_new(datameta,
                    reduction = "umap",
                    pt.size = 1, label = T,
                    group.by = c("sample"))

name = "figure/03_metacell/umap_sample.png"
ggsave(name, plot, dpi = 300, width = 6, height = 5)

# cytotrace

cytotrace <- cytotrace_metacell(metacell = datameta,
                                rawcell = datamodi,
                                enableFast = TRUE)

datameta <- AddMetaData(datameta, cytotrace)
info <- datameta@meta.data
input <- data.frame(type = info$sample,
                    value = info$cytotrace)

plot <- histogram_boxplot(input = input,
                          bins = 30,
                          outlier.rm = FALSE)

name = "figure/03_metacell/hist_boxplot_cytotrace.pdf"
ggsave(name, plot, dpi = 300, width = 5, height = 3)

# ITGA6

exp_distribution_ecdf(datameta,
                      assay = "RNA",
                      select = "ITGA6",
                      group = "sample",
                      method = "wilcox.test",
                      ncol = 1)

data <- GetAssayData(datameta, slot = "data", assay = "RNA")["ITGA6",]
data <- data.frame(as.matrix(data))

info <- datameta@meta.data
input <- data.frame(value = data[,1],
                    type = info[,"sample"])

plot <- histogram_boxplot(input = input,
                          bins = 30,
                          outlier.rm = FALSE)

name = "figure/03_metacell/hist_boxplot_ITGA6.pdf"
ggsave(name, plot, dpi = 300, width = 5, height = 3)


# PVR

exp_distribution_ecdf(datameta,
                      assay = "RNA",
                      select = "PVR",
                      group = "sample",
                      method = "wilcox.test",
                      ncol = 1)



# cosine similarity to CD49f ===================================================

CD49f_cosine <- do.call(rbind, lapply(names(allexpmc), function(i){
  
  library(proxy)
  library(Seurat)
  library(plyr)
  library(dplyr)
  
  subexp <- allexpmc[[i]]
  subpro <- allpromc[[i]]
  
  vargene <- FindVariableFeatures(subexp, selection.method = "vst",
                                  nfeatures = 2000) %>% VariableFeatures()
  
  subexp <- as.matrix(GetAssayData(subexp, slot = "data"))
  metadata <- as.matrix(GetAssayData(datameta, slot = "data", assay = "RNA"))
  
  common <- intersect(vargene, rownames(datameta))
  subexp <- subexp[common,]
  metadata <- metadata[common,]
  
  dist <- proxy::dist(x = metadata, y = subexp,
                      method = "cosine", by_rows = FALSE)
  dist <- 1-as.matrix(as.data.frame.matrix(dist))
  
  subpro <- t(as.matrix(GetAssayData(subpro, slot = "data")))
  test <- dist %*% subpro
  
  input <- data.frame(id = rownames(test), value = test[,"CD49f"])
  type <- data.frame(id = colnames(datameta),
                     type = datameta$sample)
  input <- merge(type, input, by = "id")
  input <- input[,-1]
  
  normlize <- function(x){(x - min(x))/(max(x)-min(x))}
  input$value <- normlize(input$value)
  
  input <- input  %>%
    group_by(type) %>%
    dplyr::summarise(cor = mean(value))
  
  input$id <- i
  input
}))

input <- acast(CD49f_cosine, id ~ type, value.var = "cor",
             fun.aggregate = mean, na.rm = TRUE)

input <- data.frame(input)
input <- rownames_to_column(input, var = "id")

plot <- ggpaired(input, cond1 = "tag2", cond2 = "tag3",
         fill = "condition", palette = "jco", point.size = 3)+
  stat_compare_means(method = "t.test",paired = TRUE)

name = "figure/05_mapping/boxplot_CD49f_cosine.pdf"
ggsave(name, plot, width = 3.5, height = 4)



# mapping visualization ================================================================

library(Nebulosa)

allplots <- lapply(names(allexpmc), function(i){
  subpro <- allpromc[[i]]
  plot <- plot_density(subpro, "CD49f", size = 3)
  plot + ggtitle(i)
})

allplots <- plot_grid(plotlist = allplots, ncol = 5)
name = "figure/05_mapping/CD49f_exp.pdf"
ggsave(name, allplots, width = 27, height = 15)

subpro <- allpromc[["PLCPR5"]]
plot <- plot_density(subpro, "CD49f", size = 3)

name = "figure/05_mapping/CD49f_exp_PLC.pdf"
ggsave(name, plot, width = 6, height = 5)

allplots <- lapply(names(allexpmc), function(i){
  
  subexp <- allexpmc[[i]]
  ccl4exp <- datameta[,datameta$sample == "tag2"]
  
  mapdata <- refmapping(refer = subexp, query = ccl4exp,
                        method = "rpca", k.weight = 10,
                        refer.label = "seurat_clusters",
                        refer.coord = "wnn")
  
  refer <- data.frame(Embeddings(subexp, "wnn"))
  query <- data.frame(coord_x = mapdata$coord_x,
                      coord_y = mapdata$coord_y)
  
  ggplot(query, aes(x = query[,1],y = query[,2])) +
    
    stat_density_2d(aes(fill = ..density..), geom = "raster", h = 3.5, contour = FALSE) + 
    geom_density2d(size = 0.1, colour = "#FDAF9199", bins = 15, h = 3) +
    scale_fill_viridis(option="viridis") + 
    
    geom_point(data = refer, aes(x = refer[,1],y = refer[,2]),
               size = 2, color = "white", alpha = 0.3) + 
    
    theme_bw() + 
    labs(x = 'UMAP1',y= 'UMAP2',title = '') + 
    theme(strip.background = element_rect(colour = NA, fill = 'grey90'),
          strip.text.x = element_text(size = 15), 
          panel.grid = element_blank(),
          strip.placement = 'outside',
          axis.ticks = element_blank(),
          axis.text = element_blank()) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) + 
    ggtitle(i)
})

allplots <- plot_grid(plotlist = allplots, ncol = 5)
name = "figure/05_mapping/mapping_ccl4.pdf"
ggsave(name, allplots, width = 27, height = 15)

subexp <- allexpmc[["PLCPR5"]]
ccl4exp <- datameta[,datameta$sample == "tag3"]

mapdata <- refmapping(refer = subexp, query = ccl4exp,
                      method = "rpca", k.weight = 10,
                      refer.label = "seurat_clusters",
                      refer.coord = "wnn")

refer <- data.frame(Embeddings(subexp, "wnn"))
query <- data.frame(coord_x = mapdata$coord_x,
                    coord_y = mapdata$coord_y)

plot <- ggplot(query, aes(x = query[,1],y = query[,2])) +
  
  stat_density_2d(aes(fill = ..density..), geom = "raster", h = 3.5, contour = FALSE) + 
  geom_density2d(size = 0.1, colour = "#FDAF9199", bins = 15, h = 3) +
  scale_fill_viridis(option="viridis") + 
  
  geom_point(data = refer, aes(x = refer[,1],y = refer[,2]),
             size = 2, color = "white", alpha = 0.3) + 
  
  theme_bw() + 
  labs(x = 'UMAP1',y= 'UMAP2',title = '') + 
  theme(strip.background = element_rect(colour = NA, fill = 'grey90'),
        strip.text.x = element_text(size = 15), 
        panel.grid = element_blank(),
        strip.placement = 'outside',
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))

name = "figure/05_mapping/mapping_ccl4_PLC.pdf"
ggsave(name, plot, width = 5.7, height = 5)
