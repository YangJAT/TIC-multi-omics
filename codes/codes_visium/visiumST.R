library(Seurat)
library(SeuratData)
library(tibble)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Rmagic)
library(reticulate)
use_condaenv('magic')
source("citeseq_function.R")
dataspat <- readRDS("data/spatial/HCC1T.rds")


# correlation ====================

dataspat <- dataspat[,dataspat$ident == "tumor"]
select = c("ITGA6", "CCL4", "CCR5")

plot <- spatial_exp(dataspat,
                    select = select,
                    assay = "Spatial",
                    outlier.rm = F)

name = "figure/spatial/HCC1T_correlation.pdf"
ggsave(name, plot, width = 37.5, height = 12.5)


# calculation ====================

dataspat <- dataspat[,dataspat$ident == "tumor"]

dataimpu <- magic(dataspat, genes = c("ITGA6", "CCL4"),
                  npca = 50, n.jobs = 30)

data <- GetAssayData(dataimpu, slot = "data", assay = "MAGIC_Spatial")
data <- data.frame(t(data))
cor.test(data$ITGA6, data$CCL4, method = "spearman")


# feature visualization ====================

select = c("Bcell", "CD4T", "CD8T", "NKcell", "Treg",
           "Macro", "Mono", "DC", "Neutrophil", 
           "Endothelial", "Fibroblast", "tumor")

plot <- spatial_feature(dataspat = dataspat,
                        select = select,
                        outlier.rm = F,
                        ncol = 3)

name = "figure/spatial/HCC1T_cellprop.pdf"
ggsave(name, plot, width = 37.5, height = 50, limitsize = FALSE)


# dimplot visualization ====================

info <- dataspat@meta.data
info$refined_pred <- as.factor(info$refined_pred)
dataspat@meta.data <- info

select = c("refined_pred")
plot <- spatial_dimplot(dataspat = dataspat,
                        select = select,
                        color = "discrete")

name = "figure/spatial/HCC1T_cluster.pdf"
ggsave(name, plot, width = 12.5, height = 12.5, limitsize = FALSE)


# dimplot neutrophil spot ====================

info <- dataspat@meta.data
info$neu_type <- "L"
info$neu_type[info$Neutrophil > quantile(info$Neutrophil, probs = c(0.9))] <- "H"
dataspat@meta.data <- info

select = c("neu_type")
plot <- spatial_dimplot(dataspat = dataspat,
                        select = select,
                        color = "continuous")

name = "figure/spatial/HCC1T_neutrophil.pdf"
ggsave(name, plot, width = 12.5, height = 12.5, limitsize = FALSE)



# neutrophil - tumor distance ====================

dataspat <- readRDS("data/spatial/HCC1T.rds")
info <- dataspat@meta.data
select1 <- rownames(info)[info$Neutrophil > 
                          quantile(info$Neutrophil, probs = c(0.9))]

select2 <- rownames(info)[info$ident == "tumor"]
dist <- spot_distance(dataspat = dataspat,
                      select1 = select1,
                      select2 = select2)

input <- do.call(rbind, lapply(colnames(dist), function(i){
  data.frame(id = i, dist = min(dist[,i]))
}))

dataimpu <- magic(dataspat, genes = c("ITGA6"),
                  npca = 50, n.jobs = 30, assay = "Spatial")

data <- GetAssayData(dataimpu, slot = "data",
                     assay = "MAGIC_Spatial")

data <- data.frame(t(data), check.names = F)
data <- rownames_to_column(data, "id")
input <- merge(input, data, by = "id")

name = "figure/spatial/distance_cor_HCC1T.txt"
write.table(input, name, sep = "\t", quote = F,
            row.names = F, col.names = T)

corplot(input = input[,-1],
        method = "spearman")

name = "figure/spatial/distance_cor_HCC1T.txt"
input <- read.table(name, sep = "\t", header = T, check.names = F)
input <- arrange(input, ITGA6)
input$number <- 1:nrow(input)

plot <- ggplot(input, aes(y=dist,x=number))+
  geom_smooth(span = 1,  method = "loess", se = TRUE) +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black", size = 15),
        axis.title = element_text(colour = "black", size = 15),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 1))

name = "figure/spatial/HCC1T_neutro_distance1.pdf"
ggsave(name, plot, width = 7, height = 3.2, limitsize = FALSE)


plot <- ggplot(input, aes(y=dist,x=number))+
  geom_line(aes(y=ITGA6,x=number), color="#E8BF80")+
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black", size = 15),
        axis.title = element_text(colour = "black", size = 15),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 1))

name = "figure/spatial/HCC1T_neutro_distance2.pdf"
ggsave(name, plot, width = 7, height = 3.2, limitsize = FALSE)
