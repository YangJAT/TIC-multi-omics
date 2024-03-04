library(Signac)
library(Seurat)
library(BSgenome)
library(GenomeInfoDb)
library(biovizBase)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v75) # hg19
library(EnsDb.Hsapiens.v86) # hg38
library(clusterProfiler)
library(ggplot2)
library(patchwork)
library(magrittr)
library(readr)
library(Matrix)
library(tidyr)
library(dplyr)
library(stringr)
source("citeseq_function.R")


# genome ===================================================================

genome = "hg38"
species = "human"
bsgenome <- getBSgenome(genome)

if (genome == "hg19"){
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)} else {
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)}
seqlevelsStyle(annotations) <- 'UCSC'

if (species == "human") {
  genesize = 3.3e+09} else {
  genesize = 3.0e+09}

bsgenome <- getBSgenome(genome)
chr_length <- seqlengths(bsgenome)
chr_length <- data.frame(chr = names(chr_length), len = chr_length)

select <- str_detect(chr_length$chr, paste(c("fix", "alt"), collapse = '|'))
chr_length <- chr_length[!select,]

chr_length_v <- chr_length$len
names(chr_length_v) <- chr_length$chr

if (genome == "hg19") {
  blacklist = blacklist_hg19} else {
  blacklist = blacklist_hg38_unified}


# Clustering ===================================================================

# 1ST clustering ==============================

fragment = "/home/ug0302/CITEseq3/05_scATAC/data/HC006/HC006_fragments.tsv.gz"
fragment_obj <- CreateFragmentObject(fragment)

binmat <- GenomeBinMatrix(fragments = fragment_obj,
                          process_n = 50000, sep = c("-", "-"),
                          genome = chr_length_v, binsize = 2500)

path = "/home/ug0302/CITEseq3/05_scATAC/data/HC006/metadata_HC006.csv"
metadata <- read.csv(file = path, header = TRUE, row.names = 1)

atac <- CreateChromatinAssay(counts = binmat,
                             sep = c("-", "-"),
                             genome = genome, 
                             min.cells = 10, min.features = 200, 
                             annotation = annotations, fragments = fragment)

data <- CreateSeuratObject(counts = atac,
                           assay = "peaks",
                           meta.data = metadata)


# filtering --------------------

data$pct_reads_in_peaks <- data$peak_region_fragments / data$passed_filters * 100
data$blacklist_ratio <- data$blacklist_region_fragments / data$peak_region_fragments

data <- TSSEnrichment(data)
data <- NucleosomeSignal(data)

par(mfrow = c(2,3))
hist(data$peak_region_fragments, breaks = 100, prob = TRUE)
hist(data$pct_reads_in_peaks, breaks = 100, prob = TRUE)
hist(data$blacklist_ratio, breaks = 100, prob = TRUE)
hist(data$TSS.enrichment, breaks = 100, prob = TRUE)
hist(data$nucleosome_signal, breaks = 100, prob = TRUE)

filter_cell <- colnames(data)[data$peak_region_fragments > 3000 &
                                data$peak_region_fragments < 60000 &
                                data$pct_reads_in_peaks > 30 &
                                data$blacklist_ratio < 0.025 &
                                data$TSS.enrichment > 3 &
                                data$nucleosome_signal < 4]

datatac <- data[,filter_cell]


# 4 normalization --------------------

datatac <- BinarizeCounts(datatac)

datatac <- RunTFIDF(datatac)
datatac <- FindTopFeatures(datatac, min.cutoff = 'q0')
datatac <- RunSVD(object = datatac, assay = 'peaks',
                  reduction.key = 'LSI_', reduction.name = 'lsi')

datatac <- RunUMAP(object = datatac, reduction = 'lsi', dims = 2:30)
datatac <- FindNeighbors(object = datatac, reduction = 'lsi', dims = 2:30)
datatac <- FindClusters(object = datatac, verbose = FALSE, algorithm = 3)
DimPlot(object = datatac, label = TRUE) + NoLegend()


# Peak calling --------------------

type = "seurat_clusters"
macs2 = "/pub/anaconda3/bin/macs2"

# peak calling

peaks <- CallPeaks(object = datatac, group.by = type,
                   effective.genome.size = genesize,
                   macs2.path = macs2)

peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist, invert = TRUE)


# 2ST clustering ==============================

submat <- FeatureMatrix(fragments = fragment_obj,
                        cells = colnames(datatac),
                        sep = c("-", "-"),
                        features = peaks,
                        process_n = 50000)

atac <- CreateChromatinAssay(counts = submat,
                             sep = c("-", "-"),
                             genome = genome, 
                             min.cells = 10, min.features = 200, 
                             annotation = annotations, fragments = fragment)

datatac <- CreateSeuratObject(counts = atac,
                              assay = "peaks",
                              meta.data = metadata)


# 3 normalization --------------------

datatac <- BinarizeCounts(datatac)

datatac <- RunTFIDF(datatac)
datatac <- FindTopFeatures(datatac, min.cutoff = 'q0')
datatac <- RunSVD(object = datatac, assay = 'peaks',
                  reduction.key = 'LSI_', reduction.name = 'lsi')

datatac <- RunUMAP(object = datatac, reduction = 'lsi', dims = 2:30)

datatac <- FindNeighbors(object = datatac, reduction = 'lsi', dims = 2:30)
datatac <- FindClusters(object = datatac, verbose = FALSE,
                        resolution = 0.3, algorithm = 3)

DimPlot(object = datatac, label = TRUE) + NoLegend()
datatac <- RegionStats(datatac, genome = bsgenome)


# gene activity =================================================================

gene.activities <- GeneActivity(datatac)
datatac[['activity1']] <- CreateAssayObject(counts = gene.activities)

datatac <- NormalizeData(object = datatac, assay = 'activity1',
                         normalization.method = 'LogNormalize',
                         scale.factor = median(datatac$nCount_peaks))

DefaultAssay(datatac) <- 'activity1'
FeaturePlot(object = datatac,
            features = c('ITGA6', 'PVR'),
            pt.size = 0.1, max.cutoff = 'q95', ncol = 3)

save(datatac, file = "data/datatac_HC006.Rdata")


# ATAC-RNA integration =========================================================

# scRNA-seq

load("/home/ug0302/CITEseq/savedata/data1.Rdata")
datarna <- allexp[["HC006"]]
datapro <- allpro[["HC006"]]

datarna <- RunUMAP(datarna, dims = 1:20,
                   reduction = "pca", 
                   reduction.name = "umap", 
                   n.neighbors = 20, min.dist = 0.5, return.model = TRUE)

# reference mapping

anchors <- FindTransferAnchors(reference = datarna, query = datatac,
                               features = VariableFeatures(object = datarna),
                               reference.assay = "RNA", query.assay = "activity1",
                               reduction = 'cca')
# label prediction

data <- GetAssayData(datapro, slot = "data", assay = "RNA")
data <- as.matrix(data)

label <- TransferData(anchorset = anchors, refdata = data,
                      weight.reduction = datatac[['lsi']], dims = 2:30)

tranferdata <- data.frame(t(as.matrix(label@data)))
datatac <- AddMetaData(object = datatac, metadata = tranferdata)

# co-embedding

datatac <- MapQuery(anchorset = anchors,
                    reference = datarna, query = datatac, 
                    reference.reduction = "pca", reduction.model = "umap",
                    projectumap.args = list(reduction.name = "umap_mapping"))

data <- GetAssayData(datapro, slot = "data", assay = "RNA")
data <- data.frame(t(as.matrix(data)))
datarna <- AddMetaData(datarna, data)

plot1 <- featureplot_new(data = datarna, pt.size = 2,
                         reduction = "umap", color = "ybb",
                         features = "CD49f", outlier.rm = F, order = F)

plot2 <- featureplot_new(data = datatac, pt.size = 2,
                         reduction = "umap_mapping", color = "ybb",
                         features = "CD49f", outlier.rm = T, order = F)

plot <- plot1 | plot2
ggsave("figure/HC006_CD49f_prediction.pdf", plot, width = 11, height = 5.5)

refdata <- Embeddings(datarna, reduction = 'umap') %>% data.frame()
refdata$celltype <- "refer"; names(refdata)[1:2] <- c("a", "b")

quydata <- Embeddings(datatac, reduction = 'umap_mapping') %>% data.frame()
quydata$celltype <- "query"; names(quydata)[1:2] <- c("a", "b")

plotdata <- rbind(refdata, quydata)

ggplot(plotdata, aes(x = plotdata[,1],y = plotdata[,2])) + 
  geom_point(aes(color = celltype), size = 0.75) + 
  scale_color_manual(values=c("#0070b2", "grey90")) + 
  
  labs(x = 'UMAP1',y= 'UMAP2',title = '') + 
  guides(colour = guide_legend(override.aes = list(size = 3))) + 
  
  theme_bw() + 
  theme(strip.background = element_rect(colour = NA, fill = 'grey90'),
        strip.text.x = element_text(size = 15), 
        panel.grid = element_blank(),
        strip.placement = 'outside',
        axis.ticks = element_blank(),
        axis.text = element_blank())



# CD49f group ============================================================

load("data/datatac_HC006.Rdata")

info <- datatac@meta.data
info$group <- "M"
info$group[info$CD49f > quantile(info$CD49f, probs = c(0.9))] <- "H"
info$group[info$CD49f < quantile(info$CD49f, probs = c(0.1))] <- "L"

datatac <- AddMetaData(datatac, info)
datatac <- datatac[,datatac$group %in% c("H", "L")]

select  = c("ITGA6")
plot <- exp_distribution_ecdf(datatac,
                              assay = "activity1",
                              select = select,
                              group = "group",
                              method = "wilcox.test",
                              ncol = 1)

ggsave("figure/HC006_ecdf_ITGA6.pdf",
       plot, width = 5.5, height = 5)

diff <- seurat_diffall(datatac,
                       min.pct = 0.05,
                       thres.fc = 0.05,
                       group.by = "group",
                       assay = "activity1")

write.table(diff, "figure/diff.txt", sep = "\t",
            quote = F, row.names = F, col.names = T)



# diff peaks =================================================================

library(future.apply)
plan(multisession, workers = 4)
options(future.globals.maxSize = 3000*1024^2)

diffpeak <- FindAllMarkers(datatac, 
                           only.pos = TRUE, test.use = 'LR', 
                           assay = "peaks", min.pct = 0.05,
                           latent.vars = "peak_region_fragments")

plan(sequential)
diffpeak_flt <- diffpeak[diffpeak$p_val < 0.005, ]

expdata <- as.matrix(GetAssayData(datatac, slot = "data", assay = "peaks"))

clustinfo <- datatac@meta.data
clustinfo <- data.frame(id = rownames(clustinfo),
                        cluster = clustinfo$group)

diffdata <- data.frame(id = diffpeak_flt$gene,
                       cluster = diffpeak_flt$cluster,
                       logfc = diffpeak_flt$avg_log2FC,
                       pvalue = diffpeak_flt$p_val_adj)

library(ComplexHeatmap)
library(gtools)

col <- colorRampPalette(rev(c("#861b20","#A23E48","#CD6981",
                              "#FBD0C0","#F0F7F0","white")))(100)

input <- lapply(as.character(unique(clustinfo$cluster)), function(i){
  select <- clustinfo$id[clustinfo$cluster == i]
  data.frame(rowMeans(expdata[,select]))
})

input <- do.call(cbind, input)
colnames(input) <- as.character(unique(clustinfo$cluster))

input <- input[,mixedsort(colnames(input))]
select <- diffdata$id[diffdata$cluster == "H"]
subinput1 <- input[select,]; subinput1 <- subinput1[order(subinput1$H, decreasing = T),]

select <- diffdata$id[diffdata$cluster == "L"]
subinput2 <- input[select,]; subinput2 <- subinput2[order(subinput2$L, decreasing = T),]

input <- rbind(subinput1, subinput2)

pdf(file = "figure/heatmap_diff.pdf", width = 5, height = 7)
Heatmap(input, col = col, 
        show_row_names = label, 
        show_column_names = T, 
        column_names_rot = col_rot,
        cluster_rows = F, show_row_dend = F,
        cluster_columns = F, show_column_dend = F,
        row_names_gp = gpar(fontsize = font.size))
dev.off()


# CD49f-hi group --------------------

regions_h <- diffpeak_flt$gene[diffpeak_flt$cluster == "H"]
regions_h <- StringToGRanges(regions_h)

regionMatrix_h <- RegionMatrix(datatac, regions = regions_h,
                               assay = "peaks", group.by = "group",
                               idents = c("H", "L"), key = "peaks",
                               upstream = 10000, downstream = 10000)

plot <- RegionHeatmap(regionMatrix_h,
                      key = 'peaks',
                      assay = 'peaks',
                      idents = c("H", "L"),
                      normalize = TRUE,
                      upstream = 1000,
                      cols = NULL,
                      downstream = 1000,
                      max.cutoff = "q95",
                      min.counts = 1,
                      window = (50)/1,
                      order = TRUE,
                      nrow = NULL)

ggsave("figure/HC006_diff_heatmap1.pdf",
       plot, width = 5.5, height = 5)


# CD49f-lo group --------------------

regions_l <- diffpeak_flt$gene[diffpeak_flt$cluster == "L"]
regions_l <- StringToGRanges(regions_l)

regionMatrix_l <- RegionMatrix(datatac, regions = regions_l,
                               assay = "peaks", group.by = "group",
                               idents = c("H", "L"), key = "peaks",
                               upstream = 10000, downstream = 10000)

plot <- RegionHeatmap(regionMatrix_l,
                      key = 'peaks',
                      assay = 'peaks',
                      idents = c("H", "L"),
                      normalize = TRUE,
                      upstream = 1000,
                      cols = NULL,
                      downstream = 1000,
                      max.cutoff = "q95",
                      min.counts = 1,
                      window = (50)/1,
                      order = TRUE,
                      nrow = NULL)

ggsave("figure/HC006_diff_heatmap2.pdf",
       plot, width = 5.5, height = 5)



# Motif analysis ===============================================================

library(JASPAR2022)
library(TFBSTools)
library(motifmatchr)
library(ggseqlogo)
library(chromVAR)
library(BSgenome)

DefaultAssay(datatac) <- 'peaks'


# Adding motif information --------------------

pfm <- getMatrixSet(x = JASPAR2022,
                    opts = list(species = "Homo sapiens", all_versions = FALSE))

datatac <- AddMotifs(object = datatac, pfm = pfm, 
                     genome = bsgenome, assay = "peaks")

opts <- list()
opts[["tax_group"]] <- "vertebrates"
opts[["collection"]] <- "CORE"
opts[["all_versions"]] <- FALSE
motifsToScan <- getMatrixSet(JASPAR2022, opts)

anno_motif <- do.call(rbind, lapply(names(motifsToScan), function(i){
  mtf <- motifsToScan@listData[[i]]
  data.frame(id = i, tf = mtf@name,
             species = unname(mtf@tags[["species"]]))
}))

anno_motif <- anno_motif[anno_motif$species %in% c("Homo sapiens"),]


# motif enrichment analysis --------------------

diffpeak_flt <- diffpeak[diffpeak$p_val < 0.005,]
select <- as.character(unique(diffpeak_flt$cluster))

diffmotif <- do.call(rbind, lapply(select, function(i){
  subdiff <- diffpeak_flt[diffpeak_flt$cluster == i,]$gene
  
  # background peaks
  
  open.peaks <- AccessiblePeaks(datatac, assay = "peaks", idents = c(i))
  meta.feature <- GetAssayData(datatac, assay = "peaks", slot = "meta.features")
  
  peaks.matched <- MatchRegionStats(
    meta.feature = meta.feature[open.peaks, ],
    query.feature = meta.feature[subdiff, ], n = 50000)
  
  motif <- FindMotifs(object = datatac,
                      features = subdiff,
                      background = peaks.matched)
  
  motif <- data.frame(gene = motif$motif.name, motif = motif$motif,
                      pvalue = motif$pvalue, adjustp = motif$p.adjust,
                      fold_enrich = motif$fold.enrichment, cluster = i)
  
  motif <- motif[motif$adjustp < 0.01,]
}))

write.table(diffmotif, "figure/diffmotif.txt", sep = "\t",
            quote = F, row.names = F, col.names = T)

diffmotif_h <- diffmotif[diffmotif$cluster == "H",]
diffmotif_h <- diffmotif_h[diffmotif_h$gene %in% c("KLF4", "MYC", "MYCN", "STAT3"),]

diffmotif_l <- diffmotif[diffmotif$cluster == "L",]
diffmotif_l <- diffmotif_l[diffmotif_l$gene %in% c("HNF4A", "CEBPA", "TBX3"),]

plot <- MotifPlot(object = datatac, motifs = c("MA1566.2"))
ggsave("figure/motif_TBX3.pdf", plot, width = 5, height = 3)



# coverage plot ==========================================================

DefaultAssay(datatac) <- 'peaks'

datatac <- SetIdent(datatac, value = "group")
levels(datatac) <- c("H", "L")

col_setting <- c("H" = "#0070B2",
                 "L" = "#5CB3DA")


# BMI1 --------------------

region = "BMI1"
# region = "chr2-172400000-172500000"

cov_plot <- CoveragePlot(object = datatac,
                         region = region, links = FALSE,
                         tile.size = 200,
                         tile.cells = 200,
                         window = 200,
                         extend.upstream = 2000,
                         extend.downstream = 0,
                         peaks = TRUE, annotation = TRUE) + 
            scale_fill_manual(values = col_setting)

ggsave("figure/HC006_trackplot_BMI1.pdf", cov_plot, width = 5, height = 3)


# ITGA6 --------------------

region = "ITGA6"
cov_plot <- CoveragePlot(object = datatac,
                         region = region, links = FALSE,
                         tile.size = 200,
                         tile.cells = 200,
                         window = 500,
                         extend.upstream = 5000,
                         extend.downstream = 5000,
                         peaks = TRUE, annotation = TRUE) + 
  scale_fill_manual(values = col_setting)

ggsave("figure/HC006_trackplot_ITGA6.pdf", cov_plot, width = 5, height = 3)


# FSCN1 --------------------

region = "FSCN1"
cov_plot <- CoveragePlot(object = datatac,
                         region = region, links = FALSE,
                         tile.size = 200,
                         tile.cells = 200,
                         window = 300,
                         extend.upstream = 3000,
                         extend.downstream = 0,
                         peaks = TRUE, annotation = TRUE) + 
  scale_fill_manual(values = col_setting)

ggsave("figure/HC006_trackplot_FSCN1.pdf", cov_plot, width = 5, height = 3)

