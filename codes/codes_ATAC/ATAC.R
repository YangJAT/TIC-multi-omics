library(ChIPQC)
library(DiffBind)
library(Rsubread)
library(ChIPseeker)
library(clusterProfiler)
library(rtracklayer)
library(GenomicAlignments)

source("function_bulkATAC.R")

genome = "hg38"

if (genome == "hg19") {
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
  bsgenome = BSgenome.Hsapiens.UCSC.hg19
  
} else if (genome == "hg38") {
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
  bsgenome = BSgenome.Hsapiens.UCSC.hg38
  
} else if (genome == "mm9") {
  library(BSgenome.Mmusculus.UCSC.mm9)
  library(TxDb.Mmusculus.UCSC.mm9.knownGene)
  txdb = TxDb.Mmusculus.UCSC.mm9.knownGene
  bsgenome = BSgenome.Mmusculus.UCSC.mm9
  
} else if (genome == "mm10") {
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
  bsgenome = BSgenome.Mmusculus.UCSC.mm10
}



# QC =====================================================================

input = "/home/ug0302/CITEseq3/06_ATAC2/bam"
output = "/home/ug0302/CITEseq3/06_ATAC2/figure/basic"

atac_qcplot(input = input,
            output = output,
            pattern = "bam$")

atac_tss(input = input,
         output = output,
         pattern = "bam$",
         genome = genome)


# DiffBind ====================

input = "/home/ug0302/CITEseq3/06_ATAC2/diffbind/diffbind_info.csv"

dbobj <- dba(sampleSheet = input)
dbobj <- dba.count(dbobj, bUseSummarizeOverlaps = TRUE)

output = "/home/ug0302/CITEseq3/06_ATAC2/figure/basic"

diffbind_qc(dbobj = dbobj,
            output = output)

output = "/home/ug0302/CITEseq3/06_ATAC2/peakd"

diffbind_bed(dbobj = dbobj,
             method = DBA_DESEQ2,
             FDR = FALSE, 
             thres.p = 0.05,
             thres.fc = 2,
             output = output)

diff <- diffbind_dataframe(dbobj = dbobj,
                           method = DBA_DESEQ2,
                           FDR = FALSE,
                           thres.p = 0.05,
                           thres.fc = 2)

diff <- data.frame(diff)
diff$Fold <- ifelse(diff$Fold > 0, "H", "L")
stat <- data.frame(table(diff$Fold))

plot <- circ_prop(stat)
ggsave("diffbind/prop.pdf", plot, width = 5, height = 5)


# Track plot ITGA6 ====================

path.bw = "/home/ug0302/CITEseq3/06_ATAC2/bigwig_new"
path.bed = "/home/ug0302/CITEseq3/06_ATAC2/diffbind"

name = "figure/trackplot/trackplot_ITGA6_H.pdf"
pdf(name, height = 6, width = 10)

trackplot(path.bw = path.bw,
          path.bed = path.bed,
          pattern.bw = "all_",
          pattern.bed = ".bed",
          region = "ITGA6",
          highlight = NULL,
          extend = 1000,
          ylim = c(0,9), 
          genome = "hg38")
dev.off()


# Track plot PVR ====================

path.bw = "/home/ug0302/CITEseq3/06_ATAC2/bigwig_new"
path.bed = "/home/ug0302/CITEseq3/06_ATAC2/diffbind"

name = "figure/trackplot/trackplot_PVR.pdf"
pdf(name, height = 6, width = 10)

trackplot(path.bw = path.bw,
          path.bed = path.bed,
          pattern.bw = "all_",
          pattern.bed = ".bed",
          region = "PVR",
          highlight = NULL,
          extend = 0,
          ylim = c(0,4), 
          genome = "hg38")
dev.off()


# Counting reads

rg <- GRanges(seqnames = "chr2",
              ranges = IRanges(start = c(172455280),
                               end = c(172455681)))

path = "/home/ug0302/CITEseq3/06_ATAC2/bam_new"
reads <- counting_reads(bam = path, 
                        pattern = ".bam$", 
                        region = rg,
                        scale.factor = 10000000,
                        genome = "hg38")

reads$type <- c(rep("H", 3), rep("L", 3))
reads <- data.frame(value = reads[,2], type = reads$type)

plot <- barplot_stat(input = reads,
                     rotate = 0,
                     order = TRUE,
                     show.point = TRUE,
                     output = NULL)

name = "figure/trackplot/stat_ITGA6.pdf"
ggsave(name, plot, width = 3.5, height = 5)


# Annotating Peaks ====================

anno <- annotatePeak(diff, TxDb = txdb)
plotAnnoPie(anno)
anno <- data.frame(anno)
anno <- anno[abs(anno$distanceToTSS) < 500,]



# Motif analysis ===============================================================

# One motif vs All sequence ====================

library(MotifDb)
library(JASPAR2022)
library(TFBSTools)
library(seqLogo)

# select motif

mtf <- getMatrixByName(JASPAR2022, "STAT3")

input <- as.matrix(mtf) / colSums(as.matrix(mtf))
seqLogo(input, ic.scale = TRUE)

mtfloc <- match_motif2chr(mtf = mtf,
                          min.score = "75%",
                          select.chr = "chr19",
                          genome = "hg19")

name = "/home/ug0302/CITEseq3/06_ATAC2/peakd/STAT3_motif.bed"
export(mtfloc, name, format = "bed")


# All motif vs One sequence ====================

library(MotifDb)
library(JASPAR2022)
library(TFBSTools)
library(seqLogo)
library(motifmatchr)

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

name = "/home/ug0302/CITEseq3/06_ATAC2/peakg/all_CD49H_peaks.narrowPeak"
input <- import(name, format = "narrowPeak")

select = c("all_CD49H_peak_29656", "all_CD49H_peak_34162")
input <- input[mcols(input)$name %in% select]

seq <- getSeq(bsgenome, input)
names(seq) <- as.character(input)

