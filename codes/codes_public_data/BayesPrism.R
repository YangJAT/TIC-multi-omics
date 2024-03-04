library(Seurat)
library(stringr)
library(tibble)
library(BayesPrism)
source("citeseq_function.R")


# reference ===========================================================

datafilt <- readRDS("/home/ug0302/CITEseq2/data/HCC_from_tumor_Nature.rds")
datafilt <- datafilt[,sample(colnames(datafilt), 50000)]

datafilt_prism <- bayesprism_reference(datafilt, cluster = "celltype_sig2")
saveRDS(datafilt_prism, "data/bayesprism/datafilt_prism_5.rds")


# bulk purification =================================================================

datafilt_prism <- readRDS("data/bayesprism/datafilt_prism_3.rds")


# sorafenib ====================

bkdata <- read.table("expdata.txt", sep = "\t",
                     header = T, check.names = F, row.names = 1)

bpdata <- bayesprism_analysis(scdata = datafilt_prism,
                              bkdata = bkdata,
                              celltype = "celltype_sig2",
                              cellstate = "state",
                              ncores = 50)

saveRDS(bpdata, "data/bayesprism/bpdata_sorafenib.rds")

prop <- get.fraction(bpdata,
                     which.theta = "final",
                     state.or.type = "type")

puridata <- get.exp(bpdata,
                    state.or.type = "type",
                    cell.name = "tumor")


# immunotherapy datasets ====================

file <- "ICB_HCC_cohorts_exp.rds"
HCC_cohorts <- readRDS(file)

bpICB <- list()
for (i in names(HCC_cohorts)) {
  
  cohort <- HCC_cohorts[[i]]
  cohort <- 2^cohort-1
  bpdata <- bayesprism_analysis(scdata = datafilt_prism,
                                bkdata = cohort,
                                celltype = "celltype_sig2",
                                cellstate = "state",
                                ncores = 50)
  
  name = paste0("/home/ug0302/CITEseq2/data/bayesprism", "/bpdata_", i, ".rds")
  saveRDS(bpdata, name)
}

