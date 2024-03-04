library(viper)
source("viper_function.R")
source("citeseq_function.R")


# HCC treatment =====================================================================

bpdata <- readRDS("data/bayesprism/bpdata_sorafenib.rds")
data <- t(get.exp(bpdata, state.or.type = "type", cell.name = "tumor"))
data <- log2(data)

info <- read.table("info.txt", sep = "\t",
                   header = T, check.names = F, row.names = 1)

info <- info[info$treatment == "Sor",]
info <- data.frame(id = rownames(info),
                   type = info$response)

input <- data.frame(id = colnames(data),
                    value = data["ITGA6",])

input <- na.omit(input)
input <- merge(input, info, by = "id")
input <- input[,-1]

plot <- common_dotbox(input, rotate = 45,
                      method = "wilcox.test", 
                      order = TRUE, decreasing = T)

name = "figure/treatment/sorafenib_puri.pdf"
ggsave(name, plot, width = 3.5, height = 6)


# HCC ICB1 =================================================================

bpdata <- readRDS("data/bayesprism/bpdata_2023_Nanjing.rds")
data <- t(get.exp(bpdata, state.or.type = "type", cell.name = "tumor"))
data <- log2(data + 1)

name = "ICB_HCC_cohorts_clinical.rds"
info <- readRDS(name)
info <- info[["2023_Nanjing"]]

info <- data.frame(id = info$id,
                   type = info$response1)

input <- data.frame(id = colnames(data),
                    value = data["PVR",])

input <- na.omit(input)
input <- merge(input, info, by = "id")
input <- input[,-1]

plot <- common_dotbox(input, rotate = 45,
                      method = "wilcox.test", 
                      order = TRUE, decreasing = T)

name = "figure/treatment/2023_Nanjing.pdf"
ggsave(name, plot, width = 3.5, height = 6)


# HCC ICB2 =================================================================

bpdata <- readRDS("data/bayesprism/bpdata_imbrave_150.rds")
data <- t(get.exp(bpdata, state.or.type = "type", cell.name = "tumor"))
data <- log2(data + 1)

name = "ICB_HCC_cohorts_clinical.rds"
info <- readRDS(name)
info <- info[["imbrave_150"]]

info <- data.frame(id = info$id,
                   type = info$response1)

input <- data.frame(id = colnames(data),
                    value = data["PVR",])

input <- na.omit(input)
input <- merge(input, info, by = "id")
input <- input[,-1]
input <- na.omit(input)

plot <- common_dotbox(input, rotate = 45,
                      method = "wilcox.test", 
                      order = TRUE, decreasing = T)

name = "figure/treatment/2022_imbrave_150.pdf"
ggsave(name, plot, width = 3.5, height = 6)


# survival data analysis ----------

bpdata <- readRDS("data/bayesprism/bpdata_imbrave_150.rds")
data <- t(get.exp(bpdata, state.or.type = "type", cell.name = "tumor"))
data <- log2(data + 1)

name = "ICB_HCC_cohorts_clinical.rds"
info <- readRDS(name)
info <- info[["imbrave_150"]]

info <- data.frame(id = info$id,
                   time = info$os_day,
                   status = info$os_censor)

info <- na.omit(info)

input <- data.frame(id = colnames(data),
                    value = data["PVR",])

input <- merge(input, info, by = "id")
input <- input[,c(2:4)]

survival_curve_gene(input = input,
                    group.by = "optimal",
                    censor = T,
                    confid = T,
                    name = "cohort")
