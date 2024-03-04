library(tibble)

source("citeseq_function.R")
source("viper_function.R")

# scRNA ====================================================================

name <- "/home/ug0302/Cancer_sc/05_Allmerge/data/datafilt_A_HI_2022_Nature.rds"
datafilt <- readRDS(name)

datafilt <- datafilt[,datafilt$disease == "HCC"]
datafilt <- datafilt[,datafilt$origin == "tumor"]

subfilt <- datafilt[,datafilt$celltype_sig2 == "tumor"]
data <- meanexpression(subfilt, slot = "data", group = "sample")
data <- as.matrix(data)

saveRDS(data, "/home/ug0302/CITEseq2/figure/STAT3/scRNA_exp.rds")
data <- readRDS("/home/ug0302/CITEseq2/figure/STAT3/scRNA_exp.rds")
data <- RankTransform(data)

path = "/home/ug0302/CITEseq/public_data/stem_sig.csv"
scores <- sig_scores(data, method = "zscore",
                     source = path, geneset = NULL,
                     min.sz = 2,
                     max.sz = 1000)

protein <- WeightedVIPER(data, regulon_list_sub, ncores = 10,
                         regulon.size = 5, ret.weights = FALSE)

input1 <- data.frame(id = colnames(protein), STAT3 = protein["STAT3",])
input2 <- data.frame(id = colnames(scores), signature = scores["Sensitizer",])

input <- merge(input1, input2, by = "id")
input <- input[,-1]

plot <- corplot_density(input = input,
                        method = "spearman",
                        density.bin = 10,
                        density.h = 5)

name = "figure/STAT3/sc_Sensitizer.pdf"
ggsave(name, plot, width = 5.2, height = 5)


# bulk data =================================================================

select <- c("data/bayesprism/bpdata_TCGA.rds",
            "data/bayesprism/bpdata_LIRI.rds",
            "data/bayesprism/bpdata_LICA.rds",
            "data/bayesprism/bpdata_CHCC.rds")

cordata <- do.call(rbind, lapply(select, function(i){
  
  bpdata <- readRDS(i)
  data <- t(get.exp(bpdata, state.or.type = "type", cell.name = "tumor"))
  data <- log2(data + 1)
  data <- RankTransform(data)
  
  path = "/home/ug0302/CITEseq/public_data/stem_sig.csv"
  scores <- sig_scores(data, method = "zscore",
                       source = path, geneset = NULL,
                       min.sz = 2,
                       max.sz = 1000)
  
  protein <- WeightedVIPER(data, regulon_list_sub, ncores = 10,
                           regulon.size = 5, ret.weights = FALSE)
  
  select_sig <- c("Resistor", "Sensitizer")
  subcor <- do.call(rbind, lapply(select_sig, function(j){
    
    input1 <- data.frame(id = colnames(protein), STAT3 = protein["STAT3",])
    input2 <- data.frame(id = colnames(scores), signature = scores[j,])
    
    input <- merge(input1, input2, by = "id")
    stat <- cor.test(input$STAT3, input$signature, method = c("spearman"))
    data.frame(id = i,
               sig = j,
               cor = stat$estimate,
               pvalue = stat$p.value)
  }))
}))

input <- acast(cordata, id ~ sig, value.var = "cor",
               fun.aggregate = mean, na.rm = TRUE)

input <- data.frame(input)
input <- rownames_to_column(input, var = "id")

plot <- ggpaired(input, cond1 = "Resistor", cond2 = "Sensitizer",
                 fill = "condition", palette = "jco", point.size = 3)+
  stat_compare_means(method = "t.test",paired = TRUE)

name = "figure/STAT3/boxplot_cor.pdf"
ggsave(name, plot, width = 3, height = 3.5)

