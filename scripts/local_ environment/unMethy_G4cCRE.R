
set.seed(1)
options(scipen=200)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(bedtoolsr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(EnrichedHeatmap))
suppressPackageStartupMessages(library(paletteer))
suppressPackageStartupMessages(library(circlize))

hg19tohg38 <- function(data = NULL) {
  # liftOver oldFile map.chain newFile unMapped
  fwrite(data, "./tmp_hg19_Methy.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  system(paste0("liftOver ./tmp_hg19_Methy.txt ../data/ref/hg19ToHg38.over.chain.gz ./tmp_hg38_Methy.txt unMapped"))
  map.data <- fread("./tmp_hg38_Methy.txt", sep = "\t", header = FALSE) %>% data.frame()
  unlink("./tmp_hg19_Methy.txt")
  unlink("./tmp_hg38_Methy.txt")
  unlink("unMapped")
  return(map.data)
}

GrData <- function(data = NULL, center = TRUE) {
  data <- data[, 1:3]
  if (center == TRUE) {
    data[, 2] <- ceiling((data[, 2] + data[, 3]) /  2) - 1
    data[, 3] <- data[, 2] + 1
  }
  data$score <- 1
  colnames(data)[1:3] <- c("chr", "start", "end")
  data <- makeGRangesFromDataFrame(data,
                                   keep.extra.columns = TRUE,
                                   ignore.strand = TRUE)

  return(data)
}

extend.width <- 2500
width = 10

# load ccre
ccre <- fread("../data/ccre/V3/G4_cCRE_annotation.txt", sep = "\t", header = TRUE) %>% data.frame() %>% filter(chr %in% paste0("chr", c(1:22, "X")))
ccre.pls <- ccre %>% filter(encodeLabel %in% c("PLS"))
ccre.pls.G4 <- ccre.pls %>% filter(contain_G4 > 0)
ccre.pls.nonG4 <- ccre.pls %>% filter(contain_G4 == 0)
ccre.pels <- ccre %>% filter(encodeLabel %in% c("pELS"))
ccre.pels.G4 <- ccre.pels %>% filter(contain_G4 > 0)
ccre.pels.nonG4 <- ccre.pels %>% filter(contain_G4 == 0)
ccre.dels <- ccre %>% filter(encodeLabel %in% c("dELS"))
ccre.dels.G4 <- ccre.dels %>% filter(contain_G4 > 0)
ccre.dels.nonG4 <- ccre.dels %>% filter(contain_G4 == 0)

ccre.pls.G4.gr <- GrData(ccre.pls.G4[, 1:3])
ccre.pls.nonG4.gr <- GrData(ccre.pls.nonG4[, 1:3])
ccre.pels.G4.gr <- GrData(ccre.pels.G4[, 1:3])
ccre.pels.nonG4.gr <- GrData(ccre.pels.nonG4[, 1:3])
ccre.dels.G4.gr <- GrData(ccre.dels.G4[, 1:3])
ccre.dels.nonG4.gr <- GrData(ccre.dels.nonG4[, 1:3])

mat1 <- NULL
mat2 <- NULL
mat3 <- NULL

blocks <- fread("../data/methy/GSE186458_blocks.s207.hg38.bed", sep = "\t", header = FALSE) %>% data.frame() %>% filter(V1 %in% paste0("chr", c(1:22, "X")))
blocks <- blocks %>% filter((V5 - V4) >= 4)
blocks.gr <- GrData(blocks, center = FALSE)

block.pls.G4.mat <- normalizeToMatrix(blocks.gr, ccre.pls.G4.gr,
                                      extend = extend.width, mean_mode = "coverage", w = width)
gc()

block.pls.nonG4.mat <- normalizeToMatrix(blocks.gr, ccre.pls.nonG4.gr,
                                         extend = extend.width, mean_mode = "coverage", w = width)
gc()

block.pels.G4.mat <- normalizeToMatrix(blocks.gr, ccre.pels.G4.gr,
                                       extend = extend.width, mean_mode = "coverage", w = width)
gc()

block.pels.nonG4.mat <- normalizeToMatrix(blocks.gr, ccre.pels.nonG4.gr,
                                          extend = extend.width, mean_mode = "coverage", w = width)
gc()

block.dels.G4.mat <- normalizeToMatrix(blocks.gr, ccre.dels.G4.gr,
                                       extend = extend.width, mean_mode = "coverage", w = width)
gc()

block.dels.nonG4.mat <- normalizeToMatrix(blocks.gr, ccre.dels.nonG4.gr,
                                          extend = extend.width, mean_mode = "coverage", w = width)
gc()

v1 <- colMeans(block.pls.G4.mat)
v2 <- colMeans(block.pls.nonG4.mat)
v3 <- colMeans(block.pels.G4.mat)
v4 <- colMeans(block.pels.nonG4.mat)
v5 <- colMeans(block.dels.G4.mat)
v6 <- colMeans(block.dels.nonG4.mat)

methy.files <- list.files("../data/methy/41586_2022_5580_MOESM5_ESM/") 
for (tmp.methy.f in methy.files) {
  message(tmp.methy.f)
  tmp.methy <- fread(paste0("../data/methy/41586_2022_5580_MOESM5_ESM/", tmp.methy.f), sep = "\t", header = TRUE) %>% data.frame()
  tmp.methy <- tmp.methy %>% filter(chr %in% paste0("chr", c(1:22, "X")))
  tmp.methy <- hg19tohg38(tmp.methy[, 1:3])
  tmp.methy.gr <- GrData(tmp.methy, center = FALSE)

  tmp.pls.G4.mat <- normalizeToMatrix(tmp.methy.gr, ccre.pls.G4.gr,
                                      extend = extend.width, mean_mode = "coverage", w = width)
  gc()

  tmp.pls.nonG4.mat <- normalizeToMatrix(tmp.methy.gr, ccre.pls.nonG4.gr,
                                         extend = extend.width, mean_mode = "coverage", w = width)
  gc()

  tmp.pels.G4.mat <- normalizeToMatrix(tmp.methy.gr, ccre.pels.G4.gr,
                                       extend = extend.width, mean_mode = "coverage", w = width)
  gc()

  tmp.pels.nonG4.mat <- normalizeToMatrix(tmp.methy.gr, ccre.pels.nonG4.gr,
                                          extend = extend.width, mean_mode = "coverage", w = width)
  gc()

  tmp.dels.G4.mat <- normalizeToMatrix(tmp.methy.gr, ccre.dels.G4.gr,
                                       extend = extend.width, mean_mode = "coverage", w = width)
  gc()

  tmp.dels.nonG4.mat <- normalizeToMatrix(tmp.methy.gr, ccre.dels.nonG4.gr,
                                          extend = extend.width, mean_mode = "coverage", w = width)
  gc()

  tmp.mat1 <- bind_rows(data.frame(value = colMeans(tmp.pls.G4.mat) / v1, group = "G4 PLS", range = c(-(extend.width/width):(extend.width/width)), celltp = str_split(tmp.methy.f, "\\.", simplify = TRUE)[, 1]),
                        data.frame(value = colMeans(tmp.pls.nonG4.mat) / v2, group = "non-G4 PLS", range = c(-(extend.width/width):(extend.width/width)), celltp = str_split(tmp.methy.f, "\\.", simplify = TRUE)[, 1]))
  tmp.mat2 <- bind_rows(data.frame(value = colMeans(tmp.pels.G4.mat) / v3, group = "G4 pELS", range = c(-(extend.width/width):(extend.width/width)), celltp = str_split(tmp.methy.f, "\\.", simplify = TRUE)[, 1]),
                        data.frame(value = colMeans(tmp.pels.nonG4.mat) / v4, group = "non-G4 pELS", range = c(-(extend.width/width):(extend.width/width)), celltp = str_split(tmp.methy.f, "\\.", simplify = TRUE)[, 1]))
  tmp.mat3 <- bind_rows(data.frame(value = colMeans(tmp.dels.G4.mat) / v5, group = "G4 dELS", range = c(-(extend.width/width):(extend.width/width)), celltp = str_split(tmp.methy.f, "\\.", simplify = TRUE)[, 1]),
                        data.frame(value = colMeans(tmp.dels.nonG4.mat) / v6, group = "non-G4 dELS", range = c(-(extend.width/width):(extend.width/width)), celltp = str_split(tmp.methy.f, "\\.", simplify = TRUE)[, 1]))

  mat1 <- bind_rows(mat1, tmp.mat1)
  mat2 <- bind_rows(mat2, tmp.mat2)
  mat3 <- bind_rows(mat3, tmp.mat3)

}

save(list=ls(), file = "./unMethy_G4cCRE.RData")

load("./unMethy_G4cCRE.RData")
mat1.G4 <- spread(mat1 %>% filter(group == "G4 PLS") %>% select(-group) %>% data.frame(), key = "celltp", value = "value")
rownames(mat1.G4) <- mat1.G4$range
mat1.G4 <- mat1.G4[, -1]
mat1.nonG4 <- spread(mat1 %>% filter(group == "non-G4 PLS") %>% select(-group) %>% data.frame(), key = "celltp", value = "value")
rownames(mat1.nonG4) <- mat1.nonG4$range
mat1.nonG4 <- mat1.nonG4[, -1]

mat2.G4 <- spread(mat2 %>% filter(group == "G4 pELS") %>% select(-group) %>% data.frame(), key = "celltp", value = "value")
rownames(mat2.G4) <- mat2.G4$range
mat2.G4 <- mat2.G4[, -1]
mat2.nonG4 <- spread(mat2 %>% filter(group == "non-G4 pELS") %>% select(-group) %>% data.frame(), key = "celltp", value = "value")
rownames(mat2.nonG4) <- mat2.nonG4$range
mat2.nonG4 <- mat2.nonG4[, -1]

mat3.G4 <- spread(mat3 %>% filter(group == "G4 dELS") %>% select(-group) %>% data.frame(), key = "celltp", value = "value")
rownames(mat3.G4) <- mat3.G4$range
mat3.G4 <- mat3.G4[, -1]
mat3.nonG4 <- spread(mat3 %>% filter(group == "non-G4 dELS") %>% select(-group) %>% data.frame(), key = "celltp", value = "value")
rownames(mat3.nonG4) <- mat3.nonG4$range
mat3.nonG4 <- mat3.nonG4[, -1]

col_fun1 = colorRamp2(c(0, 0.5, 1), c("#2683C6", "white", "#B80E0F"))
pdf("../figure/methy/PLS_unmethy.pdf", width = 6, height = 6)
Heatmap(t(mat1.G4), cluster_rows = TRUE, cluster_columns = FALSE, show_column_names = FALSE, name = "Unmethylated intensity", col = col_fun1) +
Heatmap(t(mat1.nonG4), cluster_rows = TRUE, cluster_columns = FALSE, show_column_names = FALSE, show_heatmap_legend = FALSE, col = col_fun1)
dev.off()

col_fun2 = colorRamp2(c(0, 0.3, 0.6), c("#2683C6", "white", "#B80E0F"))
pdf("../figure/methy/pELS_unmethy.pdf", width = 6, height = 6)
Heatmap(t(mat2.G4), cluster_rows = TRUE, cluster_columns = FALSE, show_column_names = FALSE, name = "Unmethylated intensity", col = col_fun2) +
Heatmap(t(mat2.nonG4), cluster_rows = TRUE, cluster_columns = FALSE, show_column_names = FALSE, show_heatmap_legend = FALSE, col = col_fun2)
dev.off()

col_fun3 = colorRamp2(c(0, 0.075, 0.15), c("#2683C6", "white", "#B80E0F"))
pdf("../figure/methy/dELS_unmethy.pdf", width = 6, height = 6)
Heatmap(t(mat3.G4), cluster_rows = TRUE, cluster_columns = FALSE, show_column_names = FALSE, name = "Unmethylated intensity", col = col_fun3) +
Heatmap(t(mat3.nonG4), cluster_rows = TRUE, cluster_columns = FALSE, show_column_names = FALSE, show_heatmap_legend = FALSE, col = col_fun3)
dev.off()

print("done!methy")