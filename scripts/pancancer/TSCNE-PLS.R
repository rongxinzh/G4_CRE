
set.seed(1)
options(scipen=200)

message("*************************************************")
message("********************TSCNE-ELS********************")
message("*************************************************")

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(bedtoolsr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(EnrichedHeatmap))
suppressPackageStartupMessages(library(paletteer))

GrData <- function(data = NULL, center = TRUE, noscore = TRUE) {

  if (noscore == TRUE) {
    data <- data[, 1:3]
    data$score <- 1
  } else {
    data <- data[, 1:4]
    colnames(data)[4] <- "score"
  }

  if (center == TRUE) {
    data[, 2] <- ceiling((data[, 2] + data[, 3]) /  2) - 1
    data[, 3] <- data[, 2] + 1
  }
  colnames(data)[1:3] <- c("chr", "start", "end")
  data <- makeGRangesFromDataFrame(data,
                                   keep.extra.columns = TRUE,
                                   ignore.strand = TRUE)

  return(data)
}

# load ccre
ccre <- fread("../data/ccre/V3/G4_cCRE_annotation.txt", sep = "\t", header = TRUE) %>% data.frame() %>% filter(chr %in% paste0("chr", c(1:22, "X")))
ccre.prom.p1 <- ccre %>% filter(encodeLabel %in% "PLS", contain_G4 == 1) %>% GrData()
ccre.prom.p0 <- ccre %>% filter(encodeLabel %in% "PLS", contain_G4 == 0) %>% GrData()

tscre.map <- fread("../data/TSCRE/H3K4me3/TSCRE_ID_H3K4me3.txt", sep = "\t", header = FALSE) %>% data.frame()
all.files <- list.files("../data/TSCRE/H3K4me3/", recursive = TRUE, pattern = "*-Annot.txt$", full.names = TRUE)
all.files <- data.frame(path = all.files, ID = str_split(all.files, "(/|CRE_)", simplify = TRUE)[, 7])
tscre.map2 <- fread("../data/TSCRE/H3K27ac/TSCRE_ID_H3K27ac.txt", sep = "\t", header = FALSE) %>% data.frame()
all.ct <- intersect(unique(tscre.map[, 2]), tscre.map2[, 2])

extend.width <- 2500
width = 10

pls.mat <- NULL
for (tmp.ct in all.ct) {
  message(tmp.ct)
  tmp.map <- tscre.map %>% filter(V2 == tmp.ct)
  tmp.f <- all.files %>% filter(ID %in% tmp.map[, 3])
  tmp.tscre <- NULL
  for (i in 1:nrow(tmp.f)) {
    tmp.tscre <- bind_rows(tmp.tscre, fread(tmp.f[i, 1], sep = "\t", header = TRUE) %>% data.frame())
  }
  tmp.tscre <- tmp.tscre %>% filter(Type == "Gained_In_Tumor")
  tmp.tscre <- tmp.tscre %>% filter(Seqnames %in% paste0("chr", c(1:22, "X")))
  tmp.tscre <- GrData(bt.merge(bt.sort(tmp.tscre[, c(3:5)])), center = FALSE, noscore = TRUE)

  p1.mat <- normalizeToMatrix(tmp.tscre, ccre.prom.p1, value_column = "score",
                              extend = extend.width, mean_mode = "coverage", w = width)
  gc()
  p0.mat <- normalizeToMatrix(tmp.tscre, ccre.prom.p0, value_column = "score",
                              extend = extend.width, mean_mode = "coverage", w = width)
  gc()
  window <- extend.width/width
  pls.mat <- bind_rows(pls.mat, bind_rows(
                                          data.frame(value = colMeans(p1.mat), group = "G4 cCREs", cCRE = "PLS", range = c(-window:window), ct = tmp.ct),
                                          data.frame(value = colMeans(p0.mat), group = "non-G4-associated cCREs", cCRE = "PLS", range = c(-window:window), ct = tmp.ct),
                                          ))

}

save(list=ls(), file = "./TSCNE-PLS.RData")
