
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

tfbs <- fread("../data/tf/remap2022_nr_macs2_hg38_v1_0.bed.gz", sep = "\t", header = FALSE) %>% data.frame() %>% filter(V1 %in% paste0("chr", c(1:22, "X")))
tfbs.gr <- GrData(tfbs, center = FALSE)

extend.width <- 2500
width = 50
pls.G4.mat <- normalizeToMatrix(tfbs.gr, ccre.pls.G4.gr,
                                extend = extend.width, mean_mode = "coverage", w = width)
gc()

pls.nonG4.mat <- normalizeToMatrix(tfbs.gr, ccre.pls.nonG4.gr,
                                   extend = extend.width, mean_mode = "coverage", w = width)
gc()

pels.G4.mat <- normalizeToMatrix(tfbs.gr, ccre.pels.G4.gr,
                                 extend = extend.width, mean_mode = "coverage", w = width)
gc()

pels.nonG4.mat <- normalizeToMatrix(tfbs.gr, ccre.pels.nonG4.gr,
                                    extend = extend.width, mean_mode = "coverage", w = width)
gc()

dels.G4.mat <- normalizeToMatrix(tfbs.gr, ccre.dels.G4.gr,
                                 extend = extend.width, mean_mode = "coverage", w = width)
gc()

dels.nonG4.mat <- normalizeToMatrix(tfbs.gr, ccre.dels.nonG4.gr,
                                    extend = extend.width, mean_mode = "coverage", w = width)
gc()

mat1 <- bind_rows(data.frame(value = colMeans(pls.G4.mat), group = "G4 PLS", range = c(-(extend.width/width):(extend.width/width))),
                  data.frame(value = colMeans(pls.nonG4.mat), group = "non-G4 PLS", range = c(-(extend.width/width):(extend.width/width))))
mat2 <- bind_rows(data.frame(value = colMeans(pels.G4.mat), group = "G4 pELS", range = c(-(extend.width/width):(extend.width/width))),
                  data.frame(value = colMeans(pels.nonG4.mat), group = "non-G4 pELS", range = c(-(extend.width/width):(extend.width/width))))
mat3 <- bind_rows(data.frame(value = colMeans(dels.G4.mat), group = "G4 dELS", range = c(-(extend.width/width):(extend.width/width))),
                  data.frame(value = colMeans(dels.nonG4.mat), group = "non-G4 dELS", range = c(-(extend.width/width):(extend.width/width))))

p1 <- ggplot(mat1, aes(x = range, y = value, group = group)) +
      geom_line(aes(color = group)) +
      theme_classic() + rremove("xlab") + ylab("TF binding density") + rremove("legend.title") +
      scale_color_manual(values = paletteer_d("ggthemes::excel_Main_Event"))
ggsave("../figure/TF/remap_PLS_cCRE.pdf", p1, width = 4, height = 2.5)

p2 <- ggplot(mat2, aes(x = range, y = value, group = group)) +
      geom_line(aes(color = group)) +
      theme_classic() + rremove("xlab") + ylab("TF binding density") + rremove("legend.title") +
      scale_color_manual(values = paletteer_d("ggthemes::excel_Main_Event"))
ggsave("../figure/TF/remap_pELS_cCRE.pdf", p2, width = 4, height = 2.5)

p3 <- ggplot(mat3, aes(x = range, y = value, group = group)) +
      geom_line(aes(color = group)) +
      theme_classic() + rremove("xlab") + ylab("TF binding density") + rremove("legend.title") +
      scale_color_manual(values = paletteer_d("ggthemes::excel_Main_Event"))
ggsave("../figure/TF/remap_dELS_cCRE.pdf", p3, width = 4, height = 2.5)

# save(list=ls(), file = "./TFfootprint_G4cCRE.RData")
