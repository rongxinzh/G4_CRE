
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

hg19tohg38 <- function(data = NULL) {
  # liftOver oldFile map.chain newFile unMapped
  fwrite(data, "./tmp_hg19_ELS.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  system(paste0("liftOver ./tmp_hg19_ELS.txt ../data/ref/hg19ToHg38.over.chain.gz ./tmp_hg38_ELS.txt unMapped"))
  map.data <- fread("./tmp_hg38_ELS.txt", sep = "\t", header = FALSE) %>% data.frame()
  unlink("./tmp_hg19_ELS.txt")
  unlink("./tmp_hg38_ELS.txt")
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

loop <- fread("../data/TAD/hg19/41586_2020_2151_MOESM5_ESM.txt", sep = "\t", header = TRUE) %>% data.frame()
loop1 <- loop[, 1:3]
loop2 <- loop[, 4:6]
colnames(loop1) <- colnames(loop2) <- c("chr", "start", "end")
loop.anchor <- bind_rows(loop1, loop2)
loop.anchor <- hg19tohg38(loop.anchor)
loop.anchor.gr <- GrData(loop.anchor, center = FALSE)

ccre.pls.G4.gr <- GrData(ccre.pls.G4[, 1:3])
ccre.pls.nonG4.gr <- GrData(ccre.pls.nonG4[, 1:3])
ccre.pels.G4.gr <- GrData(ccre.pels.G4[, 1:3])
ccre.pels.nonG4.gr <- GrData(ccre.pels.nonG4[, 1:3])
ccre.dels.G4.gr <- GrData(ccre.dels.G4[, 1:3])
ccre.dels.nonG4.gr <- GrData(ccre.dels.nonG4[, 1:3])

extend.width <- 10000
width = 100
pls.G4.mat <- normalizeToMatrix(loop.anchor.gr, ccre.pls.G4.gr,
                                extend = extend.width, mean_mode = "coverage", w = width)
gc()

pls.nonG4.mat <- normalizeToMatrix(loop.anchor.gr, ccre.pls.nonG4.gr,
                                   extend = extend.width, mean_mode = "coverage", w = width)
gc()
pels.G4.mat <- normalizeToMatrix(loop.anchor.gr, ccre.pels.G4.gr,
                                 extend = extend.width, mean_mode = "coverage", w = width)
gc()

pels.nonG4.mat <- normalizeToMatrix(loop.anchor.gr, ccre.pels.nonG4.gr,
                                    extend = extend.width, mean_mode = "coverage", w = width)
gc()
dels.G4.mat <- normalizeToMatrix(loop.anchor.gr, ccre.dels.G4.gr,
                                 extend = extend.width, mean_mode = "coverage", w = width)
gc()

dels.nonG4.mat <- normalizeToMatrix(loop.anchor.gr, ccre.dels.nonG4.gr,
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
      theme_classic() + rremove("xlab") + ylab("Loop anchor density") + rremove("legend.title") +
      scale_color_manual(values = paletteer_d("ggthemes::excel_Main_Event"))
ggsave("../figure/loop/PLS_cCRE.pdf", p1, width = 4, height = 2.5)

p2 <- ggplot(mat2, aes(x = range, y = value, group = group)) +
      geom_line(aes(color = group)) +
      theme_classic() + rremove("xlab") + ylab("Loop anchor density") + rremove("legend.title") +
      scale_color_manual(values = paletteer_d("ggthemes::excel_Main_Event"))
ggsave("../figure/loop/pELS_cCRE.pdf", p2, width = 4, height = 2.5)

p3 <- ggplot(mat3, aes(x = range, y = value, group = group)) +
      geom_line(aes(color = group)) +
      theme_classic() + rremove("xlab") + ylab("Loop anchor density") + rremove("legend.title") +
      scale_color_manual(values = paletteer_d("ggthemes::excel_Main_Event"))
ggsave("../figure/loop/dELS_cCRE.pdf", p3, width = 4, height = 2.5)
