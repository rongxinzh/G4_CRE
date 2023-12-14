
set.seed(1)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(bedtoolsr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(EnrichedHeatmap))
suppressPackageStartupMessages(library(paletteer))
suppressPackageStartupMessages(library(colorRamp2))

hg19tohg38 <- function(data = NULL) {
  # liftOver oldFile map.chain newFile unMapped
  fwrite(data, "./tmp_hg19_ELS.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  system(paste0("liftOver ./tmp_hg19_ELS.txt ./hg19ToHg38.over.chain.gz ./tmp_hg38_ELS.txt unMapped"))
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
ccre <- fread("./G4_cCRE_annotation.txt", sep = "\t", header = TRUE) %>% data.frame() %>% filter(chr %in% paste0("chr", c(1:22, "X")))
ctcf.ccre <- ccre[str_detect(ccre$ccre, "CTCF-bound"),]
ctcfonly.ccre <- ctcf.ccre[str_detect(ctcf.ccre$ccre, "CTCF-only,CTCF-bound"),]
ctcfonly.ccre.gr <- GrData(ctcfonly.ccre, center = TRUE)
ctcfother.ccre <- ctcf.ccre[!str_detect(ctcf.ccre$ccre, "CTCF-only,CTCF-bound"),]
ctcfother.ccre.gr <- GrData(ctcfother.ccre, center = TRUE)

### chromatin loops
loop <- fread("./41586_2020_2151_MOESM5_ESM.txt", sep = "\t", header = TRUE) %>% data.frame()
loop1 <- loop[, 1:3]
loop2 <- loop[, 4:6]
colnames(loop1) <- colnames(loop2) <- c("chr", "start", "end")
loop.anchor <- bind_rows(loop1, loop2)
loop.anchor.gr <- GrData(hg19tohg38(loop.anchor), center = FALSE)

extend.width <- 10000
width = 10
ctcfonly.mat <- normalizeToMatrix(loop.anchor.gr, ctcfonly.ccre.gr,
                                	extend = extend.width, mean_mode = "coverage", w = width)
ctcfonly.mat.mean <- colMeans(ctcfonly.mat)
rm(ctcfonly.mat)
gc()
ctcfother.mat <- normalizeToMatrix(loop.anchor.gr, ctcfother.ccre.gr,
                                	 extend = extend.width, mean_mode = "coverage", w = width)
ctcfother.mat.mean <- colMeans(ctcfother.mat)
rm(ctcfother.mat)
gc()

mat <- bind_rows(data.frame(value = ctcfonly.mat.mean, group = "CTCF-only", range = c(-(extend.width/width):(extend.width/width))),
                 data.frame(value = ctcfother.mat.mean, group = "CTCF-hybrid", range = c(-(extend.width/width):(extend.width/width))))

p <- ggplot(mat, aes(x = range, y = value, group = group)) +
     geom_line(aes(color = group)) +
     theme_classic() + rremove("xlab") + ylab("Loop anchor density") + rremove("legend.title") +
     scale_color_manual(values = paletteer_d("ggthemes::excel_Main_Event"))
ggsave("./chromatinloop_G4_ctcf.pdf", p, width = 4, height = 2.5)

