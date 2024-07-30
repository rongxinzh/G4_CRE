
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
ctcf.ccre <- ccre[str_detect(ccre$ccre, "CTCF-bound"),]
ctcfonly.ccre <- ctcf.ccre[str_detect(ctcf.ccre$ccre, "CTCF-only,CTCF-bound"),]
ctcfonly.ccre.gr <- GrData(ctcfonly.ccre, center = TRUE)
ctcfother.ccre <- ctcf.ccre[!str_detect(ctcf.ccre$ccre, "CTCF-only,CTCF-bound"),]
ctcfother.ccre.gr <- GrData(ctcfother.ccre, center = TRUE)

G4 <- fread("../data/G4/cCRE_G4_annotation.txt", sep = "\t", header = TRUE) %>% data.frame()
G4.gr <- GrData(G4, center = FALSE)

extend.width <- 1000
width = 10

ctcfonly.mat <- normalizeToMatrix(G4.gr, ctcfonly.ccre.gr,
                                  extend = extend.width, mean_mode = "coverage", w = width)
ctcfonly.mat.mean <- colMeans(ctcfonly.mat)
rm(ctcfonly.mat)
gc()


ctcfother.mat <- normalizeToMatrix(G4.gr, ctcfother.ccre.gr,
                                   extend = extend.width, mean_mode = "coverage", w = width)
ctcfother.mat.mean <- colMeans(ctcfother.mat)
rm(ctcfother.mat)
gc()

mat <- bind_rows(data.frame(value = ctcfonly.mat.mean, group = "CTCF-only", range = c(-(extend.width/width):(extend.width/width))),
                 data.frame(value = ctcfother.mat.mean, group = "CTCF-hybrid", range = c(-(extend.width/width):(extend.width/width))))

p <- ggplot(mat, aes(x = range, y = value, group = group)) +
     geom_line(aes(color = group)) +
     theme_classic() + rremove("xlab") + ylab("G4 density") + rremove("legend.title") +
     scale_color_manual(values = paletteer_d("ggthemes::excel_Main_Event"))
ggsave("./G4_ctcf_hybrid_only.pdf", p, width = 5, height = 2.5)
