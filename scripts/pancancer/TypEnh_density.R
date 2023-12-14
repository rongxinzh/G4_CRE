
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
suppressPackageStartupMessages(library(ggh4x))

GetEnh <- function(files = NULL, atac = NULL) {
  all.enh <- NULL
  for (tmp.f in files) {
    tmp.enh <- fread(paste0("../data/pancEnh/typical_enhancer/", tmp.f), sep = "\t", header = TRUE) %>% data.frame() %>% filter(Chr %in% paste0("chr", c(1:22, "X")))
    tmp.enh <- tmp.enh[, c(2:4, 6)]
    tmp.enh <- bt.intersect(a = tmp.enh, b = atac, wa = TRUE) %>% unique() %>% data.frame()
    colnames(tmp.enh) <- c("chr", "start", "end", "score")
    tmp.enh <- tmp.enh %>% filter(chr %in% paste0("chr", c(1:22, "X")))
    all.enh <- bind_rows(all.enh, tmp.enh)
  }
  return(all.enh)
}

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

# ct used
ct.used <- fread("../data/PANC_CT_used.txt", sep = "\t", header = FALSE) %>% unlist() %>% c()

# load ccre
ccre <- fread("../data/ccre/V3/G4_cCRE_annotation.txt", sep = "\t", header = TRUE) %>% data.frame() %>% filter(chr %in% paste0("chr", c(1:22, "X")))
ccre.enh <- ccre %>% filter(encodeLabel %in% c("pELS", "dELS"))

# load typical enhancers
all.enh.files <- list.files("../data/pancEnh/typical_enhancer/")

pan.atac.files <- list.files("../data/pancATAC/TCGA-ATAC_Cancer_Type-specific_PeakCalls/")

#tissue.use <- intersect(str_split(pan.atac.files, "_",simplify = TRUE)[, 1], unique(str_split(all.enh.files, "_", simplify = TRUE)[, 1]))
tissue.use <- ct.used

mat <- NULL
for (tmp.atac.file in pan.atac.files) {
  tmp.ct <- str_split(tmp.atac.file, "_", simplify = TRUE)[1]
  if (tmp.ct %in% tissue.use) {
    message(tmp.ct)

    tmp.atac <- fread(paste0("../data/pancATAC/TCGA-ATAC_Cancer_Type-specific_PeakCalls/", tmp.atac.file), sep = "\t", header = TRUE) %>% data.frame()
    tmp.ccre.enh <- bt.intersect(a = ccre.enh, b = tmp.atac, wa = TRUE) %>% unique() %>% data.frame()
    tmp.ccre.enh <- ccre.enh

    tmp.ccre.enh.p1 <- tmp.ccre.enh %>% filter(encodeLabel == "pELS", contain_G4 == 1) %>% GrData()
    tmp.ccre.enh.p0 <- tmp.ccre.enh %>% filter(encodeLabel == "pELS", contain_G4 == 0) %>% GrData()
    tmp.ccre.enh.d1 <- tmp.ccre.enh %>% filter(encodeLabel == "dELS", contain_G4 == 1) %>% GrData()
    tmp.ccre.enh.d0 <- tmp.ccre.enh %>% filter(encodeLabel == "dELS", contain_G4 == 0) %>% GrData()

    tmp.enh.files <- all.enh.files[str_starts(all.enh.files, paste0(tmp.ct, "_"))]
    tmp.enh <- GetEnh(tmp.enh.files, tmp.atac)
    tmp.enh.gr <- GrData(bt.sort(tmp.enh), center = FALSE, noscore = FALSE)

    extend.width <- 2500
    width = 10
    p1.mat <- normalizeToMatrix(tmp.enh.gr, tmp.ccre.enh.p1, value_column = "score",
                                extend = extend.width, mean_mode = "coverage", w = width)
    gc()
    p0.mat <- normalizeToMatrix(tmp.enh.gr, tmp.ccre.enh.p0, value_column = "score",
                                extend = extend.width, mean_mode = "coverage", w = width)
    gc()
    d1.mat <- normalizeToMatrix(tmp.enh.gr, tmp.ccre.enh.d1, value_column = "score",
                                extend = extend.width, mean_mode = "coverage", w = width)
    gc()
    d0.mat <- normalizeToMatrix(tmp.enh.gr, tmp.ccre.enh.d0, value_column = "score",
                                extend = extend.width, mean_mode = "coverage", w = width)
    gc()
    window <- extend.width/width
    mat <- bind_rows(mat, bind_rows(
                                    data.frame(value = colMeans(p1.mat), group = "G4 cCREs", cCRE = "pELS", range = c(-window:window), ct = tmp.ct),
                                    data.frame(value = colMeans(p0.mat), group = "Other cCREs", cCRE = "pELS", range = c(-window:window), ct = tmp.ct),
                                    data.frame(value = colMeans(d1.mat), group = "G4 cCREs", cCRE = "dELS", range = c(-window:window), ct = tmp.ct),
                                    data.frame(value = colMeans(d0.mat), group = "Other cCREs", cCRE = "dELS", range = c(-window:window), ct = tmp.ct)
                                    ))
  }
}

save(list=ls(), file = "./TypEnh_density.RData")

p.f <- ggplot(mat %>% filter(cCRE == "pELS"), aes(x = range, y = value, group = group)) +
       geom_line(aes(color = group)) +
       theme_classic() + rremove("xlab") + ylab("Typical enhancer density") + rremove("legend.title") +
       theme(legend.position="top") + 
       scale_color_manual(values = paletteer_d("ggthemes::excel_Main_Event")) +
       facet_wrap(vars(ct), scale = "free_y", ncol = 5) 
       #facet_grid2(cCRE ~ ct, independent = "y", scales = "free") 
ggsave("../figure/pancancer/typical_enhancer_density_pELS.pdf", p.f, width = 10, height = 4.5)

d.f <- ggplot(mat %>% filter(cCRE == "dELS"), aes(x = range, y = value, group = group)) +
       geom_line(aes(color = group)) +
       theme_classic() + rremove("xlab") + ylab("Typical enhancer density") + rremove("legend.title") +
       theme(legend.position="top") + 
       scale_color_manual(values = paletteer_d("ggthemes::excel_Main_Event")) +
       facet_wrap(vars(ct), scale = "free_y", ncol = 5) 
       #facet_grid2(cCRE ~ ct, independent = "y", scales = "free") 
ggsave("../figure/pancancer/typical_enhancer_density_dELS.pdf", d.f, width = 10, height = 4.5)




