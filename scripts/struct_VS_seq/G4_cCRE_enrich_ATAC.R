
set.seed(1)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(bedtoolsr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(colorRamp2))
suppressPackageStartupMessages(library(paletteer))
suppressPackageStartupMessages(library(random))
suppressPackageStartupMessages(library(EnrichedHeatmap))

GetBW <- function(elements = NULL, bw.file = NULL) {

  elements <- elements[, 1:4]
  elements[, 2] <- ceiling((elements[, 2] + elements[, 3]) / 2)
  elements[, 2] <- elements[, 2] - 100
  elements[, 3] <- elements[, 2] + 200 
  #
  elements[, 2] <- elements[, 2] - 1

  rand.f <- randomStrings(n = 1, len = 20) %>% as.character
  fwrite(elements, rand.f, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  system(paste0("bigWigAverageOverBed ", bw.file," ", rand.f, " ", rand.f, ".out"))
  bw.file.score <- fread(paste0(rand.f, ".out"), sep = "\t", header = FALSE) %>% data.frame
  colnames(bw.file.score)[1] <- colnames(elements)[4] <- "id"
  bw.file.score <- elements %>% left_join(bw.file.score, by = "id")

  unlink(rand.f)
  unlink(paste0(rand.f, ".out"))

  bw.file.score[, 2] <- bw.file.score[, 2] + 1
  bw.file.score <- bw.file.score[, c(1:4, 9)]
  colnames(bw.file.score) <- c("chr", "start", "end", "id", "score")
  return(bw.file.score)
}


G4 <- fread("../data/G4/G4Hunter_w25_s1.5_hg38.txt", sep = "\t", header = FALSE) %>% data.frame()
G4[, 4] <- paste0("G4_", 1:dim(G4)[1])

k562.G4 <- fread("../data/G4/K562_G4_hg38.bed", sep = "\t", header = FALSE) %>% 
	data.frame() %>% filter(V1 %in% paste0("chr", c(1:22, "X")))
hepg2.G4 <- fread("../data/G4/HepG2_G4_hg38.bed", sep = "\t", header = FALSE) %>% 
	data.frame() %>% filter(V1 %in% paste0("chr", c(1:22, "X")))

k562.ccre <- fread("../data/ccre/V3/ENCFF414OGC_ENCFF806YEZ_ENCFF849TDM_ENCFF736UDR.7group.bed", sep = "\t", header = FALSE) %>% 
	data.frame() %>% filter(V1 %in% paste0("chr", c(1:22, "X")))
hepg2.ccre <- fread("../data/ccre/V3/ENCFF546MZK_ENCFF732PJK_ENCFF795ONN_ENCFF357NFO.7group.bed", sep = "\t", header = FALSE) %>% 
	data.frame() %>% filter(V1 %in% paste0("chr", c(1:22, "X")))

k562.pG4 <- bt.intersect(a = G4, b = k562.G4, wa = TRUE) %>% unique()
nonk562.pG4 <- bt.intersect(a = G4, b = k562.G4, wa = TRUE, v = TRUE) %>% unique()
hepg2.pG4 <- bt.intersect(a = G4, b = hepg2.G4, wa = TRUE) %>% unique()
nonhepg2.pG4 <- bt.intersect(a = G4, b = hepg2.G4, wa = TRUE, v = TRUE) %>% unique()

k562.pG4.pls <- bt.intersect(a = k562.pG4, b = k562.ccre %>% filter(V10 %in% c("PLS", "PLS,CTCF-bound")), wa = TRUE) %>% unique()
nonk562.pG4.pls <- bt.intersect(a = nonk562.pG4, b = k562.ccre %>% filter(V10 %in% c("PLS", "PLS,CTCF-bound")), wa = TRUE) %>% unique()

k562.pG4.pels <- bt.intersect(a = k562.pG4, b = k562.ccre %>% filter(V10 %in% c("pELS", "pELS,CTCF-bound")), wa = TRUE) %>% unique()
nonk562.pG4.pels <- bt.intersect(a = nonk562.pG4, b = k562.ccre %>% filter(V10 %in% c("pELS", "pELS,CTCF-bound")), wa = TRUE) %>% unique()
  
k562.pG4.dels <- bt.intersect(a = k562.pG4, b = k562.ccre %>% filter(V10 %in% c("dELS", "dELS,CTCF-bound")), wa = TRUE) %>% unique()
nonk562.pG4.dels <- bt.intersect(a = nonk562.pG4, b = k562.ccre %>% filter(V10 %in% c("dELS", "dELS,CTCF-bound")), wa = TRUE) %>% unique()

k562.pG4.pls.atac <- GetBW(k562.pG4.pls, "../data/cellline_valid/atac/GSE170378_ENCFF102ARJ_fold_change_over_control_GRCh38.bigWig")
nonk562.pG4.pls.atac <- GetBW(nonk562.pG4.pls, "../data/cellline_valid/atac/GSE170378_ENCFF102ARJ_fold_change_over_control_GRCh38.bigWig")
k562.pG4.pels.atac <- GetBW(k562.pG4.pels, "../data/cellline_valid/atac/GSE170378_ENCFF102ARJ_fold_change_over_control_GRCh38.bigWig")
nonk562.pG4.pels.atac <- GetBW(nonk562.pG4.pels, "../data/cellline_valid/atac/GSE170378_ENCFF102ARJ_fold_change_over_control_GRCh38.bigWig")
k562.pG4.dels.atac <- GetBW(k562.pG4.dels, "../data/cellline_valid/atac/GSE170378_ENCFF102ARJ_fold_change_over_control_GRCh38.bigWig")
nonk562.pG4.dels.atac <- GetBW(nonk562.pG4.dels, "../data/cellline_valid/atac/GSE170378_ENCFF102ARJ_fold_change_over_control_GRCh38.bigWig")

hepg2.pG4.pls <- bt.intersect(a = hepg2.pG4, b = hepg2.ccre %>% filter(V10 %in% c("PLS", "PLS,CTCF-bound")), wa = TRUE) %>% unique()
nonhepg2.pG4.pls <- bt.intersect(a = nonhepg2.pG4, b = hepg2.ccre %>% filter(V10 %in% c("PLS", "PLS,CTCF-bound")), wa = TRUE) %>% unique()

hepg2.pG4.pels <- bt.intersect(a = hepg2.pG4, b = hepg2.ccre %>% filter(V10 %in% c("pELS", "pELS,CTCF-bound")), wa = TRUE) %>% unique()
nonhepg2.pG4.pels <- bt.intersect(a = nonhepg2.pG4, b = hepg2.ccre %>% filter(V10 %in% c("pELS", "pELS,CTCF-bound")), wa = TRUE) %>% unique()
  
hepg2.pG4.dels <- bt.intersect(a = hepg2.pG4, b = hepg2.ccre %>% filter(V10 %in% c("dELS", "dELS,CTCF-bound")), wa = TRUE) %>% unique()
nonhepg2.pG4.dels <- bt.intersect(a = nonhepg2.pG4, b = hepg2.ccre %>% filter(V10 %in% c("dELS", "dELS,CTCF-bound")), wa = TRUE) %>% unique()

hepg2.pG4.pls.atac <- GetBW(hepg2.pG4.pls, "../data/cellline_valid/atac/GSE170251_ENCFF645OHB_fold_change_over_control_GRCh38.bigWig")
nonhepg2.pG4.pls.atac <- GetBW(nonhepg2.pG4.pls, "../data/cellline_valid/atac/GSE170251_ENCFF645OHB_fold_change_over_control_GRCh38.bigWig")
hepg2.pG4.pels.atac <- GetBW(hepg2.pG4.pels, "../data/cellline_valid/atac/GSE170251_ENCFF645OHB_fold_change_over_control_GRCh38.bigWig")
nonhepg2.pG4.pels.atac <- GetBW(nonhepg2.pG4.pels, "../data/cellline_valid/atac/GSE170251_ENCFF645OHB_fold_change_over_control_GRCh38.bigWig")
hepg2.pG4.dels.atac <- GetBW(hepg2.pG4.dels, "../data/cellline_valid/atac/GSE170251_ENCFF645OHB_fold_change_over_control_GRCh38.bigWig")
nonhepg2.pG4.dels.atac <- GetBW(nonhepg2.pG4.dels, "../data/cellline_valid/atac/GSE170251_ENCFF645OHB_fold_change_over_control_GRCh38.bigWig")

k562.atac.score <- bind_rows(
                      data.frame(value = log2(k562.pG4.pls.atac$score + 1), cCRE = "PLS", group = "G4 ChIP-seq supported pG4s", cell = "K562"),
                      data.frame(value = log2(nonk562.pG4.pls.atac$score + 1), cCRE = "PLS", group = "other pG4s", cell = "K562"),
                      data.frame(value = log2(k562.pG4.pels.atac$score + 1), cCRE = "pELS", group = "G4 ChIP-seq supported pG4s", cell = "K562"),
                      data.frame(value = log2(nonk562.pG4.pels.atac$score + 1), cCRE = "pELS", group = "other pG4s", cell = "K562"),
                      data.frame(value = log2(k562.pG4.dels.atac$score + 1), cCRE = "dELS", group = "G4 ChIP-seq supported pG4s", cell = "K562"),
                      data.frame(value = log2(nonk562.pG4.dels.atac$score + 1), cCRE = "dELS", group = "other pG4s", cell = "K562")
                    )

# p val
# 0 0 0
ggsave("../figure/cellline_valid/k562_atac.pdf",
  ggboxplot(k562.atac.score, "cCRE", "value", fill = "group", palette = c("#f05c3b", "#8a9097"), width = 0.5) + 
  rremove("xlab") + rremove("legend.title") + ylab("ATAC signals(log2)"), width = 4, height = 2.5)

hepg2.atac.score <- bind_rows(
                      data.frame(value = log2(hepg2.pG4.pls.atac$score + 1), cCRE = "PLS", group = "G4 ChIP-seq supported pG4s", cell = "HepG2"),
                      data.frame(value = log2(nonhepg2.pG4.pls.atac$score + 1), cCRE = "PLS", group = "other pG4s", cell = "HepG2"),
                      data.frame(value = log2(hepg2.pG4.pels.atac$score + 1), cCRE = "pELS", group = "G4 ChIP-seq supported pG4s", cell = "HepG2"),
                      data.frame(value = log2(nonhepg2.pG4.pels.atac$score + 1), cCRE = "pELS", group = "other pG4s", cell = "HepG2"),
                      data.frame(value = log2(hepg2.pG4.dels.atac$score + 1), cCRE = "dELS", group = "G4 ChIP-seq supported pG4s", cell = "HepG2"),
                      data.frame(value = log2(nonhepg2.pG4.dels.atac$score + 1), cCRE = "dELS", group = "other pG4s", cell = "HepG2")
                    )
# p val
# 0 6.708371e-295 1.943439e-179
ggsave("../figure/cellline_valid/hepg2_atac.pdf",
  ggboxplot(hepg2.atac.score, "cCRE", "value", fill = "group", palette = c("#f05c3b", "#8a9097"), width = 0.5) + 
  rremove("xlab") + rremove("legend.title") + ylab("ATAC signals(log2)"), width = 4, height = 2.5)

#wilcox.test(hepg2.pG4.dels.atac$score, nonhepg2.pG4.dels.atac$score, alternative = "less")



