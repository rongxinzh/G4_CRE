
set.seed(1)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(bedtoolsr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(colorRamp2))
suppressPackageStartupMessages(library(paletteer))
suppressPackageStartupMessages(library(EnrichedHeatmap))

GrEHmap <- function(data = NULL, center = FALSE) {
  data <- data[, 1:3]
  data$score <- 1
  if (center != FALSE) {
  	data[, 2] <- floor((data[, 2] + data[, 3]) / 2)
  	data[, 3] <- data[, 2] + 1
  }
  colnames(data)[1:3] <- c("chr", "start", "end")
  data <- makeGRangesFromDataFrame(data,
                                   keep.extra.columns = TRUE,
                                   ignore.strand = TRUE)

  return(data)
}

signalComp <- function(gr1, gr2, gr3) {
  
  mat <- NULL
  mat.gr1 <- normalizeToMatrix(gr1, gr3, value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
  mat.gr2 <- normalizeToMatrix(gr2, gr3, value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
  gc()
  mat <- bind_rows(mat,
                   data.frame(value = colMeans(data.frame(mat.gr1)), group = "active cCREs", range = c(-(1000/10):(1000/10))),
                   data.frame(value = colMeans(data.frame(mat.gr2)), group = "silent cCREs", range = c(-(1000/10):(1000/10))))

  mat$group <- factor(mat$group, levels = c("active cCREs", "silent cCREs"))

  return(mat)

}

G4 <- fread("../data/G4/G4Hunter_w25_s1.5_hg38.txt", sep = "\t", header = FALSE) %>% data.frame()
ccre <- fread("../data/ccre/V3/G4_cCRE_annotation.txt", sep = "\t", header = TRUE) %>% data.frame() %>% filter(chr %in% paste0("chr", c(1:22, "X")))

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

k562.pG4.gr <- GrEHmap(k562.pG4, center = TRUE)
nonk562.pG4.gr <- GrEHmap(nonk562.pG4, center = TRUE)
hepg2.pG4.gr <- GrEHmap(hepg2.pG4, center = TRUE)
nonhepg2.pG4.gr <- GrEHmap(nonhepg2.pG4, center = TRUE)

k562.pls.gr <- GrEHmap(k562.ccre %>% filter(V10 %in% c("PLS", "PLS,CTCF-bound")), center = FALSE)
k562.pels.gr <- GrEHmap(k562.ccre %>% filter(V10 %in% c("pELS", "pELS,CTCF-bound")), center = FALSE)
k562.dels.gr <- GrEHmap(k562.ccre %>% filter(V10 %in% c("dELS", "dELS,CTCF-bound")), center = FALSE)

hepg2.pls.gr <- GrEHmap(hepg2.ccre %>% filter(V10 %in% c("PLS", "PLS,CTCF-bound")), center = FALSE)
hepg2.pels.gr <- GrEHmap(hepg2.ccre %>% filter(V10 %in% c("pELS", "pELS,CTCF-bound")), center = FALSE)
hepg2.dels.gr <- GrEHmap(hepg2.ccre %>% filter(V10 %in% c("dELS", "dELS,CTCF-bound")), center = FALSE)

k562.nonACT.pls <- ccre %>% filter(encodeLabel == "PLS", !(ID %in% (k562.ccre %>% filter(V10 %in% c("PLS", "PLS,CTCF-bound")) %>% select(V4) %>% unlist() %>% as.vector())))
k562.nonACT.pels <- ccre %>% filter(encodeLabel == "pELS", !(ID %in% (k562.ccre %>% filter(V10 %in% c("pELS", "pELS,CTCF-bound")) %>% select(V4) %>% unlist() %>% as.vector())))
k562.nonACT.dels <- ccre %>% filter(encodeLabel == "dELS", !(ID %in% (k562.ccre %>% filter(V10 %in% c("dELS", "dELS,CTCF-bound")) %>% select(V4) %>% unlist() %>% as.vector())))

hepg2.nonACT.pls <- ccre %>% filter(encodeLabel == "PLS", !(ID %in% (hepg2.ccre %>% filter(V10 %in% c("PLS", "PLS,CTCF-bound")) %>% select(V4) %>% unlist() %>% as.vector())))
hepg2.nonACT.pels <- ccre %>% filter(encodeLabel == "pELS", !(ID %in% (hepg2.ccre %>% filter(V10 %in% c("pELS", "pELS,CTCF-bound")) %>% select(V4) %>% unlist() %>% as.vector())))
hepg2.nonACT.dels <- ccre %>% filter(encodeLabel == "dELS", !(ID %in% (hepg2.ccre %>% filter(V10 %in% c("dELS", "dELS,CTCF-bound")) %>% select(V4) %>% unlist() %>% as.vector())))

k562.nonACT.pls.gr <- GrEHmap(k562.nonACT.pls, center = FALSE)
k562.nonACT.pels.gr <- GrEHmap(k562.nonACT.pels, center = FALSE)
k562.nonACT.dels.gr <- GrEHmap(k562.nonACT.dels, center = FALSE)

hepg2.nonACT.pls.gr <- GrEHmap(hepg2.nonACT.pls, center = FALSE)
hepg2.nonACT.pels.gr <- GrEHmap(hepg2.nonACT.pels, center = FALSE)
hepg2.nonACT.dels.gr <- GrEHmap(hepg2.nonACT.dels, center = FALSE)

k562.pls.mat <- signalComp(k562.pls.gr, k562.nonACT.pls.gr, k562.pG4.gr)
nonk562.pls.mat <- signalComp(k562.pls.gr, k562.nonACT.pls.gr, nonk562.pG4.gr)

k562.pels.mat <- signalComp(k562.pels.gr, k562.nonACT.pels.gr, k562.pG4.gr)
nonk562.pels.mat <- signalComp(k562.pels.gr, k562.nonACT.pels.gr, nonk562.pG4.gr)

k562.dels.mat <- signalComp(k562.dels.gr, k562.nonACT.dels.gr, k562.pG4.gr)
nonk562.dels.mat <- signalComp(k562.dels.gr, k562.nonACT.dels.gr, nonk562.pG4.gr)

hepg2.pls.mat <- signalComp(hepg2.pls.gr, hepg2.nonACT.pls.gr, hepg2.pG4.gr)
nonhepg2.pls.mat <- signalComp(hepg2.pls.gr, hepg2.nonACT.pls.gr, nonhepg2.pG4.gr)

hepg2.pels.mat <- signalComp(hepg2.pels.gr, hepg2.nonACT.pels.gr, hepg2.pG4.gr)
nonhepg2.pels.mat <- signalComp(hepg2.pels.gr, hepg2.nonACT.pels.gr, nonhepg2.pG4.gr)

hepg2.dels.mat <- signalComp(hepg2.dels.gr, hepg2.nonACT.dels.gr, hepg2.pG4.gr)
nonhepg2.dels.mat <- signalComp(hepg2.dels.gr, hepg2.nonACT.dels.gr, nonhepg2.pG4.gr)

save(list=ls(), file = "G4_cCRE_enrich_v2.RData")

all.data <- bind_rows(
  data.frame(k562.pls.mat, cell = "K562", cCRE = "PLS", support = "G4 ChIP-seq"),
  data.frame(nonk562.pls.mat, cell = "K562", cCRE = "PLS", support = "No G4 ChIP-seq"),
  data.frame(k562.pels.mat, cell = "K562", cCRE = "pELS", support = "G4 ChIP-seq"),
  data.frame(nonk562.pels.mat, cell = "K562", cCRE = "pELS", support = "No G4 ChIP-seq"),
  data.frame(k562.dels.mat, cell = "K562", cCRE = "dELS", support = "G4 ChIP-seq"),
  data.frame(nonk562.dels.mat, cell = "K562", cCRE = "dELS", support = "No G4 ChIP-seq"),
  data.frame(hepg2.pls.mat, cell = "HepG2", cCRE = "PLS", support = "G4 ChIP-seq"),
  data.frame(nonhepg2.pls.mat, cell = "HepG2", cCRE = "PLS", support = "No G4 ChIP-seq"),
  data.frame(hepg2.pels.mat, cell = "HepG2", cCRE = "pELS", support = "G4 ChIP-seq"),
  data.frame(nonhepg2.pels.mat, cell = "HepG2", cCRE = "pELS", support = "No G4 ChIP-seq"),
  data.frame(hepg2.dels.mat, cell = "HepG2", cCRE = "dELS", support = "G4 ChIP-seq"),
  data.frame(nonhepg2.dels.mat, cell = "HepG2", cCRE = "dELS", support = "No G4 ChIP-seq")
)

all.data$group <- factor(all.data$group, levels = unique(all.data$group))
all.data$cCRE <- factor(all.data$cCRE, levels = unique(all.data$cCRE))
all.data$support <- factor(all.data$support, levels = unique(all.data$support))

p1 <- ggplot(all.data %>% filter(cell == "K562"), aes(x = range, y = value, group = group, color = group)) +
      geom_line() +
      geom_vline(xintercept = 0, linetype = "dotted", color = "grey", linewidth = 0.35) + 
      theme_classic() + rremove("xlab") + theme(legend.position="top", legend.direction = "vertical") + 
      ylab("cCRE density") + rremove("legend.title") + facet_wrap(vars(support, cCRE), scales = "free") + 
      scale_color_manual(values = paletteer_d("ggsci::default_aaas"))

p2 <- ggplot(all.data %>% filter(cell == "HepG2"), aes(x = range, y = value, group = group, color = group)) +
      geom_line() +
      geom_vline(xintercept = 0, linetype = "dotted", color = "grey", linewidth = 0.35) + 
      theme_classic() + rremove("xlab") + theme(legend.position="top", legend.direction = "vertical") + 
      ylab("cCRE density") + rremove("legend.title") + facet_wrap(vars(support, cCRE), scales = "free") + 
      scale_color_manual(values = paletteer_d("ggsci::default_aaas"))

ggsave(filename = "../figure/cellline_valid/K562_pG4_cCRE.pdf", p1, device = "pdf", width = 8, height = 5)
ggsave(filename = "../figure/cellline_valid/HepG2_pG4_cCRE.pdf", p2, device = "pdf", width = 8, height = 5)



p1 <- ggplot(k562.pls.mat, aes(x = range, y = value, group = group, color = group)) +
      geom_line() +
      geom_vline(xintercept = 0, linetype = "dotted", color = "grey", linewidth = 0.35) + 
      theme_classic() + rremove("xlab") + theme(legend.position="top", legend.direction = "vertical") + 
      ylab("cCRE density") + rremove("legend.title") +
      scale_color_manual(values = paletteer_d("ggsci::default_aaas"))

p2 <- ggplot(nonk562.pls.mat, aes(x = range, y = value, group = group, color = group)) +
      geom_line() +
      geom_vline(xintercept = 0, linetype = "dotted", color = "grey", linewidth = 0.35) + 
      theme_classic() + rremove("xlab") + theme(legend.position="top", legend.direction = "vertical") + 
      ylab("cCRE density") + rremove("legend.title") +
      scale_color_manual(values = paletteer_d("ggsci::default_aaas"))

p3 <- ggplot(k562.pels.mat, aes(x = range, y = value, group = group, color = group)) +
      geom_line() +
      geom_vline(xintercept = 0, linetype = "dotted", color = "grey", linewidth = 0.35) + 
      theme_classic() + rremove("xlab") + theme(legend.position="top", legend.direction = "vertical") + 
      ylab("cCRE density") + rremove("legend.title") +
      scale_color_manual(values = paletteer_d("ggsci::default_aaas"))

p4 <- ggplot(nonk562.pels.mat, aes(x = range, y = value, group = group, color = group)) +
      geom_line() +
      geom_vline(xintercept = 0, linetype = "dotted", color = "grey", linewidth = 0.35) + 
      theme_classic() + rremove("xlab") + theme(legend.position="top", legend.direction = "vertical") + 
      ylab("cCRE density") + rremove("legend.title") +
      scale_color_manual(values = paletteer_d("ggsci::default_aaas"))

p5 <- ggplot(k562.dels.mat, aes(x = range, y = value, group = group, color = group)) +
      geom_line() +
      geom_vline(xintercept = 0, linetype = "dotted", color = "grey", linewidth = 0.35) + 
      theme_classic() + rremove("xlab") + theme(legend.position="top", legend.direction = "vertical") + 
      ylab("cCRE density") + rremove("legend.title") +
      scale_color_manual(values = paletteer_d("ggsci::default_aaas"))

p6 <- ggplot(nonk562.dels.mat, aes(x = range, y = value, group = group, color = group)) +
      geom_line() +
      geom_vline(xintercept = 0, linetype = "dotted", color = "grey", linewidth = 0.35) + 
      theme_classic() + rremove("xlab") + theme(legend.position="top", legend.direction = "vertical") + 
      ylab("cCRE density") + rremove("legend.title") +
      scale_color_manual(values = paletteer_d("ggsci::default_aaas"))

p7 <- ggplot(hepg2.pls.mat, aes(x = range, y = value, group = group, color = group)) +
      geom_line() +
      geom_vline(xintercept = 0, linetype = "dotted", color = "grey", linewidth = 0.35) + 
      theme_classic() + rremove("xlab") + theme(legend.position="top", legend.direction = "vertical") + 
      ylab("cCRE density") + rremove("legend.title") +
      scale_color_manual(values = paletteer_d("ggsci::default_aaas"))

p8 <- ggplot(nonhepg2.pls.mat, aes(x = range, y = value, group = group, color = group)) +
      geom_line() +
      geom_vline(xintercept = 0, linetype = "dotted", color = "grey", linewidth = 0.35) + 
      theme_classic() + rremove("xlab") + theme(legend.position="top", legend.direction = "vertical") + 
      ylab("cCRE density") + rremove("legend.title") +
      scale_color_manual(values = paletteer_d("ggsci::default_aaas"))

p9 <- ggplot(hepg2.pels.mat, aes(x = range, y = value, group = group, color = group)) +
      geom_line() +
      geom_vline(xintercept = 0, linetype = "dotted", color = "grey", linewidth = 0.35) + 
      theme_classic() + rremove("xlab") + theme(legend.position="top", legend.direction = "vertical") + 
      ylab("cCRE density") + rremove("legend.title") +
      scale_color_manual(values = paletteer_d("ggsci::default_aaas"))

p10 <- ggplot(nonhepg2.pels.mat, aes(x = range, y = value, group = group, color = group)) +
      geom_line() +
      geom_vline(xintercept = 0, linetype = "dotted", color = "grey", linewidth = 0.35) + 
      theme_classic() + rremove("xlab") + theme(legend.position="top", legend.direction = "vertical") + 
      ylab("cCRE density") + rremove("legend.title") +
      scale_color_manual(values = paletteer_d("ggsci::default_aaas"))

p11 <- ggplot(hepg2.dels.mat, aes(x = range, y = value, group = group, color = group)) +
      geom_line() +
      geom_vline(xintercept = 0, linetype = "dotted", color = "grey", linewidth = 0.35) + 
      theme_classic() + rremove("xlab") + theme(legend.position="top", legend.direction = "vertical") + 
      ylab("cCRE density") + rremove("legend.title") +
      scale_color_manual(values = paletteer_d("ggsci::default_aaas"))

p12 <- ggplot(nonhepg2.dels.mat, aes(x = range, y = value, group = group, color = group)) +
      geom_line() +
      geom_vline(xintercept = 0, linetype = "dotted", color = "grey", linewidth = 0.35) + 
      theme_classic() + rremove("xlab") + theme(legend.position="top", legend.direction = "vertical") + 
      ylab("cCRE density") + rremove("legend.title") +
      scale_color_manual(values = paletteer_d("ggsci::default_aaas"))

ggsave(filename = "dqw1212dqw.pdf", ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, nrow = 2, ncol = 6), device = "pdf", width = 15, height = 5)








