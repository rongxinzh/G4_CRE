
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

G4 <- fread("../data/G4/G4Hunter_w25_s1.5_hg38.txt", sep = "\t", header = FALSE) %>% data.frame()

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

tf.ensg <- fread("../data/cellline_valid/tf/TFs_Ensembl_v_1.01.txt", sep = "\t", header = FALSE) %>% unlist()
tf.entrez.ensg <- fread("../data/cellline_valid/tf/tf_ensg.txt", sep = "\t", header = FALSE) %>% data.frame()
tf.entrez.ensg <- tf.entrez.ensg %>% filter(V2 %in% tf.ensg)

tfbs <- fread("../data/cellline_valid/tf/encRegTfbsClusteredWithCells.hg38.bed", sep = "\t", header = FALSE) %>% data.frame()
tfbs <- tfbs %>% filter(V4 %in% tf.entrez.ensg$V1)
k562.tfbs <- tfbs %>% filter(str_detect(V6, "K562"))
k562.tf <- unique(k562.tfbs[, 4])
hepg2.tfbs <- tfbs %>% filter(str_detect(V6, "HepG2"))
hepg2.tf <- unique(hepg2.tfbs[, 4])
tf.use <- intersect(hepg2.tf, k562.tf)

k562.mat <- NULL
hepg2.mat <- NULL

for (tmp.tf in tf.use) {

  message(tmp.tf)

  tmp.k562.tfbs <- k562.tfbs %>% filter(V4 == tmp.tf)
  # pls, pels, dels
  k562.pG4.pls <- bt.intersect(a = k562.pG4, b = k562.ccre %>% filter(V10 %in% c("PLS", "PLS,CTCF-bound")), wa = TRUE) %>% unique()
  nonk562.pG4.pls <- bt.intersect(a = nonk562.pG4, b = k562.ccre %>% filter(V10 %in% c("PLS", "PLS,CTCF-bound")), wa = TRUE) %>% unique()
  tmp.k562.pG4.tfbs.pls.cvg <- bt.coverage(a = k562.pG4.pls, b = tmp.k562.tfbs)
  tmp.nonk562.pG4.tfbs.pls.cvg <- bt.coverage(a = nonk562.pG4.pls, b = tmp.k562.tfbs)

  k562.pG4.pels <- bt.intersect(a = k562.pG4, b = k562.ccre %>% filter(V10 %in% c("pELS", "pELS,CTCF-bound")), wa = TRUE) %>% unique()
  nonk562.pG4.pels <- bt.intersect(a = nonk562.pG4, b = k562.ccre %>% filter(V10 %in% c("pELS", "pELS,CTCF-bound")), wa = TRUE) %>% unique()
  tmp.k562.pG4.tfbs.pels.cvg <- bt.coverage(a = k562.pG4.pels, b = tmp.k562.tfbs)
  tmp.nonk562.pG4.tfbs.pels.cvg <- bt.coverage(a = nonk562.pG4.pels, b = tmp.k562.tfbs)

  k562.pG4.dels <- bt.intersect(a = k562.pG4, b = k562.ccre %>% filter(V10 %in% c("dELS", "dELS,CTCF-bound")), wa = TRUE) %>% unique()
  nonk562.pG4.dels <- bt.intersect(a = nonk562.pG4, b = k562.ccre %>% filter(V10 %in% c("dELS", "dELS,CTCF-bound")), wa = TRUE) %>% unique()
  tmp.k562.pG4.tfbs.dels.cvg <- bt.coverage(a = k562.pG4.dels, b = tmp.k562.tfbs)
  tmp.nonk562.pG4.tfbs.dels.cvg <- bt.coverage(a = nonk562.pG4.dels, b = tmp.k562.tfbs)

  k562.mat <- bind_rows(k562.mat, 
    data.frame(tf = tmp.tf, cCRE = "PLS", v1 = mean(tmp.k562.pG4.tfbs.pls.cvg[, 13]), v2 = mean(tmp.nonk562.pG4.tfbs.pls.cvg[, 13])), 
    data.frame(tf = tmp.tf, cCRE = "pELS", v1 = mean(tmp.k562.pG4.tfbs.pels.cvg[, 13]), v2 = mean(tmp.nonk562.pG4.tfbs.pels.cvg[, 13])), 
    data.frame(tf = tmp.tf, cCRE = "dELS", v1 = mean(tmp.k562.pG4.tfbs.dels.cvg[, 13]), v2 = mean(tmp.nonk562.pG4.tfbs.dels.cvg[, 13]))
    )

  tmp.hepg2.tfbs <- hepg2.tfbs %>% filter(V4 == tmp.tf)
  # pls, pels, dels
  hepg2.pG4.pls <- bt.intersect(a = hepg2.pG4, b = hepg2.ccre %>% filter(V10 %in% c("PLS", "PLS,CTCF-bound")), wa = TRUE) %>% unique()
  nonhepg2.pG4.pls <- bt.intersect(a = nonhepg2.pG4, b = hepg2.ccre %>% filter(V10 %in% c("PLS", "PLS,CTCF-bound")), wa = TRUE) %>% unique()
  tmp.hepg2.pG4.tfbs.pls.cvg <- bt.coverage(a = hepg2.pG4.pls, b = tmp.hepg2.tfbs)
  tmp.nonhepg2.pG4.tfbs.pls.cvg <- bt.coverage(a = nonhepg2.pG4.pls, b = tmp.hepg2.tfbs)

  hepg2.pG4.pels <- bt.intersect(a = hepg2.pG4, b = hepg2.ccre %>% filter(V10 %in% c("pELS", "pELS,CTCF-bound")), wa = TRUE) %>% unique()
  nonhepg2.pG4.pels <- bt.intersect(a = nonhepg2.pG4, b = hepg2.ccre %>% filter(V10 %in% c("pELS", "pELS,CTCF-bound")), wa = TRUE) %>% unique()
  tmp.hepg2.pG4.tfbs.pels.cvg <- bt.coverage(a = hepg2.pG4.pels, b = tmp.hepg2.tfbs)
  tmp.nonhepg2.pG4.tfbs.pels.cvg <- bt.coverage(a = nonhepg2.pG4.pels, b = tmp.hepg2.tfbs)

  hepg2.pG4.dels <- bt.intersect(a = hepg2.pG4, b = hepg2.ccre %>% filter(V10 %in% c("dELS", "dELS,CTCF-bound")), wa = TRUE) %>% unique()
  nonhepg2.pG4.dels <- bt.intersect(a = nonhepg2.pG4, b = hepg2.ccre %>% filter(V10 %in% c("dELS", "dELS,CTCF-bound")), wa = TRUE) %>% unique()
  tmp.hepg2.pG4.tfbs.dels.cvg <- bt.coverage(a = hepg2.pG4.dels, b = tmp.hepg2.tfbs)
  tmp.nonhepg2.pG4.tfbs.dels.cvg <- bt.coverage(a = nonhepg2.pG4.dels, b = tmp.hepg2.tfbs)

  hepg2.mat <- bind_rows(hepg2.mat, 
    data.frame(tf = tmp.tf, cCRE = "PLS", v1 = mean(tmp.hepg2.pG4.tfbs.pls.cvg[, 13]), v2 = mean(tmp.nonhepg2.pG4.tfbs.pls.cvg[, 13])), 
    data.frame(tf = tmp.tf, cCRE = "pELS", v1 = mean(tmp.hepg2.pG4.tfbs.pels.cvg[, 13]), v2 = mean(tmp.nonhepg2.pG4.tfbs.pels.cvg[, 13])), 
    data.frame(tf = tmp.tf, cCRE = "dELS", v1 = mean(tmp.hepg2.pG4.tfbs.dels.cvg[, 13]), v2 = mean(tmp.nonhepg2.pG4.tfbs.dels.cvg[, 13]))
    )

}

save(list=ls(), file = "G4_cCRE_enrich_tf.RData")
load("G4_cCRE_enrich_tf.RData")


k562.mat$col <- ifelse(k562.mat$v1 >= k562.mat$v2, "A", "B")
k562.mat$col <- ifelse((k562.mat$v1 >= k562.mat$v2) & (k562.mat$v1 > 0.5), "C", k562.mat$col)
k562.mat$col <- ifelse((k562.mat$v1 < k562.mat$v2) & (k562.mat$v2 > 0.5), "D", k562.mat$col)

p1 <- ggscatter(k562.mat %>% filter(cCRE == "PLS"), x = "v1", y = "v2", color = "col", palette = c("#F39C12", "#B03A2E", "#2E86C1"), label = "tf", 
        repel = TRUE, label.select = k562.mat %>% filter(cCRE == "PLS") %>% filter(col == "C" | col == "D") %>% select(tf) %>% unlist()) + rremove("legend") + 
      xlim(0, 1) + ylim(0, 1) + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") + rremove("xlab") + rremove("ylab")

p2 <- ggscatter(k562.mat %>% filter(cCRE == "pELS"), x = "v1", y = "v2", color = "col", palette = c("#F39C12", "#B03A2E", "#2E86C1"), label = "tf", 
        repel = TRUE, label.select = k562.mat %>% filter(cCRE == "pELS") %>% filter(col == "C" | col == "D") %>% select(tf) %>% unlist()) + rremove("legend") + 
      xlim(0, 1) + ylim(0, 1) + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") + rremove("xlab") + rremove("ylab")

p3 <- ggscatter(k562.mat %>% filter(cCRE == "dELS"), x = "v1", y = "v2", color = "col", palette = c("#F39C12", "#239B56", "#B03A2E", "#2E86C1"), label = "tf", 
        repel = TRUE, label.select = k562.mat %>% filter(cCRE == "dELS") %>% filter(col == "C" | col == "D") %>% select(tf) %>% unlist()) + rremove("legend") + 
      xlim(0, 1) + ylim(0, 1) + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") + rremove("xlab") + rremove("ylab")

hepg2.mat$col <- ifelse(hepg2.mat$v1 >= hepg2.mat$v2, "A", "B")
hepg2.mat$col <- ifelse((hepg2.mat$v1 >= hepg2.mat$v2) & (hepg2.mat$v1 > 0.5), "C", hepg2.mat$col)
hepg2.mat$col <- ifelse((hepg2.mat$v1 < hepg2.mat$v2) & (hepg2.mat$v2 > 0.5), "D", hepg2.mat$col)

p4 <- ggscatter(hepg2.mat %>% filter(cCRE == "PLS"), x = "v1", y = "v2", color = "col", palette = c("#F39C12", "#B03A2E", "#2E86C1"), label = "tf", 
        repel = TRUE, label.select = hepg2.mat %>% filter(cCRE == "PLS") %>% filter(col == "C" | col == "D") %>% select(tf) %>% unlist()) + rremove("legend") + 
      xlim(0, 1) + ylim(0, 1) + geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.8, color = "grey") + rremove("xlab") + rremove("ylab")

p5 <- ggscatter(hepg2.mat %>% filter(cCRE == "pELS"), x = "v1", y = "v2", color = "col", palette = c("#F39C12", "#239B56", "#B03A2E", "#2E86C1"), label = "tf", 
        repel = TRUE, label.select = hepg2.mat %>% filter(cCRE == "pELS") %>% filter(col == "C" | col == "D") %>% select(tf) %>% unlist()) + rremove("legend") + 
      xlim(0, 1) + ylim(0, 1) + geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.8, color = "grey") + rremove("xlab") + rremove("ylab")

p6 <- ggscatter(hepg2.mat %>% filter(cCRE == "dELS"), x = "v1", y = "v2", color = "col", palette = c("#F39C12", "#B03A2E", "#2E86C1"), label = "tf", 
        repel = TRUE, label.select = hepg2.mat %>% filter(cCRE == "dELS") %>% filter(col == "C" | col == "D") %>% select(tf) %>% unlist()) + rremove("legend") + 
      xlim(0, 1) + ylim(0, 1) + geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.8, color = "grey") + rremove("xlab") + rremove("ylab")

ggsave("../figure/cellline_valid/k562_tf.pdf", ggarrange(p1, p2, p3, ncol = 3), width = 8, height = 2.5)
ggsave("../figure/cellline_valid/hepg2_tf.pdf", ggarrange(p4, p5, p6, ncol = 3), width = 8, height = 2.5)
