
set.seed(1)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(bedtoolsr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(colorRamp2))
suppressPackageStartupMessages(library(paletteer))

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
#nonk562.pG4 <- bt.intersect(a = G4, b = k562.G4, wa = TRUE, v = TRUE) %>% unique()
hepg2.pG4 <- bt.intersect(a = G4, b = hepg2.G4, wa = TRUE) %>% unique()
#nonhepg2.pG4 <- bt.intersect(a = G4, b = hepg2.G4, wa = TRUE, v = TRUE) %>% unique()

k562.pls <- k562.ccre %>% filter(V10 %in% c("PLS", "PLS,CTCF-bound"))
k562.pels <- k562.ccre %>% filter(V10 %in% c("pELS", "pELS,CTCF-bound"))
k562.dels <- k562.ccre %>% filter(V10 %in% c("dELS", "dELS,CTCF-bound"))

hepg2.pls <- hepg2.ccre %>% filter(V10 %in% c("PLS", "PLS,CTCF-bound"))
hepg2.pels <- hepg2.ccre %>% filter(V10 %in% c("pELS", "pELS,CTCF-bound"))
hepg2.dels <- hepg2.ccre %>% filter(V10 %in% c("dELS", "dELS,CTCF-bound"))

k562.silent.pls <- ccre %>% filter(encodeLabel == "PLS", !(ID %in% (k562.ccre %>% filter(V10 %in% c("PLS", "PLS,CTCF-bound")) %>% select(V4) %>% unlist() %>% as.vector())))
k562.silent.pels <- ccre %>% filter(encodeLabel == "pELS", !(ID %in% (k562.ccre %>% filter(V10 %in% c("pELS", "pELS,CTCF-bound")) %>% select(V4) %>% unlist() %>% as.vector())))
k562.silent.dels <- ccre %>% filter(encodeLabel == "dELS", !(ID %in% (k562.ccre %>% filter(V10 %in% c("dELS", "dELS,CTCF-bound")) %>% select(V4) %>% unlist() %>% as.vector())))

hepg2.silent.pls <- ccre %>% filter(encodeLabel == "PLS", !(ID %in% (hepg2.ccre %>% filter(V10 %in% c("PLS", "PLS,CTCF-bound")) %>% select(V4) %>% unlist() %>% as.vector())))
hepg2.silent.pels <- ccre %>% filter(encodeLabel == "pELS", !(ID %in% (hepg2.ccre %>% filter(V10 %in% c("pELS", "pELS,CTCF-bound")) %>% select(V4) %>% unlist() %>% as.vector())))
hepg2.silent.dels <- ccre %>% filter(encodeLabel == "dELS", !(ID %in% (hepg2.ccre %>% filter(V10 %in% c("dELS", "dELS,CTCF-bound")) %>% select(V4) %>% unlist() %>% as.vector())))

k562.pls.p1 <- bt.intersect(a = k562.pG4, b = k562.pls, wa = TRUE) %>% unique() %>% nrow()
k562.pls.p2 <- bt.intersect(a = k562.pG4, b = k562.silent.pls, wa = TRUE) %>% unique() %>% nrow()
k562.pels.p1 <- bt.intersect(a = k562.pG4, b = k562.pels, wa = TRUE) %>% unique() %>% nrow()
k562.pels.p2 <- bt.intersect(a = k562.pG4, b = k562.silent.pels, wa = TRUE) %>% unique() %>% nrow()
k562.dels.p1 <- bt.intersect(a = k562.pG4, b = k562.dels, wa = TRUE) %>% unique() %>% nrow()
k562.dels.p2 <- bt.intersect(a = k562.pG4, b = k562.silent.dels, wa = TRUE) %>% unique() %>% nrow()

hepg2.pls.p1 <- bt.intersect(a = hepg2.pG4, b = hepg2.pls, wa = TRUE) %>% unique() %>% nrow()
hepg2.pls.p2 <- bt.intersect(a = hepg2.pG4, b = hepg2.silent.pls, wa = TRUE) %>% unique() %>% nrow()
hepg2.pels.p1 <- bt.intersect(a = hepg2.pG4, b = hepg2.pels, wa = TRUE) %>% unique() %>% nrow()
hepg2.pels.p2 <- bt.intersect(a = hepg2.pG4, b = hepg2.silent.pels, wa = TRUE) %>% unique() %>% nrow()
hepg2.dels.p1 <- bt.intersect(a = hepg2.pG4, b = hepg2.dels, wa = TRUE) %>% unique() %>% nrow()
hepg2.dels.p2 <- bt.intersect(a = hepg2.pG4, b = hepg2.silent.dels, wa = TRUE) %>% unique() %>% nrow()

k.pls.p <- k562.pls.p1 / (k562.pls.p1 + k562.pls.p2)
k.pels.p <- k562.pels.p1 / (k562.pels.p1 + k562.pels.p2)
k.dels.p <- k562.dels.p1 / (k562.dels.p1 + k562.dels.p2)

h.pls.p <- hepg2.pls.p1 / (hepg2.pls.p1 + hepg2.pls.p2)
h.pels.p <- hepg2.pels.p1 / (hepg2.pels.p1 + hepg2.pels.p2)
h.dels.p <- hepg2.dels.p1 / (hepg2.dels.p1 + hepg2.dels.p2)

# k562
k.pls.rand <- NULL
k.pels.rand <- NULL
k.dels.rand <- NULL

for (tmp.round in 1:1000) {
  
  message(tmp.round)
  tmp.G4 <- G4[sample(1:nrow(G4), nrow(k562.pG4)), ]

  tmp.pls.p1 <- bt.intersect(a = tmp.G4, b = k562.pls, wa = TRUE) %>% unique() %>% nrow()
  tmp.pls.p2 <- bt.intersect(a = tmp.G4, b = k562.silent.pls, wa = TRUE) %>% unique() %>% nrow()
  tmp.pels.p1 <- bt.intersect(a = tmp.G4, b = k562.pels, wa = TRUE) %>% unique() %>% nrow()
  tmp.pels.p2 <- bt.intersect(a = tmp.G4, b = k562.silent.pels, wa = TRUE) %>% unique() %>% nrow()
  tmp.dels.p1 <- bt.intersect(a = tmp.G4, b = k562.dels, wa = TRUE) %>% unique() %>% nrow()
  tmp.dels.p2 <- bt.intersect(a = tmp.G4, b = k562.silent.dels, wa = TRUE) %>% unique() %>% nrow()

  k.pls.rand <- c(k.pls.rand, tmp.pls.p1 / (tmp.pls.p1 + tmp.pls.p2))
  k.pels.rand <- c(k.pels.rand, tmp.pels.p1 / (tmp.pels.p1 + tmp.pels.p2))
  k.dels.rand <- c(k.dels.rand, tmp.dels.p1 / (tmp.dels.p1 + tmp.dels.p2))

}
message("K562 done!")

# hepg2
h.pls.rand <- NULL
h.pels.rand <- NULL
h.dels.rand <- NULL

for (tmp.round in 1:1000) {
  
  message(tmp.round)
  tmp.G4 <- G4[sample(1:nrow(G4), nrow(hepg2.pG4)), ]

  tmp.pls.p1 <- bt.intersect(a = tmp.G4, b = hepg2.pls, wa = TRUE) %>% unique() %>% nrow()
  tmp.pls.p2 <- bt.intersect(a = tmp.G4, b = hepg2.silent.pls, wa = TRUE) %>% unique() %>% nrow()
  tmp.pels.p1 <- bt.intersect(a = tmp.G4, b = hepg2.pels, wa = TRUE) %>% unique() %>% nrow()
  tmp.pels.p2 <- bt.intersect(a = tmp.G4, b = hepg2.silent.pels, wa = TRUE) %>% unique() %>% nrow()
  tmp.dels.p1 <- bt.intersect(a = tmp.G4, b = hepg2.dels, wa = TRUE) %>% unique() %>% nrow()
  tmp.dels.p2 <- bt.intersect(a = tmp.G4, b = hepg2.silent.dels, wa = TRUE) %>% unique() %>% nrow()

  h.pls.rand <- c(h.pls.rand, tmp.pls.p1 / (tmp.pls.p1 + tmp.pls.p2))
  h.pels.rand <- c(h.pels.rand, tmp.pels.p1 / (tmp.pels.p1 + tmp.pels.p2))
  h.dels.rand <- c(h.dels.rand, tmp.dels.p1 / (tmp.dels.p1 + tmp.dels.p2))

}
message("HepG2 done!")

save(list=ls(), file = "G4_cCRE_enrich_fold.RData")

mat <- bind_rows(
  data.frame(fold = k.pls.rand, group = "PLS", cell = "K562"),
  data.frame(fold = k.pels.rand, group = "pELS", cell = "K562"),
  data.frame(fold = k.dels.rand, group = "dELS", cell = "K562"),
  data.frame(fold = h.pls.rand, group = "PLS", cell = "HepG2"),
  data.frame(fold = h.pels.rand, group = "pELS", cell = "HepG2"),
  data.frame(fold = h.dels.rand, group = "dELS", cell = "HepG2")
)
mat$group <- factor(mat$group, levels = unique(mat$group))
mat$cell <- factor(mat$cell, levels = unique(mat$cell))

v.mat <- bind_rows(
  data.frame(fold = k.pls.p, group = "PLS", cell = "K562"),
  data.frame(fold = k.pels.p, group = "pELS", cell = "K562"),
  data.frame(fold = k.dels.p, group = "dELS", cell = "K562"),
  data.frame(fold = h.pls.p, group = "PLS", cell = "HepG2"),
  data.frame(fold = h.pels.p, group = "pELS", cell = "HepG2"),
  data.frame(fold = h.dels.p, group = "dELS", cell = "HepG2")
)

pdf("../figure/cellline_valid/G4_proportion_density_fold.pdf", height = 3, width = 6)
ggplot(mat, aes(x=fold, fill=group)) + ylab('Density') +
  geom_density()  +
  geom_vline(data = v.mat, aes(xintercept = fold), linetype = 'dashed', show.legend = NA) + 
  facet_wrap(.~cell, scales = "free_x") +
  scale_fill_manual(values = c("#e64b35", "#4ebad5", "#049f87", "#3c5488", "#f49b7f", "#7d6148")) +
  rremove("legend.title") + xlab("Proportion of pG4s located in activated cCREs")
dev.off()

#mat <- bind_rows(
#  data.frame(fold = k.pls.p / k.pls.rand, group = "PLS", cell = "K562"),
#  data.frame(fold = k.pels.p / k.pels.rand, group = "pELS", cell = "K562"),
#  data.frame(fold = k.dels.p / k.dels.rand, group = "dELS", cell = "K562"),
#  data.frame(fold = h.pls.p / h.pls.rand, group = "PLS", cell = "HepG2"),
#  data.frame(fold = h.pels.p / h.pels.rand, group = "pELS", cell = "HepG2"),
#  data.frame(fold = h.dels.p / h.dels.rand, group = "dELS", cell = "HepG2")
#)
#mat$group <- factor(mat$group, levels = unique(mat$group))
#mat$cell <- factor(mat$cell, levels = unique(mat$cell))

#pdf("../figure/cellline_valid/odds_ratio.pdf", height = 3, width = 6)
#ggplot(mat, aes(x=fold, fill=group)) + ylab('Density') +
#  geom_density(alpha=.7)  + facet_wrap(.~cell, scales = "free_x") +
#  scale_fill_manual(values = c("#e64b35", "#4ebad5", "#049f87", "#3c5488", "#f49b7f", "#7d6148")) +
#  rremove("legend.title") + xlab("Odds ratio")
#dev.off()


