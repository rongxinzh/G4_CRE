
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

k.pls <- NULL
k.pels <- NULL
k.dels <- NULL

for (tmp.round in 1:1000) {
  
  message(tmp.round)
  tmp.G4 <- k562.pG4[sample(1:nrow(k562.pG4), 1000), ]

  tmp.pls.p1 <- bt.intersect(a = tmp.G4, b = k562.pls, wa = TRUE) %>% unique() %>% nrow()
  tmp.pls.p2 <- bt.intersect(a = tmp.G4, b = k562.silent.pls, wa = TRUE) %>% unique() %>% nrow()
  tmp.pels.p1 <- bt.intersect(a = tmp.G4, b = k562.pels, wa = TRUE) %>% unique() %>% nrow()
  tmp.pels.p2 <- bt.intersect(a = tmp.G4, b = k562.silent.pels, wa = TRUE) %>% unique() %>% nrow()
  tmp.dels.p1 <- bt.intersect(a = tmp.G4, b = k562.dels, wa = TRUE) %>% unique() %>% nrow()
  tmp.dels.p2 <- bt.intersect(a = tmp.G4, b = k562.silent.dels, wa = TRUE) %>% unique() %>% nrow()

  k.pls <- c(k.pls, tmp.pls.p1 / (tmp.pls.p1 + tmp.pls.p2))
  k.pels <- c(k.pels, tmp.pels.p1 / (tmp.pels.p1 + tmp.pels.p2))
  k.dels <- c(k.dels, tmp.dels.p1 / (tmp.dels.p1 + tmp.dels.p2))

}
message("1.K562 done!")

# hepg2
h.pls <- NULL
h.pels <- NULL
h.dels <- NULL

for (tmp.round in 1:1000) {
  
  message(tmp.round)
  tmp.G4 <- hepg2.pG4[sample(1:nrow(hepg2.pG4), 1000), ]

  tmp.pls.p1 <- bt.intersect(a = tmp.G4, b = hepg2.pls, wa = TRUE) %>% unique() %>% nrow()
  tmp.pls.p2 <- bt.intersect(a = tmp.G4, b = hepg2.silent.pls, wa = TRUE) %>% unique() %>% nrow()
  tmp.pels.p1 <- bt.intersect(a = tmp.G4, b = hepg2.pels, wa = TRUE) %>% unique() %>% nrow()
  tmp.pels.p2 <- bt.intersect(a = tmp.G4, b = hepg2.silent.pels, wa = TRUE) %>% unique() %>% nrow()
  tmp.dels.p1 <- bt.intersect(a = tmp.G4, b = hepg2.dels, wa = TRUE) %>% unique() %>% nrow()
  tmp.dels.p2 <- bt.intersect(a = tmp.G4, b = hepg2.silent.dels, wa = TRUE) %>% unique() %>% nrow()

  h.pls <- c(h.pls, tmp.pls.p1 / (tmp.pls.p1 + tmp.pls.p2))
  h.pels <- c(h.pels, tmp.pels.p1 / (tmp.pels.p1 + tmp.pels.p2))
  h.dels <- c(h.dels, tmp.dels.p1 / (tmp.dels.p1 + tmp.dels.p2))

}
message("1.HepG2 done!")

# k562
k.pls.rand <- NULL
k.pels.rand <- NULL
k.dels.rand <- NULL

for (tmp.round in 1:1000) {
  
  message(tmp.round)
  tmp.G4 <- G4[sample(1:nrow(G4), 1000), ]

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
message("2.K562 done!")

# hepg2
h.pls.rand <- NULL
h.pels.rand <- NULL
h.dels.rand <- NULL

for (tmp.round in 1:1000) {
  
  message(tmp.round)
  tmp.G4 <- G4[sample(1:nrow(G4), 1000), ]

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
message("2.HepG2 done!")

save(list=ls(), file = "G4_cCRE_enrich_fold_2.RData")

mat <- bind_rows(
  data.frame(fold = k.pls, group = "PLS", cell = "K562", anno = "ChIP-seq supported G4s"),
  data.frame(fold = k.pels, group = "pELS", cell = "K562", anno = "ChIP-seq supported G4s"),
  data.frame(fold = k.dels, group = "dELS", cell = "K562", anno = "ChIP-seq supported G4s"),
  data.frame(fold = h.pls, group = "PLS", cell = "HepG2", anno = "ChIP-seq supported G4s"),
  data.frame(fold = h.pels, group = "pELS", cell = "HepG2", anno = "ChIP-seq supported G4s"),
  data.frame(fold = h.dels, group = "dELS", cell = "HepG2", anno = "ChIP-seq supported G4s"),
  data.frame(fold = k.pls.rand, group = "PLS", cell = "K562", anno = "All G4s"),
  data.frame(fold = k.pels.rand, group = "pELS", cell = "K562", anno = "All G4s"),
  data.frame(fold = k.dels.rand, group = "dELS", cell = "K562", anno = "All G4s"),
  data.frame(fold = h.pls.rand, group = "PLS", cell = "HepG2", anno = "All G4s"),
  data.frame(fold = h.pels.rand, group = "pELS", cell = "HepG2", anno = "All G4s"),
  data.frame(fold = h.dels.rand, group = "dELS", cell = "HepG2", anno = "All G4s")
)
mat$group <- factor(mat$group, levels = unique(mat$group))
mat$cell <- factor(mat$cell, levels = unique(mat$cell))
mat$anno <- factor(mat$anno, levels = unique(mat$anno))

pdf("../figure/cellline_valid/G4_proportion_density_fold.pdf", height = 6, width = 5)
ggplot(mat, aes(x = group, y = fold, fill = anno)) + ylab('Proportion of G4s in activated cCREs') +
  geom_boxplot() +
  facet_wrap(vars(cell), nrow = 2) +
  scale_fill_manual(values = c("#f05c3b", "#8a9097")) +
  labs(fill = 'Estimation for') + rremove("xlab")# + theme(legend.position = "top", legend.direction = "vertical")
dev.off()

pdf("./G4_proportion_density_fold.pdf", height = 4.5, width = 8)
ggplot(mat, aes(x = fold, fill = anno)) + 
geom_histogram(colour = "black",
               lwd = 0.75,
               linetype = 1,
               position = "identity") +
scale_fill_manual(values = c("#f05c3b", "#8a9097")) +
xlab("Proportion of G4s located in activated cCREs") + ylab("Count") +
theme(legend.position = "top", legend.direction = "horizontal") +
rremove("legend.title") + 
facet_grid(cell ~ group, scales='free')
dev.off()






