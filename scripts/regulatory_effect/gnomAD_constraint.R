
set.seed(1)
options(scipen=200)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(bedtoolsr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyr))

# load ccre
ccre <- fread("../data/ccre/V3/G4_cCRE_annotation.txt", sep = "\t", header = TRUE) %>% data.frame() %>% filter(chr %in% paste0("chr", c(1:22, "X")))

# load G4
G4 <- fread("../data/G4/cCRE_G4_annotation.txt", sep = "\t", header = TRUE) %>% data.frame()

gnomad.cst <- fread("../data/gnomad/constraint_z_genome_1kb.qc.download.txt", sep = "\t", header = TRUE) %>% data.frame()

G4.gnomad <- bt.intersect(a = G4 %>% filter(overlap_count <= 1), b = gnomad.cst, wo = TRUE) %>% unique() %>% data.frame()
G4.gnomad <- G4.gnomad %>% group_by(V10) %>% summarise(gnomad.score = sum(V25*V26) / sum(V26)) %>% data.frame()
colnames(G4.gnomad)[1] <- "G4_ID"
G4.gnomad <- G4.gnomad %>% left_join(G4, by = "G4_ID")

G4.zscore <- bind_rows(
    data.frame(zscore = G4.gnomad %>% filter(PLS == 1) %>% select(gnomad.score) %>% unlist(), group = "PLS G4s"), 
    data.frame(zscore = G4.gnomad %>% filter(pELS == 1) %>% select(gnomad.score) %>% unlist(), group = "pELS G4s"),
    data.frame(zscore = G4.gnomad %>% filter(dELS == 1) %>% select(gnomad.score) %>% unlist(), group = "dELS G4s"),
    data.frame(zscore = G4.gnomad %>% filter(CTCF.only == 1) %>% select(gnomad.score) %>% unlist(), group = "CTCF-only G4s"),
    data.frame(zscore = G4.gnomad %>% filter(DNase.H3K4me3 == 1) %>% select(gnomad.score) %>% unlist(), group = "DNase-H3K4me3 G4s"),
    data.frame(zscore = G4.gnomad %>% filter(overlap_count == 0) %>% select(gnomad.score) %>% unlist(), group = "Other G4s"))
G4.zscore$group <- factor(G4.zscore$group, levels = unique(G4.zscore$group))

pdf("../figure/gnomAD/G4_zscores.pdf", width = 3, height = 4)
ggplot(G4.zscore, aes(group, zscore)) + 
  geom_boxplot(aes(fill = group), show.legend = FALSE, width = 0.5) +
  labs(y="gnomAD constraint z-scores") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_brewer(palette = "Pastel1") + rremove("xlab") + rotate_x_text(90)
dev.off()

ccre.gnomad <- bt.intersect(a = ccre, b = gnomad.cst, wo = TRUE) %>% unique() %>% data.frame()
ccre.gnomad <- ccre.gnomad %>% group_by(V5) %>% summarise(gnomad.score = sum(V17*V18) / sum(V18)) %>% data.frame()
colnames(ccre.gnomad)[1] <- "ID"
ccre.gnomad <- ccre.gnomad %>% left_join(ccre, by = "ID")

ccre.zscore <- data.frame(zscore = ccre.gnomad$gnomad.score, cCRE = ccre.gnomad$encodeLabel, group = ccre.gnomad$contain_G4)
ccre.zscore$group <- ifelse(ccre.zscore$group > 0, "G4 cCREs", "Other cCREs")
ccre.zscore$cCRE <- factor(ccre.zscore$cCRE, levels = c("PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3"))

#
wilcox.test(ccre.zscore %>% filter((cCRE == "DNase-H3K4me3") & (group == "G4 cCREs")) %>% select(zscore) %>% unlist(), 
            ccre.zscore %>% filter((cCRE == "DNase-H3K4me3") & (group == "Other cCREs")) %>% select(zscore) %>% unlist())$p.value

pdf("../figure/gnomAD/cCRE_zscores.pdf", width = 5, height = 4)
ggplot(ccre.zscore, aes(cCRE, zscore)) + 
  geom_boxplot(aes(fill = group), width = 0.6) +
  labs(y = "gnomAD constraint z-scores") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_brewer(palette = "Pastel1") + rremove("xlab") + rotate_x_text(90)
dev.off()


