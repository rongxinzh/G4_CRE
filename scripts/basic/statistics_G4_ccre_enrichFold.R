
set.seed(1)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(bedtoolsr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(gridExtra))

chr.use <- paste0("chr", c(1:22, "X"))
# load ccre
ccre <- fread("../data/ccre/V3/G4_cCRE_annotation.txt", sep = "\t", header = TRUE) %>% data.frame() %>% filter(chr %in% chr.use)

# load G4
G4 <- fread("../data/G4/cCRE_G4_annotation.txt", sep = "\t", header = TRUE) %>% data.frame()

ccre.obs <- bt.nuc(fi = "../data/ref/hg38.fa", bed = ccre)
ccre.obs <- bt.coverage(a = ccre.obs, b = G4)

obs.cvg <- ccre.obs %>% group_by(V7) %>% summarise(G4_covered = sum(V19), cCRE_len = sum(V20), coverage = G4_covered / cCRE_len) %>% data.frame()
obs.dens <- ccre.obs %>% group_by(V7) %>% summarise(G4_count = sum(V18), cCRE_len = sum(V20), density = G4_count / cCRE_len) %>% data.frame()
obs.gc <- ccre.obs %>% group_by(V7) %>% summarise(GC_covered = sum(V12) + sum(V13), cCRE_len = sum(V17), ratio = GC_covered / cCRE_len) %>% data.frame()

g.use <- bt.subtract(a = "../data/ref/hg38.chromused.bed", b = "../data/ref/gap_hg38.txt")
g.use.info <- bt.nuc(fi = "../data/ref/hg38.fa", bed = g.use)
g.use.info <- bt.coverage(a = g.use.info, b = G4)

g.use.info.cvg <- sum(g.use.info[, "V14"]) / sum(g.use.info[, "V15"])
g.use.info.dens <- sum(g.use.info[, "V13"]) / sum(g.use.info[, "V15"])
g.use.info.gc <- sum(g.use.info[, "V7"] + g.use.info[, "V8"]) / sum(g.use.info[, "V12"])

fold.data <- bind_rows(data.frame(group = obs.cvg[, 1], value = obs.cvg[, 4] / g.use.info.cvg, label = "Without GC% correction", subgroup = "Coverage"),
                       data.frame(group = obs.cvg[, 1], value = (obs.cvg[, 4] / obs.gc[, 4]) / (g.use.info.cvg / g.use.info.gc), label = "With GC% correction", subgroup = "Coverage"),
                       data.frame(group = obs.cvg[, 1], value = obs.dens[, 4] / g.use.info.dens, label = "Without GC% correction", subgroup = "Density"),
                       data.frame(group = obs.cvg[, 1], value = (obs.dens[, 4] / obs.gc[, 4]) / (g.use.info.dens / g.use.info.gc), label = "With GC% correction", subgroup = "Density"))
fold.data$group <- factor(fold.data$group, levels = c("PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3"))
fold.data$label <- factor(fold.data$label, levels = c("Without GC% correction", "With GC% correction"))

p1 <- ggplot(fold.data %>% filter(subgroup == "Coverage"), aes(x = group, y = value, fill = label)) +
  geom_bar(stat = "identity", position = "dodge") + ylim(0, 7.5) + 
  scale_fill_manual(values = c("Without GC% correction" = "#ce1a21", "With GC% correction" = "#e89699")) +
  labs(x = "",
       y = "Fold enrichment",
       fill = "G4 coverage enrichment") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

p2 <- ggplot(fold.data %>% filter(subgroup == "Density"), aes(x = group, y = value, fill = label)) +
  geom_bar(stat = "identity", position = "dodge") + ylim(0, 7.5) + 
  scale_fill_manual(values = c("Without GC% correction" = "#3a5898", "With GC% correction" = "#afbed4")) +
  labs(x = "",
       y = "Fold enrichment",
       fill = "G4 density enrichment") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

combine.p <- grid.arrange(p1, p2, ncol = 2)
ggsave("../figure/GC_correction.pdf", combine.p, width = 9, height = 3.5)

