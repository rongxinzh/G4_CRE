
set.seed(1)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(bedtoolsr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggstatsplot))
suppressPackageStartupMessages(library(EnrichedHeatmap))

GrEHmap <- function(data = NULL, center = FALSE, score = 1) {
  data <- data[, 1:3]
  data$score <- score
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

# load ccre
ccre <- fread("../data/ccre/V3/G4_cCRE_annotation.txt", sep = "\t", header = TRUE) %>% data.frame() %>% filter(chr %in% paste0("chr", c(1:22, "X")))

# load G4
G4 <- fread("../data/G4/cCRE_G4_annotation.txt", sep = "\t", header = TRUE) %>% data.frame()
G4.2 <- fread("../data/G4/G4Hunter_w25_s1.2_hg38.txt", sep = "\t", header = FALSE) %>% data.frame() %>% filter(abs(V7) < 1.5)
G4.low <- bind_rows(bt.intersect(a = G4.2 %>% filter(V5 == "+"), b = G4 %>% filter(strand == "+"), wa = TRUE, v = TRUE) %>% unique(),
                    bt.intersect(a = G4.2 %>% filter(V5 == "-"), b = G4 %>% filter(strand == "-"), wa = TRUE, v = TRUE) %>% unique())

G4.low <- G4.low[, 1:3]
G4.stable <- G4[, 1:3]

corr.F.low <- max(sum(G4.low[, 3] - G4.low[, 2] + 1), sum(G4.stable[, 3] - G4.stable[, 2] + 1)) / sum(G4.low[, 3] - G4.low[, 2] + 1)
corr.F.stable <- max(sum(G4.low[, 3] - G4.low[, 2] + 1), sum(G4.stable[, 3] - G4.stable[, 2] + 1)) / sum(G4.stable[, 3] - G4.stable[, 2] + 1)

G4.low.gr <- GrEHmap(G4.low, score = corr.F.low)
G4.stable.gr <- GrEHmap(G4.stable, score = corr.F.stable)

pls.ccre.gr <- GrEHmap(ccre %>% filter(encodeLabel == "PLS"), center = TRUE)
pels.ccre.gr <- GrEHmap(ccre %>% filter(encodeLabel == "pELS"), center = TRUE)
dels.ccre.gr <- GrEHmap(ccre %>% filter(encodeLabel == "dELS"), center = TRUE)
ctcf.ccre.gr <- GrEHmap(ccre %>% filter(encodeLabel == "CTCF-only"), center = TRUE)
dh.ccre.gr <- GrEHmap(ccre %>% filter(encodeLabel == "DNase-H3K4me3"), center = TRUE)

low.mat.pls <- normalizeToMatrix(G4.low.gr, pls.ccre.gr,
                                 value_column = "score", extend = 5000, mean_mode = "coverage", w = 10)
low.mat.pels <- normalizeToMatrix(G4.low.gr, pels.ccre.gr,
                                  value_column = "score", extend = 5000, mean_mode = "coverage", w = 10)
low.mat.dels <- normalizeToMatrix(G4.low.gr, dels.ccre.gr,
                                  value_column = "score", extend = 5000, mean_mode = "coverage", w = 10)
low.mat.ctcf <- normalizeToMatrix(G4.low.gr, ctcf.ccre.gr,
                                  value_column = "score", extend = 5000, mean_mode = "coverage", w = 10)
low.mat.dh <- normalizeToMatrix(G4.low.gr, dh.ccre.gr,
                                value_column = "score", extend = 5000, mean_mode = "coverage", w = 10)

low.mod.mat <- bind_rows(data.frame(density = colMeans(data.frame(low.mat.pls)) * 1000, group = "PLS", window = (-(5000/10):(5000/10))),
                         data.frame(density = colMeans(data.frame(low.mat.pels)) * 1000, group = "pELS", window = (-(5000/10):(5000/10))),
                         data.frame(density = colMeans(data.frame(low.mat.dels)) * 1000, group = "dELS", window = (-(5000/10):(5000/10))),
                         data.frame(density = colMeans(data.frame(low.mat.ctcf)) * 1000, group = "CTCF-only", window = (-(5000/10):(5000/10))),
                         data.frame(density = colMeans(data.frame(low.mat.dh)) * 1000, group = "DNase-H3K4me3", window = (-(5000/10):(5000/10))))
low.mod.mat$group <- factor(low.mod.mat$group, levels = unique(low.mod.mat$group))

stable.mat.pls <- normalizeToMatrix(G4.stable.gr, pls.ccre.gr,
                                 value_column = "score", extend = 5000, mean_mode = "coverage", w = 10)
stable.mat.pels <- normalizeToMatrix(G4.stable.gr, pels.ccre.gr,
                                  value_column = "score", extend = 5000, mean_mode = "coverage", w = 10)
stable.mat.dels <- normalizeToMatrix(G4.stable.gr, dels.ccre.gr,
                                  value_column = "score", extend = 5000, mean_mode = "coverage", w = 10)
stable.mat.ctcf <- normalizeToMatrix(G4.stable.gr, ctcf.ccre.gr,
                                  value_column = "score", extend = 5000, mean_mode = "coverage", w = 10)
stable.mat.dh <- normalizeToMatrix(G4.stable.gr, dh.ccre.gr,
                                value_column = "score", extend = 5000, mean_mode = "coverage", w = 10)

stable.mod.mat <- bind_rows(data.frame(density = colMeans(data.frame(stable.mat.pls)) * 1000, group = "PLS", window = (-(5000/10):(5000/10))),
                         data.frame(density = colMeans(data.frame(stable.mat.pels)) * 1000, group = "pELS", window = (-(5000/10):(5000/10))),
                         data.frame(density = colMeans(data.frame(stable.mat.dels)) * 1000, group = "dELS", window = (-(5000/10):(5000/10))),
                         data.frame(density = colMeans(data.frame(stable.mat.ctcf)) * 1000, group = "CTCF-only", window = (-(5000/10):(5000/10))),
                         data.frame(density = colMeans(data.frame(stable.mat.dh)) * 1000, group = "DNase-H3K4me3", window = (-(5000/10):(5000/10))))
stable.mod.mat$group <- factor(stable.mod.mat$group, levels = unique(stable.mod.mat$group))

all.density <- bind_rows(data.frame(low.mod.mat, subgroup = "Less stable"), 
                         data.frame(stable.mod.mat, subgroup = "Stable"))
all.density$subgroup <- factor(all.density$subgroup, levels = c("Stable", "Less stable"))

fig <- ggplot(all.density, aes(x = window, y = density)) +
       geom_line(aes(color = subgroup)) +
       scale_color_manual(values = c("indianred4", "lightcoral", "darkcyan"), name = "G4 stability") + 
       theme_classic() + facet_wrap(vars(group), scales = "free", ncol = 5)
ggsave("../figure/G4_stability_cCRE_distn.pdf", fig, width = 10, height = 1.6)

save(list = ls(), file = "statistics_G4_stability.RData")

