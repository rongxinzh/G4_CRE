
set.seed(1)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(bedtoolsr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(colorRamp2))
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

# load ccre
ccre <- fread("../data/ccre/V3/G4_cCRE_annotation.txt", sep = "\t", header = TRUE) %>% data.frame() %>% filter(chr %in% paste0("chr", c(1:22, "X")))
ccre$ctcf <- ifelse(str_detect(ccre$ccre, "CTCF-bound"), ifelse(str_detect(ccre$ccre, "CTCF-only,CTCF-bound"), 1, 2), 0)

ccre.02 <- ccre %>% filter(ctcf != 1)

m1 <- as.matrix(ccre.02 %>% filter(ctcf == 2) %>% select(encodeLabel, contain_G4) %>% table())
m1 <- data.frame(ccre = rownames(m1), prop = as.numeric(m1[, 2] / (m1[, 1] + m1[, 2])), group = "CTCF-bound")
m2 <- as.matrix(ccre.02 %>% filter(ctcf == 0) %>% select(encodeLabel, contain_G4) %>% table())
m2 <- data.frame(ccre = rownames(m2), prop = as.numeric(m2[, 2] / (m2[, 1] + m2[, 2])), group = "CTCF-free")
m3 <- bind_rows(m1, m2)
m3$ccre <- factor(m3$ccre, levels = c("PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3"))
m3$group <- factor(m3$group, levels = c("CTCF-bound", "CTCF-free"))

fig <- ggbarplot(m3, "ccre", "prop", fill = "group", color = "group", palette = c("#0F6B99", "#A20056"), position = position_dodge(0.9)) + 
	   rremove("xlab") + rremove("legend.title") + ylab("Proportion (contain G4s)")
ggsave(paste0("../figure/enrich/CTCF_bound_prop.pdf"), fig, width = 5, height = 3.8)

G4 <- fread("../data/G4/cCRE_G4_annotation.txt", sep = "\t", header = TRUE) %>% data.frame()
G4.gr <- GrEHmap(G4)

pls.ctcf.gr <- GrEHmap(ccre %>% filter(ccre == "PLS,CTCF-bound"), center = TRUE)
pls.gr <- GrEHmap(ccre %>% filter(ccre == "PLS"), center = TRUE)
pels.ctcf.gr <- GrEHmap(ccre %>% filter(ccre == "pELS,CTCF-bound"), center = TRUE)
pels.gr <- GrEHmap(ccre %>% filter(ccre == "pELS"), center = TRUE)
dels.ctcf.gr <- GrEHmap(ccre %>% filter(ccre == "dELS,CTCF-bound"), center = TRUE)
dels.gr <- GrEHmap(ccre %>% filter(ccre == "dELS"), center = TRUE)
dh.ctcf.gr <- GrEHmap(ccre %>% filter(ccre == "DNase-H3K4me3,CTCF-bound"), center = TRUE)
dh.gr <- GrEHmap(ccre %>% filter(ccre == "DNase-H3K4me3"), center = TRUE)
ctcfonly.gr <- GrEHmap(ccre %>% filter(ccre == "CTCF-only,CTCF-bound"), center = TRUE)

mat.pls.ctcf <- normalizeToMatrix(G4.gr, pls.ctcf.gr,
                             		  value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
mat.pls <- normalizeToMatrix(G4.gr, pls.gr,
                             value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
mat.pels.ctcf <- normalizeToMatrix(G4.gr, pels.ctcf.gr,
                             		   value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
mat.pels <- normalizeToMatrix(G4.gr, pels.gr,
                             	value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
mat.dels.ctcf <- normalizeToMatrix(G4.gr, dels.ctcf.gr,
                             			 value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
mat.dels <- normalizeToMatrix(G4.gr, dels.gr,
                             	value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
mat.dh.ctcf <- normalizeToMatrix(G4.gr, dh.ctcf.gr,
                             		 value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
mat.dh <- normalizeToMatrix(G4.gr, dh.gr,
                            value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
mat.ctcfonly <- normalizeToMatrix(G4.gr, ctcfonly.gr,
                             		  value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)

mat <- bind_rows(data.frame(density = colMeans(data.frame(mat.pls.ctcf)), ctcf = "CTCF-bound", group = "PLS,CTCF-bound", subgroup = "PLS", window = (-(1000/10):(1000/10))),
					  data.frame(density = colMeans(data.frame(mat.pls)), ctcf = "CTCF-free", group = "PLS", subgroup = "PLS", window = (-(1000/10):(1000/10))),
					  data.frame(density = colMeans(data.frame(mat.pels.ctcf)), ctcf = "CTCF-bound", group = "pELS,CTCF-bound", subgroup = "pELS", window = (-(1000/10):(1000/10))),
					  data.frame(density = colMeans(data.frame(mat.pels)), ctcf = "CTCF-free", group = "pELS", subgroup = "pELS", window = (-(1000/10):(1000/10))),
					  data.frame(density = colMeans(data.frame(mat.dels.ctcf)), ctcf = "CTCF-bound", group = "dELS,CTCF-bound", subgroup = "dELS", window = (-(1000/10):(1000/10))),
					  data.frame(density = colMeans(data.frame(mat.dels)), ctcf = "CTCF-free", group = "dELS", subgroup = "dELS", window = (-(1000/10):(1000/10))),
					  data.frame(density = colMeans(data.frame(mat.dh.ctcf)), ctcf = "CTCF-bound", group = "DNase-H3K4me3,CTCF-bound", subgroup = "DNase-H3K4me3", window = (-(1000/10):(1000/10))),
					  data.frame(density = colMeans(data.frame(mat.dh)), ctcf = "CTCF-free", group = "DNase-H3K4me3", subgroup = "DNase-H3K4me3", window = (-(1000/10):(1000/10))),
					  data.frame(density = colMeans(data.frame(mat.ctcfonly)), ctcf = "CTCF-only", group = "CTCF-only", subgroup = "PLS", window = (-(1000/10):(1000/10))),
					  data.frame(density = colMeans(data.frame(mat.ctcfonly)), ctcf = "CTCF-only", group = "CTCF-only", subgroup = "pELS", window = (-(1000/10):(1000/10))),
					  data.frame(density = colMeans(data.frame(mat.ctcfonly)), ctcf = "CTCF-only", group = "CTCF-only", subgroup = "dELS", window = (-(1000/10):(1000/10))),
					  data.frame(density = colMeans(data.frame(mat.ctcfonly)), ctcf = "CTCF-only", group = "CTCF-only", subgroup = "DNase-H3K4me3", window = (-(1000/10):(1000/10)))
					 )
mat$group <- factor(mat$group, levels = unique(mat$group))
mat$subgroup <- factor(mat$subgroup, levels = unique(mat$subgroup))
mat$ctcf <- factor(mat$ctcf, levels = unique(mat$ctcf))

fig <- ggplot(mat, aes(x = window, y = density, group = group)) +
  	   geom_line(aes(color = ctcf), linewidth = 0.8) +
  	   scale_color_manual(values = c("#F05C3B", "#197EC0", "#D2AF81")) + 
  	   theme_classic() + facet_wrap(vars(subgroup), scales = "free", ncol = 4) + 
  	   rremove("xlab") + rremove("legend.title") + ylab("Density (pG4)") 
ggsave("../figure/enrich/CTCF_bound_density.pdf", fig, width = 12, height = 2.2)

# cell line analysis
# K562 and HepG2
k562.G4 <- fread("../data/G4/K562_G4_hg38.bed", sep = "\t", header = FALSE) %>% 
		   data.frame() %>% filter(V1 %in% paste0("chr", c(1:22, "X")))
hepg2.G4 <- fread("../data/G4/HepG2_G4_hg38.bed", sep = "\t", header = FALSE) %>% 
			data.frame() %>% filter(V1 %in% paste0("chr", c(1:22, "X")))

k562.ccre <- fread("../data/ccre/V3/ENCFF414OGC_ENCFF806YEZ_ENCFF849TDM_ENCFF736UDR.7group.bed", sep = "\t", header = FALSE) %>% 
			 data.frame() %>% filter(V1 %in% paste0("chr", c(1:22, "X")))
hepg2.ccre <- fread("../data/ccre/V3/ENCFF546MZK_ENCFF732PJK_ENCFF795ONN_ENCFF357NFO.7group.bed", sep = "\t", header = FALSE) %>% 
			  data.frame() %>% filter(V1 %in% paste0("chr", c(1:22, "X")))

k562.G4.gr <- GrEHmap(k562.G4)
hepg2.G4.gr <- GrEHmap(hepg2.G4)

k562.pls.ctcf.gr <- GrEHmap(k562.ccre %>% filter(V10 == "PLS,CTCF-bound"), center = TRUE)
k562.pls.gr <- GrEHmap(k562.ccre %>% filter(V10 == "PLS"), center = TRUE)
k562.pels.ctcf.gr <- GrEHmap(k562.ccre %>% filter(V10 == "pELS,CTCF-bound"), center = TRUE)
k562.pels.gr <- GrEHmap(k562.ccre %>% filter(V10 == "pELS"), center = TRUE)
k562.dels.ctcf.gr <- GrEHmap(k562.ccre %>% filter(V10 == "dELS,CTCF-bound"), center = TRUE)
k562.dels.gr <- GrEHmap(k562.ccre %>% filter(V10 == "dELS"), center = TRUE)
k562.dh.ctcf.gr <- GrEHmap(k562.ccre %>% filter(V10 == "DNase-H3K4me3,CTCF-bound"), center = TRUE)
k562.dh.gr <- GrEHmap(k562.ccre %>% filter(V10 == "DNase-H3K4me3"), center = TRUE)
k562.ctcfonly.gr <- GrEHmap(k562.ccre %>% filter(V10 == "CTCF-only,CTCF-bound"), center = TRUE)

hepg2.pls.ctcf.gr <- GrEHmap(hepg2.ccre %>% filter(V10 == "PLS,CTCF-bound"), center = TRUE)
hepg2.pls.gr <- GrEHmap(hepg2.ccre %>% filter(V10 == "PLS"), center = TRUE)
hepg2.pels.ctcf.gr <- GrEHmap(hepg2.ccre %>% filter(V10 == "pELS,CTCF-bound"), center = TRUE)
hepg2.pels.gr <- GrEHmap(hepg2.ccre %>% filter(V10 == "pELS"), center = TRUE)
hepg2.dels.ctcf.gr <- GrEHmap(hepg2.ccre %>% filter(V10 == "dELS,CTCF-bound"), center = TRUE)
hepg2.dels.gr <- GrEHmap(hepg2.ccre %>% filter(V10 == "dELS"), center = TRUE)
hepg2.dh.ctcf.gr <- GrEHmap(hepg2.ccre %>% filter(V10 == "DNase-H3K4me3,CTCF-bound"), center = TRUE)
hepg2.dh.gr <- GrEHmap(hepg2.ccre %>% filter(V10 == "DNase-H3K4me3"), center = TRUE)
hepg2.ctcfonly.gr <- GrEHmap(hepg2.ccre %>% filter(V10 == "CTCF-only,CTCF-bound"), center = TRUE)

mat.k562.pls.ctcf <- normalizeToMatrix(k562.G4.gr, k562.pls.ctcf.gr,
                             		   value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
mat.k562.pls <- normalizeToMatrix(k562.G4.gr, k562.pls.gr,
                             	  value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
mat.k562.pels.ctcf <- normalizeToMatrix(k562.G4.gr, k562.pels.ctcf.gr,
                             			value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
mat.k562.pels <- normalizeToMatrix(k562.G4.gr, k562.pels.gr,
                             	   value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
mat.k562.dels.ctcf <- normalizeToMatrix(k562.G4.gr, k562.dels.ctcf.gr,
                             			value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
mat.k562.dels <- normalizeToMatrix(k562.G4.gr, k562.dels.gr,
                             	   value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
mat.k562.dh.ctcf <- normalizeToMatrix(k562.G4.gr, k562.dh.ctcf.gr,
                             		  value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
mat.k562.dh <- normalizeToMatrix(k562.G4.gr, k562.dh.gr,
                             	 value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
mat.k562.ctcfonly <- normalizeToMatrix(k562.G4.gr, k562.ctcfonly.gr,
                             		   value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)

mat.hepg2.pls.ctcf <- normalizeToMatrix(hepg2.G4.gr, hepg2.pls.ctcf.gr,
                             		   			value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
mat.hepg2.pls <- normalizeToMatrix(hepg2.G4.gr, hepg2.pls.gr,
                             	  	 value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
mat.hepg2.pels.ctcf <- normalizeToMatrix(hepg2.G4.gr, hepg2.pels.ctcf.gr,
                             						 value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
mat.hepg2.pels <- normalizeToMatrix(hepg2.G4.gr, hepg2.pels.gr,
                             	   		value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
mat.hepg2.dels.ctcf <- normalizeToMatrix(hepg2.G4.gr, hepg2.dels.ctcf.gr,
                             						 value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
mat.hepg2.dels <- normalizeToMatrix(hepg2.G4.gr, hepg2.dels.gr,
                             	   	  value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
mat.hepg2.dh.ctcf <- normalizeToMatrix(hepg2.G4.gr, hepg2.dh.ctcf.gr,
                             		       value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
mat.hepg2.dh <- normalizeToMatrix(hepg2.G4.gr, hepg2.dh.gr,
                             	 	  value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
mat.hepg2.ctcfonly <- normalizeToMatrix(hepg2.G4.gr, hepg2.ctcfonly.gr,
                             		   			value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)

k562.mat <- bind_rows(data.frame(density = colMeans(data.frame(mat.k562.pls.ctcf)), ctcf = "CTCF-bound", group = "PLS,CTCF-bound", subgroup = "PLS", window = (-(1000/10):(1000/10))),
					  data.frame(density = colMeans(data.frame(mat.k562.pls)), ctcf = "CTCF-free", group = "PLS", subgroup = "PLS", window = (-(1000/10):(1000/10))),
					  data.frame(density = colMeans(data.frame(mat.k562.pels.ctcf)), ctcf = "CTCF-bound", group = "pELS,CTCF-bound", subgroup = "pELS", window = (-(1000/10):(1000/10))),
					  data.frame(density = colMeans(data.frame(mat.k562.pels)), ctcf = "CTCF-free", group = "pELS", subgroup = "pELS", window = (-(1000/10):(1000/10))),
					  data.frame(density = colMeans(data.frame(mat.k562.dels.ctcf)), ctcf = "CTCF-bound", group = "dELS,CTCF-bound", subgroup = "dELS", window = (-(1000/10):(1000/10))),
					  data.frame(density = colMeans(data.frame(mat.k562.dels)), ctcf = "CTCF-free", group = "dELS", subgroup = "dELS", window = (-(1000/10):(1000/10))),
					  data.frame(density = colMeans(data.frame(mat.k562.dh.ctcf)), ctcf = "CTCF-bound", group = "DNase-H3K4me3,CTCF-bound", subgroup = "DNase-H3K4me3", window = (-(1000/10):(1000/10))),
					  data.frame(density = colMeans(data.frame(mat.k562.dh)), ctcf = "CTCF-free", group = "DNase-H3K4me3", subgroup = "DNase-H3K4me3", window = (-(1000/10):(1000/10))),
					  data.frame(density = colMeans(data.frame(mat.k562.ctcfonly)), ctcf = "CTCF-only", group = "CTCF-only", subgroup = "PLS", window = (-(1000/10):(1000/10))),
					  data.frame(density = colMeans(data.frame(mat.k562.ctcfonly)), ctcf = "CTCF-only", group = "CTCF-only", subgroup = "pELS", window = (-(1000/10):(1000/10))),
					  data.frame(density = colMeans(data.frame(mat.k562.ctcfonly)), ctcf = "CTCF-only", group = "CTCF-only", subgroup = "dELS", window = (-(1000/10):(1000/10))),
					  data.frame(density = colMeans(data.frame(mat.k562.ctcfonly)), ctcf = "CTCF-only", group = "CTCF-only", subgroup = "DNase-H3K4me3", window = (-(1000/10):(1000/10)))
					 )
k562.mat$group <- factor(k562.mat$group, levels = unique(k562.mat$group))
k562.mat$subgroup <- factor(k562.mat$subgroup, levels = unique(k562.mat$subgroup))
k562.mat$ctcf <- factor(k562.mat$ctcf, levels = unique(k562.mat$ctcf))

fig <- ggplot(k562.mat, aes(x = window, y = density, group = group)) +
  	   geom_line(aes(color = ctcf), linewidth = 0.8) +
  	   scale_color_manual(values = c("#F05C3B", "#197EC0", "#D2AF81")) + 
  	   theme_classic() + facet_wrap(vars(subgroup), scales = "free", ncol = 4) + 
  	   rremove("xlab") + rremove("legend.title") + ylab("Density (K562 G4 ChIP-seq)") 
ggsave("../figure/enrich/CTCF_bound_density_K562.pdf", fig, width = 12, height = 2.2)

hepg2.mat <- bind_rows(data.frame(density = colMeans(data.frame(mat.hepg2.pls.ctcf)), ctcf = "CTCF-bound", group = "PLS,CTCF-bound", subgroup = "PLS", window = (-(1000/10):(1000/10))),
					  data.frame(density = colMeans(data.frame(mat.hepg2.pls)), ctcf = "CTCF-free", group = "PLS", subgroup = "PLS", window = (-(1000/10):(1000/10))),
					  data.frame(density = colMeans(data.frame(mat.hepg2.pels.ctcf)), ctcf = "CTCF-bound", group = "pELS,CTCF-bound", subgroup = "pELS", window = (-(1000/10):(1000/10))),
					  data.frame(density = colMeans(data.frame(mat.hepg2.pels)), ctcf = "CTCF-free", group = "pELS", subgroup = "pELS", window = (-(1000/10):(1000/10))),
					  data.frame(density = colMeans(data.frame(mat.hepg2.dels.ctcf)), ctcf = "CTCF-bound", group = "dELS,CTCF-bound", subgroup = "dELS", window = (-(1000/10):(1000/10))),
					  data.frame(density = colMeans(data.frame(mat.hepg2.dels)), ctcf = "CTCF-free", group = "dELS", subgroup = "dELS", window = (-(1000/10):(1000/10))),
					  data.frame(density = colMeans(data.frame(mat.hepg2.dh.ctcf)), ctcf = "CTCF-bound", group = "DNase-H3K4me3,CTCF-bound", subgroup = "DNase-H3K4me3", window = (-(1000/10):(1000/10))),
					  data.frame(density = colMeans(data.frame(mat.hepg2.dh)), ctcf = "CTCF-free", group = "DNase-H3K4me3", subgroup = "DNase-H3K4me3", window = (-(1000/10):(1000/10))),
					  data.frame(density = colMeans(data.frame(mat.hepg2.ctcfonly)), ctcf = "CTCF-only", group = "CTCF-only", subgroup = "PLS", window = (-(1000/10):(1000/10))),
					  data.frame(density = colMeans(data.frame(mat.hepg2.ctcfonly)), ctcf = "CTCF-only", group = "CTCF-only", subgroup = "pELS", window = (-(1000/10):(1000/10))),
					  data.frame(density = colMeans(data.frame(mat.hepg2.ctcfonly)), ctcf = "CTCF-only", group = "CTCF-only", subgroup = "dELS", window = (-(1000/10):(1000/10))),
					  data.frame(density = colMeans(data.frame(mat.hepg2.ctcfonly)), ctcf = "CTCF-only", group = "CTCF-only", subgroup = "DNase-H3K4me3", window = (-(1000/10):(1000/10)))
					 )
hepg2.mat$group <- factor(hepg2.mat$group, levels = unique(hepg2.mat$group))
hepg2.mat$subgroup <- factor(hepg2.mat$subgroup, levels = unique(hepg2.mat$subgroup))
hepg2.mat$ctcf <- factor(hepg2.mat$ctcf, levels = unique(hepg2.mat$ctcf))

fig <- ggplot(hepg2.mat, aes(x = window, y = density, group = group)) +
  	   geom_line(aes(color = ctcf), linewidth = 0.8) +
  	   scale_color_manual(values = c(c("#F05C3B", "#197EC0", "#D2AF81"))) + 
  	   theme_classic() + facet_wrap(vars(subgroup), scales = "free", ncol = 4) + 
  	   rremove("xlab") + rremove("legend.title") + ylab("Density (HepG2 G4 ChIP-seq)") 
ggsave("../figure/enrich/CTCF_bound_density_HepG2.pdf", fig, width = 12, height = 2.2)



p1 <- ggplot(data.frame(density = colMeans(data.frame(mat.k562.ctcfonly)), window = (-(1000/10):(1000/10))), 
						 aes(x = window, y = density)) +
  	  geom_line(linewidth = 0.8) +
  	  theme_classic() + ylab("Density (K562 G4 ChIP-seq)") +
  	  rremove("xlab") + rremove("legend")
ggsave("../figure/enrich/CTCF_only_density_K562.pdf", p1, width = 4.5, height = 3)

p2 <- ggplot(data.frame(density = colMeans(data.frame(mat.hepg2.ctcfonly)), window = (-(1000/10):(1000/10))), 
						 aes(x = window, y = density)) +
  	  geom_line(linewidth = 0.8) +
  	  theme_classic() + ylab("Density (HepG2 G4 ChIP-seq)") +
  	  rremove("xlab") + rremove("legend")
ggsave("../figure/enrich/CTCF_only_density_HepG2.pdf", p2, width = 4.5, height = 3)
