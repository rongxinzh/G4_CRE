
set.seed(1)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(bedtoolsr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggstatsplot))
suppressPackageStartupMessages(library(EnrichedHeatmap))

# load ccre
ccre <- fread("../data/ccre/V3/G4_cCRE_annotation.txt", sep = "\t", header = TRUE) %>% data.frame() %>% filter(chr %in% paste0("chr", c(1:22, "X")))

# load G4
G4 <- fread("../data/G4/cCRE_G4_annotation.txt", sep = "\t", header = TRUE) %>% data.frame()
G4.rand <- bt.shuffle(G4[, 1:3], "../data/ref/hg38.chrom.sizes", chrom = TRUE, noOverlapping = TRUE, excl = "../data/ref/hg38-blacklist.v2.bed")

# 
ccre.G4.cvg <- bt.coverage(a = ccre, b = bt.merge(bt.sort(G4[, 1:3])))
ccre.randG4.cvg <- bt.coverage(a = ccre, b = bt.merge(bt.sort(G4.rand[, 1:3])))

# proportion
prop.cvg <- data.frame(ccre.G4.cvg[, c(7, 11, 12)])
prop.cvg.rand <- data.frame(ccre.randG4.cvg[, c(7, 11, 12)])
colnames(prop.cvg) <- colnames(prop.cvg.rand) <- c("group", "width", "prop")
prop.cvg.all <- bind_rows(data.frame(prop.cvg, class = "G4"), data.frame(prop.cvg.rand, class = "Shuffled"))
prop.cvg.all$group <- factor(prop.cvg.all$group, level = unique(prop.cvg.all$group))

fig <- ggplot(prop.cvg.all, aes(prop, after_stat(scaled), colour = class)) + 
  		 geom_density() +
  		 geom_vline(data = prop.cvg.all %>% group_by(group, class) %>% summarise(mean = mean(prop)), 
  			 		aes(xintercept = mean, color = class), linetype = "dashed", linewidth = 0.5) +
  	   #theme_classic() +
  	   scale_color_manual(values=c("#ed2200", "#808180")) +
  	   rremove("legend.title") +
  	   xlab("G4 Coverage") +
  	   ylab("Density") +
  	   facet_grid(. ~ group)
ggsave("../figure/enrich/G4_coverage_all.pdf", fig, width = 10, height = 2)

fig <- ggplot(prop.cvg.all %>% filter(prop > 0), aes(prop, after_stat(scaled), colour = class)) + 
  		 geom_density() +
  		 geom_vline(data = prop.cvg.all %>% filter(prop > 0) %>% group_by(group, class) %>% summarise(mean = mean(prop)), 
  			 		aes(xintercept = mean, color = class), linetype = "dashed", linewidth = 0.5) +
  	   #theme_classic() +
  	   scale_color_manual(values=c("#ed2200", "#808180")) +
  	   rremove("legend.title") +
  	   xlab("G4 Coverage") +
  	   ylab("Density") + 
  	   facet_grid(. ~ group)
ggsave("../figure/enrich/G4_coverage_gt0.pdf", fig, width = 10, height = 2)

# G4 width
G4.width <- G4$end - G4$start + 1
G4.width[G4.width < 20] <- 19
G4.width[G4.width > 50] <- 51
pdf("../figure/enrich/G4_width.pdf", width = 6, height = 5)
hist(G4.width, breaks = seq(19, 51, by = 1), col = "steelblue", main = "",
     xlab = "G4 width", ylab = "Frequency", border = "black")
dev.off()

# cCRE width
ccre.width <- ccre
ccre.width$width <- ccre.width[, 3] - ccre.width[, 2] + 1
ccre.width$encodeLabel <- factor(ccre.width$encodeLabel, levels = unique(ccre.width$encodeLabel))
pdf("../figure/enrich/cCRE_width.pdf", width = 15, height = 3)
ggplot(data = ccre.width, aes(x = width, fill = encodeLabel)) + 
geom_histogram(color = "black") + xlab("cCRE width") + rremove("legend.title") + ylab("Frequency") +
scale_fill_manual(values = c("#e64b35", "#4ebad5", "#049f87", "#3c5488", "#f49b7f", "#7d6148")) + 
facet_wrap(~ encodeLabel, scales = 'free_y', ncol = 5)
dev.off()

prop.cvg.all$iscovg <- ifelse(prop.cvg.all$prop > 0, paste0("Contain ", prop.cvg.all$class), "None")
prop.cvg.all$iscovg <- factor(prop.cvg.all$iscovg, level = unique(prop.cvg.all$iscovg))
#prop.cvg.n <- prop.cvg.all %>% group_by(group, class) %>% count(iscovg) %>% data.frame

fig <- ggpiestats(data = prop.cvg.all %>% filter(class == "G4"), x = iscovg, y = group, results.subtitle = FALSE, perc.k = 1) +
	   #theme(text=element_text(family="sans")) + 
	   rremove("legend.title")
ggsave("../figure/enrich/G4_cCRE_prop.pdf", facet(set_palette(fig, c("#ed2200", "#808180")), facet.by = "group", ncol = 5), width = 10, height = 4)

fig <- ggpiestats(data = prop.cvg.all %>% filter(class == "Shuffled"), x = iscovg, y = group, results.subtitle = FALSE, perc.k = 1) + 
	   #theme(text=element_text(family="sans")) + 
	   rremove("legend.title")
ggsave("../figure/enrich/Rand_cCRE_prop.pdf", facet(set_palette(fig, c("#ed2200", "#808180")), facet.by = "group", ncol = 5), width = 10, height = 4)

prop.cvg.gt0 <- prop.cvg.all %>% filter(prop > 0)
t.test(prop.cvg.gt0 %>% filter(group == "DNase-H3K4me3", class == "G4") %>% select(prop) %>% unlist(),
			 prop.cvg.gt0 %>% filter(group == "DNase-H3K4me3", class == "Shuffled") %>% select(prop) %>% unlist())$p.value
#0, 0, 0, 0.002403741, 1.649501e-86
fig <- ggboxplot(data = prop.cvg.gt0, x = 'class', y = 'prop', fill = 'class', palette = c("#ed2200", "#808180"), 
		  		 ylab = "Coverage", width = .5) + facet_wrap(. ~ group, ncol = 5) + 
	   rremove("xlab") + rremove("legend.title")
ggsave("../figure/enrich/G4Hunter_cvggt0_diff.pdf", fig, width = 8, height = 3)

#
G4.nonccre <- G4 %>% filter(overlap_count == 0)
G4.pls <- G4 %>% filter(overlap_count == 1, PLS > 0) 
G4.pels <- G4 %>% filter(overlap_count == 1, pELS > 0) 
G4.dels <- G4 %>% filter(overlap_count == 1, dELS > 0) 
G4.ctcf <- G4 %>% filter(overlap_count == 1, CTCF.only > 0) 
G4.dh <- G4 %>% filter(overlap_count == 1, DNase.H3K4me3 > 0)

G4.score.grp <- bind_rows(
				data.frame(group = "PLS", score = abs(G4.pls[, "max_score"])),
				data.frame(group = "pELS", score = abs(G4.pels[, "max_score"])),
				data.frame(group = "dELS", score = abs(G4.dels[, "max_score"])),
				data.frame(group = "CTCF-only", score = abs(G4.ctcf[, "max_score"])),
				data.frame(group = "DNase-H3K4me3", score = abs(G4.dh[, "max_score"])),
				data.frame(group = "Others", score = abs(G4.nonccre[, "max_score"])),
				)

fig <- ggboxplot(data = G4.score.grp, x = 'group', y = 'score', fill = 'group', 
								 palette = c("#e64b35", "#4ebad5", "#049f87", "#3c5488", "#f49b7f", "#7d6148"), 
		  		 			 ylab = "|G4Hunter score|", width = .5) + 
			 rremove("xlab") + rremove("legend") + rotate_x_text(90)
ggsave("../figure/enrich/G4Hunter_score_diff.pdf", fig, width = 3, height = 4)

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

G4.gr <- GrEHmap(G4)
G4.rand.gr <- GrEHmap(G4.rand)

pls.ccre.gr <- GrEHmap(ccre %>% filter(encodeLabel == "PLS"), center = TRUE)
pels.ccre.gr <- GrEHmap(ccre %>% filter(encodeLabel == "pELS"), center = TRUE)
dels.ccre.gr <- GrEHmap(ccre %>% filter(encodeLabel == "dELS"), center = TRUE)
ctcf.ccre.gr <- GrEHmap(ccre %>% filter(encodeLabel == "CTCF-only"), center = TRUE)
dh.ccre.gr <- GrEHmap(ccre %>% filter(encodeLabel == "DNase-H3K4me3"), center = TRUE)

mat.pls <- normalizeToMatrix(G4.gr, pls.ccre.gr,
                             value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
mat.pels <- normalizeToMatrix(G4.gr, pels.ccre.gr,
                              value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
mat.dels <- normalizeToMatrix(G4.gr, dels.ccre.gr,
                              value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
mat.ctcf <- normalizeToMatrix(G4.gr, ctcf.ccre.gr,
                              value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
mat.dh <- normalizeToMatrix(G4.gr, dh.ccre.gr,
                            value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)

all.mat <- bind_rows(data.frame(density = colMeans(data.frame(mat.pls)) * 1000, group = "PLS", window = (-(1000/10):(1000/10))),
										 data.frame(density = colMeans(data.frame(mat.pels)) * 1000, group = "pELS", window = (-(1000/10):(1000/10))),
										 data.frame(density = colMeans(data.frame(mat.dels)) * 1000, group = "dELS", window = (-(1000/10):(1000/10))),
										 data.frame(density = colMeans(data.frame(mat.ctcf)) * 1000, group = "CTCF-only", window = (-(1000/10):(1000/10))),
										 data.frame(density = colMeans(data.frame(mat.dh)) * 1000, group = "DNase-H3K4me3", window = (-(1000/10):(1000/10)))
					 )
all.mat$group <- factor(all.mat$group, levels = unique(all.mat$group))

fig <- ggplot(all.mat, aes(x = window, y = density)) +
  	 	 geom_line(aes(color = group)) +
  		 scale_color_manual(values = c("#e64b35", "#4ebad5", "#049f87", "#3c5488", "#f49b7f", "#7d6148")) + 
  		 theme_classic() + facet_wrap(vars(group), scales = "free", ncol = 5) + rremove("legend")
ggsave("../figure/enrich/cCRE_G4_distribution.pdf", fig, width = 8, height = 1.6)

mat.pls.rand <- normalizeToMatrix(G4.rand.gr, pls.ccre.gr,
                                  value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
mat.pels.rand <- normalizeToMatrix(G4.rand.gr, pels.ccre.gr,
                                   value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
mat.dels.rand <- normalizeToMatrix(G4.rand.gr, dels.ccre.gr,
                                   value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
mat.ctcf.rand <- normalizeToMatrix(G4.rand.gr, ctcf.ccre.gr,
                                   value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
mat.dh.rand <- normalizeToMatrix(G4.rand.gr, dh.ccre.gr,
                                 value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)

all.mat.rand <- bind_rows(data.frame(density = colMeans(data.frame(mat.pls.rand)) * 1000, group = "PLS", window = (-(1000/10):(1000/10))),
										 			data.frame(density = colMeans(data.frame(mat.pels.rand)) * 1000, group = "pELS", window = (-(1000/10):(1000/10))),
										 			data.frame(density = colMeans(data.frame(mat.dels.rand)) * 1000, group = "dELS", window = (-(1000/10):(1000/10))),
										 			data.frame(density = colMeans(data.frame(mat.ctcf.rand)) * 1000, group = "CTCF-only", window = (-(1000/10):(1000/10))),
										 			data.frame(density = colMeans(data.frame(mat.dh.rand)) * 1000, group = "DNase-H3K4me3", window = (-(1000/10):(1000/10)))
					 )
all.mat.rand$group <- factor(all.mat.rand$group, levels = unique(all.mat.rand$group))

fig <- ggplot(all.mat.rand, aes(x = window, y = density)) +
  	 	 geom_line(aes(color = group)) +
  		 scale_color_manual(values = c("#e64b35", "#4ebad5", "#049f87", "#3c5488", "#f49b7f", "#7d6148")) + 
  		 theme_classic() + facet_wrap(vars(group), scales = "free", ncol = 5) + rremove("legend")
ggsave("../figure/enrich/cCRE_randG4_distribution.pdf", fig, width = 8, height = 1.6)


