
set.seed(1)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(bedtoolsr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggstatsplot))

# load ccre
ccre <- fread("../data/ccre/V3/G4_cCRE_annotation.txt", sep = "\t", header = TRUE) %>% data.frame() %>% filter(chr %in% paste0("chr", c(1:22, "X")))

# load G4
G4 <- fread("../data/G4/cCRE_G4_annotation.txt", sep = "\t", header = TRUE) %>% data.frame()

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

fig <- ggpiestats(data = ccre, x = contain_G4, y = encodeLabel, results.subtitle = FALSE, perc.k = 1) +
	   #theme(text=element_text(family="sans")) + 
	   rremove("legend.title")
ggsave("../figure/enrich/G4_cCRE_prop.pdf", facet(set_palette(fig, c("#ed2200", "#808180")), facet.by = "contain_G4", ncol = 5), width = 10, height = 4)

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
G4.score.grp$group <- factor(G4.score.grp$group, levels = c("PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3", "Others"))
kruskal.test(score ~ group, data = G4.score.grp)$p.value

pairwise.wilcox.test(G4.score.grp$score, G4.score.grp$group)$p.value

fig <- ggboxplot(data = G4.score.grp, x = 'group', y = 'score', fill = 'group', 
								 palette = c("#e64b35", "#4ebad5", "#049f87", "#3c5488", "#f49b7f", "#7d6148"), 
		  		 			 ylab = "|G4Hunter score|", width = .5) + 
			 rremove("xlab") + rremove("legend") + rotate_x_text(90)
ggsave("../figure/enrich/G4Hunter_score_diff.pdf", fig, width = 3, height = 4)

G4.score.grp[, 1] <- factor(G4.score.grp[, 1], levels = unique(G4.score.grp[, 1]))
anova.res <- aov(score ~ group, data = G4.score.grp)
print(summary(anova.res))

tukey.res <- TukeyHSD(anova.res)
print(tukey.res)



