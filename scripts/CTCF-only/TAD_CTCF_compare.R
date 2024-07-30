
set.seed(1)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(bedtoolsr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggstatsplot))
suppressPackageStartupMessages(library(EnrichedHeatmap))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(colorRamp2))

# load ccre
ccre <- fread("../data/ccre/V3/G4_cCRE_annotation.txt", sep = "\t", header = TRUE) %>% data.frame() %>% filter(chr %in% paste0("chr", c(1:22, "X")))
ctcf.ccre <- ccre[str_detect(ccre$ccre, "CTCF-bound"),]
ctcfonly.ccre <- ctcf.ccre[str_detect(ctcf.ccre$ccre, "CTCF-only,CTCF-bound"),]
ctcfother.ccre <- ctcf.ccre[!str_detect(ctcf.ccre$ccre, "CTCF-only,CTCF-bound"),]

# load G4
G4 <- fread("../data/G4/cCRE_G4_annotation.txt", sep = "\t", header = TRUE) %>% data.frame()

G4.ctcfonly <- bt.intersect(a = G4, b = ctcfonly.ccre, wa = TRUE, c = TRUE) %>% unique() %>% filter(V17 > 0)
G4.ctcfother <- bt.intersect(a = G4, b = ctcfother.ccre, wa = TRUE, c = TRUE) %>% unique() %>% filter(V17 > 0)

G4.score.grp <- bind_rows(
				data.frame(group = "CTCF-hybrid", score = abs(G4.ctcfother[, "V7"])),
				data.frame(group = "CTCF-only", score = abs(G4.ctcfonly[, "V7"]))
				)

fig <- ggboxplot(data = G4.score.grp, x = 'group', y = 'score', fill = 'group', 
		  		 			 palette = c("#fdaf91", "#9632b8"), ylab = "|G4Hunter score|", width = .5) + 
	   	 rremove("xlab") + rremove("legend")
ggsave("../figure/enrich/G4Hunter_score_CTCF_compare_diff.pdf", fig, width = 3, height = 3)

### TADs
all.f <- list.files("../data/TAD/hg38/")
extend <- c(5000, 10000, 20000, 50000)
tad.res <- NULL
for (tmp.tad in all.f) {
  message(tmp.tad)
  tad.coord <- fread(paste0("../data/TAD/hg38/", tmp.tad), sep = "\t", header = FALSE) %>% data.frame  
  if (!str_detect(tad.coord[1,1], "chr")) {
  	tad.coord[, 1] <- paste0("chr", tad.coord[, 1])
  }
  tad.coord <- tad.coord %>% filter(V1 %in% paste0("chr", c(1:22, "X")))

  for (tmp.extend in extend) {
  	message(tmp.extend)
  	tad.bdy <- bind_rows(data.frame(chr = tad.coord[, 1], start = tad.coord[, 2] - tmp.extend, end = tad.coord[, 2]),
  						 					 data.frame(chr = tad.coord[, 1], start = tad.coord[, 3], end = tad.coord[, 3] + tmp.extend))
  	tad.bdy[, 2] <- ifelse(tad.bdy[, 2] < 0, 1, tad.bdy[, 2])
  	tad.bdy <- tad.bdy %>% filter(start < end)

  	tmp.oly.int <- bt.intersect(a = ctcfonly.ccre, b = tad.bdy, wa = TRUE) %>% unique()
  	oly.ccre.c <- nrow(tmp.oly.int)

  	otc.ccre.c <- NULL
  	for (i in 1:1000) {
  	  tmp.otr.ccre <- ctcfother.ccre[sample(1:dim(ctcfother.ccre)[1], dim(ctcfonly.ccre)[1]), ]
   	  tmp.otr.int <- bt.intersect(a = tmp.otr.ccre, b = tad.bdy, wa = TRUE) %>% unique()
  	  otc.ccre.c <- c(otc.ccre.c, nrow(tmp.otr.int))
  	}
  	tad.res <- bind_rows(tad.res, data.frame(file = tmp.tad, gt = (oly.ccre.c - mean(otc.ccre.c)) / sd(otc.ccre.c), extend = tmp.extend))
  }
}

save(list = ls(), file = "./CTCF_compare.RData")

tad.mat <- tad.res
tad.mat <- spread(tad.mat, key = "file", value = "gt")
rownames(tad.mat) <- tad.mat[, 1]
tad.mat <- tad.mat[, -1]
colnames(tad.mat) <- str_split(colnames(tad.mat), "(\\.txt|\\.domains)", simplify = TRUE)[, 1]

pdf(paste0("../figure/enrich/TAD_zscore.pdf"), width = 12, height = 8)
ht <- Heatmap(tad.mat, cluster_rows = FALSE, cluster_columns = TRUE, name = "Z-score",
							col = colorRamp2(c(5, 0, -10), c("firebrick3", "white", "navy")),
							column_names_max_height = unit(10, "cm"),
              show_column_names = TRUE, rect_gp = gpar(col = "white", lwd = 2))
draw(ht)
dev.off()
