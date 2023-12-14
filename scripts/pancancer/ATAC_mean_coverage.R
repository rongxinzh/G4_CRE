
set.seed(1)
options(scipen=200)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(bedtoolsr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(EnrichedHeatmap))
suppressPackageStartupMessages(library(paletteer))

# ct used
ct.used <- fread("../data/PANC_CT_used.txt", sep = "\t", header = FALSE) %>% unlist() %>% c()

# load ccre
ccre <- fread("../data/ccre/V3/G4_cCRE_annotation.txt", sep = "\t", header = TRUE) %>% data.frame() %>% filter(chr %in% paste0("chr", c(1:22, "X")))

all.cvg <- NULL
pan.atac.files <- list.files("../data/pancATAC/TCGA-ATAC_Cancer_Type-specific_PeakCalls/")
for (tmp.atac.file in pan.atac.files) {

  tmp.ct <- str_split(tmp.atac.file, "_", simplify = TRUE)[1]
  
  if (!(tmp.ct %in% ct.used)) {
    next
  }
  
  message(tmp.ct)
  tmp.atac <- fread(paste0("../data/pancATAC/TCGA-ATAC_Cancer_Type-specific_PeakCalls/", tmp.atac.file), sep = "\t", header = TRUE) %>% data.frame()
  
  ccre.cvg <- bt.coverage(a = ccre, b = tmp.atac) %>% unique() %>% data.frame()

  message(wilcox.test(ccre.cvg[(ccre.cvg[, 8] == 1) & (ccre.cvg[, 7] == "PLS"), 12], ccre.cvg[(ccre.cvg[, 8] == 0) & (ccre.cvg[, 7] == "PLS"), 12], alternative = "greater")$p.value)
  message(wilcox.test(ccre.cvg[(ccre.cvg[, 8] == 1) & (ccre.cvg[, 7] == "pELS"), 12], ccre.cvg[(ccre.cvg[, 8] == 0) & (ccre.cvg[, 7] == "pELS"), 12], alternative = "greater")$p.value)
  message(wilcox.test(ccre.cvg[(ccre.cvg[, 8] == 1) & (ccre.cvg[, 7] == "dELS"), 12], ccre.cvg[(ccre.cvg[, 8] == 0) & (ccre.cvg[, 7] == "dELS"), 12], alternative = "greater")$p.value)
  message(wilcox.test(ccre.cvg[(ccre.cvg[, 8] == 1) & (ccre.cvg[, 7] == "CTCF-only"), 12], ccre.cvg[(ccre.cvg[, 8] == 0) & (ccre.cvg[, 7] == "CTCF-only"), 12], alternative = "greater")$p.value)
  message(wilcox.test(ccre.cvg[(ccre.cvg[, 8] == 1) & (ccre.cvg[, 7] == "DNase-H3K4me3"), 12], ccre.cvg[(ccre.cvg[, 8] == 0) & (ccre.cvg[, 7] == "DNase-H3K4me3"), 12], alternative = "greater")$p.value)

  all.cvg <- bind_rows(all.cvg, data.frame(cCRE = ccre.cvg[, 7],
                                           Anno = ccre.cvg[, 8],
                                           coverage = ccre.cvg[, 12],
                                           ct = tmp.ct))
}

all.cvg$cCRE <- factor(all.cvg$cCRE, levels = c("PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3"))
all.cvg$Anno <- ifelse(all.cvg$Anno == 1, "Contain G4s", "Others")
all.cvg$Anno <- factor(all.cvg$Anno, levels = c("Contain G4s", "Others"))
all.cvg$ct <- factor(all.cvg$ct, levels = unique(all.cvg$ct))

all.cvg <- all.cvg %>% group_by(cCRE, Anno, ct) %>% summarise(mean = mean(coverage)) %>% data.frame()

pdf("../figure/pancancer/mean_coverage.pdf", width = 6, height = 6)
ggplot(data = all.cvg, aes(y = mean, x = Anno, fill = Anno, color = cCRE)) + 
    geom_line(aes(group = ct), col = "black") +
    geom_point(aes(shape = ct), size = 1.5, 
    position = position_dodge(width = 0.2), show.legend = TRUE) +
    facet_wrap(~ cCRE, nrow = 1, strip.position = "top") +
    scale_color_manual(values=paletteer_d("ggthemes::calc")) +
    scale_shape_manual(values = 0:22) +
    theme_classic() +
    theme(strip.background = element_blank(), strip.text.x = element_blank()) +
    rremove("xlab") + ylab("Mean coverage") + rotate_x_text(90)
dev.off()

