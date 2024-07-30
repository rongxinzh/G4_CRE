
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

all.c <- NULL
pan.atac.files <- list.files("../data/pancATAC/TCGA-ATAC_Cancer_Type-specific_PeakCalls/")
for (tmp.atac.file in pan.atac.files) {

  tmp.ct <- str_split(tmp.atac.file, "_", simplify = TRUE)[1]
  
  if (!(tmp.ct %in% ct.used)) {
    next
  }

  message(tmp.ct)
  tmp.atac <- fread(paste0("../data/pancATAC/TCGA-ATAC_Cancer_Type-specific_PeakCalls/", tmp.atac.file), sep = "\t", header = TRUE) %>% data.frame()
  
  ccre.c <- bt.intersect(a = ccre, b = tmp.atac, wa = TRUE, c = TRUE) %>% unique() %>% data.frame()
  ccre.c[, ncol(ccre.c)] <- ifelse(ccre.c[, ncol(ccre.c)] > 0, 1, 0)
  ccre.c <- ccre.c %>% group_by(V7, V8, V9) %>% summarise(count = n()) %>% mutate(freq = count / sum(count)) %>% data.frame()

  all.c <- bind_rows(all.c, data.frame(cCRE = ccre.c[, 1],
                                       Anno = ccre.c[, 2],
                                       Open = ccre.c[, 3],
                                       freq = ccre.c[, 5],
                                       ct = tmp.ct))
}

all.c.open <- all.c %>% filter(Open == 1)
all.c.open$cCRE <- factor(all.c.open$cCRE, levels = c("PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3"))
all.c.open$Anno <- ifelse(all.c.open$Anno == 1, "Contain G4s", "Others")

ggsave("../figure/pancancer/openregion_proportion.pdf", 
  ggbarplot(all.c.open, "ct", "freq", fill = "Anno", color = "Anno", 
    palette = paletteer_d("ggthemes::Classic_Traffic_Light"), position = position_dodge()) + 
  rotate_x_text(90) + ylab("Proportion\n(overlapped with ATAC-seq peaks)") + rremove("xlab") + rremove("legend.title") +
  facet_wrap(vars(cCRE), ncol = 5, scales = "free"), 
  width = 20, height = 3.5)

