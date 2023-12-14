
set.seed(1)
options(scipen=200)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(bedtoolsr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(paletteer))

# ct used
ct.used <- fread("../data/PANC_CT_used.txt", sep = "\t", header = FALSE) %>% unlist() %>% c()

raw.prom <- fread("../data/pancExp/alternativePromoters_rawPromoterActivity.tsv", sep = "\t", header = FALSE) %>% data.frame()
sample <- fread("../data/pancExp/pcawg_sample_sheet.tsv", sep = "\t", header = TRUE) %>% data.frame()
raw.prom.sample <- raw.prom[1, 3:ncol(raw.prom)] %>% t()
raw.prom.sample <- data.frame(aliquot_id = c(raw.prom.sample)) %>% inner_join(sample, by = "aliquot_id")
raw.prom.sample <- raw.prom.sample[!(str_starts(raw.prom.sample$dcc_specimen_type, "Normal")),]
raw.prom.sample$ct <- str_split(raw.prom.sample$dcc_project_code, "-", simplify = TRUE)[, 1]

pan.atac.files <- list.files("../data/pancATAC/TCGA-ATAC_Cancer_Type-specific_PeakCalls/")
#
pan.atac.ct <- str_split(pan.atac.files, "_",simplify=T)[, 1]
#

raw.prom.sample <- raw.prom.sample %>% filter(ct %in% ct.used)
unlink("../data/pancExp/pcawg_prom_exp.txt")
for (i in 3:ncol(raw.prom)) {
  message(i)
  tmp.id <- raw.prom[1, i]
  tmp.mat <- raw.prom[2:nrow(raw.prom), c(1, 2, i)]
  tmp.mat[, 3] <- as.numeric(tmp.mat[, 3])
  tmp.ct <- raw.prom.sample %>% filter(aliquot_id == tmp.id) %>% select(ct) %>% unlist()
  if (length(tmp.ct) != 0) {
    tmp.mat$ct <- tmp.ct
    fwrite(tmp.mat, "../data/pancExp/pcawg_prom_exp.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
  }
}

# load ccre
ccre <- fread("../data/ccre/V3/G4_cCRE_annotation.txt", sep = "\t", header = TRUE) %>% data.frame() %>% filter(chr %in% paste0("chr", c(1:22, "X")))
ccre.prom <- ccre %>% filter(encodeLabel %in% "PLS")

#fwrite(ccre.prom[, 1:3], "../data/gse_kegg/PLS.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE) 
tss <- fread("../data/ref/tss_regions_200To200.bed", sep = "\t", header = FALSE) %>% data.frame()
tss[, 1] <- paste0("chr", tss[, 1])
tss[, 5] <- str_split(tss[, 5], "\\.", simplify = TRUE)[, 1]
colnames(tss)[5] <- "ENST"

# load AP promoter info
prom.info <- fread("../data/pancExp/1-s2.0-S0092867419309067-mmc1.txt", sep = "\t", header = TRUE) %>% data.frame()
prom.info[, 1] <- str_split(prom.info[, 1], "\\.", simplify = TRUE)[, 1]
colnames(prom.info)[1] <- "ENST"
prom.info <- prom.info %>% inner_join(tss, by = "ENST")
prom.info <- prom.info %>% select(V1:V7, ENST, promoterId)

# 
ccre.prom.tss <- bt.closest(a = bt.sort(ccre.prom), b = bt.sort(prom.info), d = TRUE) %>% unique() %>% data.frame()
ccre.prom.tss <- ccre.prom.tss %>% filter(V17 == 0)
ccre.prom.tss <- ccre.prom.tss %>% select(V16, V8, V5) %>% unique()
colnames(ccre.prom.tss) <- c("promoterId", "annoG4", "ID")
ccre.prom.fltr <- ccre.prom.tss %>% group_by(promoterId) %>% mutate(desc = length(unique(annoG4)) > 1) %>% 
  filter(desc == FALSE) %>% select(promoterId, annoG4, ID) %>% data.frame() %>% unique()


# raw prom activity
raw.prom <- fread("../data/pancExp/pcawg_prom_exp.txt", sep = "\t", header = FALSE) %>% data.frame()
colnames(raw.prom) <- c("promoterId", "ENSG", "exp", "ct")
raw.prom.avg <- raw.prom %>% group_by(promoterId, ct) %>% summarise(avgExp = mean(exp)) %>% data.frame()

all.ts <- unique(raw.prom.avg$ct)
exp.mat <- NULL
act.mat <- NULL
for (tmp.ts in all.ts) {
  message(tmp.ts)
  tmp.raw.avg.exp <- raw.prom.avg %>% filter(ct == tmp.ts)
  tmp.atac <- fread(paste0("../data/pancATAC/TCGA-ATAC_Cancer_Type-specific_PeakCalls/", tmp.ts, "_peakCalls.txt"), sep = "\t", header = TRUE) %>% data.frame()
  tmp.ccre.prom <- bt.intersect(a = ccre.prom, b = tmp.atac, wa = TRUE) %>% data.frame() 
  tmp.ccre.prom.fltr <- ccre.prom.fltr %>% filter(ID %in% unique(tmp.ccre.prom$V5))

  tmp.G4ccre.exp <- tmp.raw.avg.exp[tmp.raw.avg.exp$promoterId %in% tmp.ccre.prom.fltr[tmp.ccre.prom.fltr$annoG4 == 1, "promoterId"], "avgExp"]
  tmp.nonG4ccre.exp <- tmp.raw.avg.exp[tmp.raw.avg.exp$promoterId %in% tmp.ccre.prom.fltr[tmp.ccre.prom.fltr$annoG4 == 0, "promoterId"], "avgExp"]

  exp.mat <- bind_rows(exp.mat, data.frame(exp = tmp.G4ccre.exp, ct = tmp.ts, group = "G4 PLS"),
                       data.frame(exp = tmp.nonG4ccre.exp, ct = tmp.ts, group = "nonG4 PLS"))

  act.mat <- bind_rows(act.mat, 
                       data.frame(p = length(tmp.G4ccre.exp[tmp.G4ccre.exp > 1.5 ]) / length(tmp.G4ccre.exp), ct = tmp.ts, group = "G4 PLS"),
                       data.frame(p = length(tmp.nonG4ccre.exp[tmp.nonG4ccre.exp > 1.5 ]) / length(tmp.nonG4ccre.exp), ct = tmp.ts, group = "nonG4 PLS")
                       )

  wilcox.test(tmp.G4ccre.exp , tmp.nonG4ccre.exp)$p.value %>% print()

}

pdf("../figure/pancancer/PROM_exp.pdf", height = 3, width = 8)
ggplot(exp.mat, aes(x=ct,y=log2(exp+1))) +
geom_boxplot(aes(fill = group), width=.5) +
scale_fill_manual(values = paletteer_d("ggthemes::Classic_Traffic_Light")) +
theme_bw() + rotate_x_text(90) + ylab("Promoter activity (log2)") + rremove("xlab")
dev.off()

ggsave("../figure/pancancer/PROM_active_prop.pdf", 
ggbarplot(act.mat, "ct", "p", color = "group", fill = "group",
          palette = paletteer_d("ggthemes::Classic_Traffic_Light"), position = position_dodge()) + 
rotate_x_text(90) + ylab("Proportion (PROM activity > 1.5)") + rremove("xlab") + rremove("legend.title"), 
width = 8, height = 3)

