
set.seed(1)
options(scipen=200)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(bedtoolsr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(random))

# load ccre
pls.ccre <- fread("../data/ccre/V3/GRCh38-cCREs.PLS.bed", sep = "\t", header = FALSE) %>% data.frame() %>% filter(V1 %in% paste0("chr", c(1:22, "X")))
pels.ccre <- fread("../data/ccre/V3/GRCh38-cCREs.pELS.bed", sep = "\t", header = FALSE) %>% data.frame() %>% filter(V1 %in% paste0("chr", c(1:22, "X")))
dels.ccre <- fread("../data/ccre/V3/GRCh38-cCREs.dELS.bed", sep = "\t", header = FALSE) %>% data.frame() %>% filter(V1 %in% paste0("chr", c(1:22, "X")))
ctcf.ccre <- fread("../data/ccre/V3/GRCh38-cCREs.CTCF-only.bed", sep = "\t", header = FALSE) %>% data.frame() %>% filter(V1 %in% paste0("chr", c(1:22, "X")))
dh.ccre <- fread("../data/ccre/V3/GRCh38-cCREs.DNase-H3K4me3.bed", sep = "\t", header = FALSE) %>% data.frame() %>% filter(V1 %in% paste0("chr", c(1:22, "X")))

all.ccre <- bind_rows(data.frame(pls.ccre, group = "PLS"),
                      data.frame(pels.ccre, group = "pELS"),
                      data.frame(dels.ccre, group = "dELS"),
                      data.frame(ctcf.ccre, group = "CTCF-only"),
                      data.frame(dh.ccre, group = "DNase-H3K4me3"))

# load G4
G4 <- fread("../data/G4/G4Hunter_w25_s1.5_hg38.txt", sep = "\t", header = FALSE) %>% data.frame()

all.ccre.G4 <- bt.intersect(a = all.ccre, b = G4, wa = TRUE, c = TRUE) %>% unique() %>% data.frame()
all.ccre.G4$V8 <- ifelse(all.ccre.G4$V8 > 0, 1, 0)

colnames(all.ccre.G4) <- c("chr", "start", "end", "IDD", "ID", "ccre", "encodeLabel", "contain_G4")

fwrite(all.ccre.G4, "../data/ccre/V3/G4_cCRE_annotation.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

