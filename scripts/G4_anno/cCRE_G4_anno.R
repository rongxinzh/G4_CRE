
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

# load G4
G4 <- fread("../data/G4/G4Hunter_w25_s1.5_hg38.txt", sep = "\t", header = FALSE) %>% data.frame()
G4$id <- paste0("G4_", 1:dim(G4)[1])

G4 <- bt.intersect(a = G4, b = pls.ccre, wa = TRUE, c = TRUE) %>% unique() %>% data.frame()
G4 <- bt.intersect(a = G4, b = pels.ccre, wa = TRUE, c = TRUE) %>% unique() %>% data.frame()
G4 <- bt.intersect(a = G4, b = dels.ccre, wa = TRUE, c = TRUE) %>% unique() %>% data.frame()
G4 <- bt.intersect(a = G4, b = ctcf.ccre, wa = TRUE, c = TRUE) %>% unique() %>% data.frame()
G4 <- bt.intersect(a = G4, b = dh.ccre, wa = TRUE, c = TRUE) %>% unique() %>% data.frame()

G4$V11 <- ifelse(G4$V11 > 0, 1, 0)
G4$V12 <- ifelse(G4$V12 > 0, 1, 0)
G4$V13 <- ifelse(G4$V13 > 0, 1, 0)
G4$V14 <- ifelse(G4$V14 > 0, 1, 0)
G4$V15 <- ifelse(G4$V15 > 0, 1, 0)

G4$V16 <- rowSums(G4[, 11:15])

colnames(G4) <- c("chr", "start", "end", "width", "strand", "score", "max_score", "threshold", "window", "G4_ID", 
                  "PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3", "overlap_count")

fwrite(G4, "../data/G4/cCRE_G4_annotation.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
