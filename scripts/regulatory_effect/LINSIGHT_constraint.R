
set.seed(1)
options(scipen=200)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(bedtoolsr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(random))
suppressPackageStartupMessages(library(paletteer))

GetCSTR <- function(elements = NULL, CSTR = NULL) {

  elements <- elements[, 1:4]
  #
  elements[, 2] <- elements[, 2] - 1

  rand.f <- randomStrings(n = 1, len = 20) %>% as.character
  fwrite(elements, rand.f, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  system(paste0("bigWigAverageOverBed ", CSTR," ", rand.f, " ", rand.f, ".out"))
  CSTR.score <- fread(paste0(rand.f, ".out"), sep = "\t", header = FALSE) %>% data.frame
  colnames(CSTR.score)[1] <- colnames(elements)[4] <- "id"
  CSTR.score <- elements %>% left_join(CSTR.score, by = "id")

  unlink(rand.f)
  unlink(paste0(rand.f, ".out"))

  CSTR.score[, 2] <- CSTR.score[, 2] + 1
  CSTR.score <- CSTR.score[, c(1:4, 9)]
  colnames(CSTR.score) <- c("chr", "start", "end", "id", "score")
  return(CSTR.score)
}

# load ccre
ccre <- fread("../data/ccre/V3/G4_cCRE_annotation.txt", sep = "\t", header = TRUE) %>% data.frame() %>% filter(chr %in% paste0("chr", c(1:22, "X")))

# load G4
G4 <- fread("../data/G4/cCRE_G4_annotation.txt", sep = "\t", header = TRUE) %>% data.frame()

lins.path <- "../data/LINSIGHT/LINSIGHT_hg38.bw"

G4.value <- GetCSTR(G4[G4[, "overlap_count"] <= 1, c(1:3, 10)], lins.path)
colnames(G4.value)[4] <- "G4_ID"
G4.value <- G4.value %>% left_join(G4, by = "G4_ID")

G4.lins.score <- bind_rows(data.frame(score = G4.value %>% filter(PLS == 1) %>% select(score.x) %>% unlist(), group = "PLS-G4s", G4score = G4.value %>% filter(PLS == 1) %>% select(max_score) %>% unlist() %>% abs()), 
                           data.frame(score = G4.value %>% filter(pELS == 1) %>% select(score.x) %>% unlist(), group = "pELS-G4s", G4score = G4.value %>% filter(pELS == 1) %>% select(max_score) %>% unlist() %>% abs()),
                           data.frame(score = G4.value %>% filter(dELS == 1) %>% select(score.x) %>% unlist(), group = "dELS-G4s", G4score = G4.value %>% filter(dELS == 1) %>% select(max_score) %>% unlist() %>% abs()),
                           data.frame(score = G4.value %>% filter(CTCF.only == 1) %>% select(score.x) %>% unlist(), group = "CTCF-only-G4s", G4score = G4.value %>% filter(CTCF.only == 1) %>% select(max_score) %>% unlist() %>% abs()),
                           data.frame(score = G4.value %>% filter(DNase.H3K4me3 == 1) %>% select(score.x) %>% unlist(), group = "DNase-H3K4me3-G4s", G4score = G4.value %>% filter(DNase.H3K4me3 == 1) %>% select(max_score) %>% unlist() %>% abs()),
                           data.frame(score = G4.value %>% filter(overlap_count == 0) %>% select(score.x) %>% unlist(), group = "Other-G4s", G4score = G4.value %>% filter(overlap_count == 0) %>% select(max_score) %>% unlist() %>% abs()))
G4.lins.score$group <- factor(G4.lins.score$group, level = unique(G4.lins.score$group))

kruskal.test(score ~ group, data = G4.lins.score)$p.value

pairwise.wilcox.test(G4.lins.score$score, G4.lins.score$group)$p.value

fig <- ggerrorplot(G4.lins.score, x = "group", y = "score", add = "mean_sd", color = "group",  ylab = "LINSIGHT score",
                   palette = c("#e64b35", "#4ebad5", "#049f87", "#3c5488", "#f49b7f", "#7d6148")) +
       rremove("xlab") + rremove("legend") + rotate_x_text(90)
ggsave("../figure/LINSIGHT/LINSIGHTscore_cCREG4s_VS_otherG4s.pdf", fig, width = 5, height = 6)

ccre.value <- GetCSTR(ccre[, c(1:3, 5)], lins.path)
colnames(ccre.value)[4] <- "ID"
ccre.value <- ccre.value %>% left_join(ccre, by = "ID")
ccre.value <- ccre.value[, c("score", "encodeLabel", "contain_G4")]
ccre.value$contain_G4 <- ifelse(ccre.value$contain_G4 > 0, "G4-associated", "Others")
#
wilcox.test(ccre.value %>% filter(encodeLabel == "DNase-H3K4me3", contain_G4 == "G4-associated") %>% select(score) %>% unlist(),
            ccre.value %>% filter(encodeLabel == "DNase-H3K4me3", contain_G4 == "Others") %>% select(score) %>% unlist())$p.value

fig <- ggerrorplot(ccre.value, x = "encodeLabel", y = "score", add = "mean_sd", color = "contain_G4",  ylab = "LINSIGHT score",
                   palette = c("firebrick", "goldenrod1")) +
       rremove("xlab") + rremove("legend.title") + rotate_x_text(90)
ggsave("../figure/LINSIGHT/LINSIGHTscore_G4cCREs_VS_othercCREs.pdf", fig, width = 5, height = 5)

G4.Runs <- fread("../data/G4/G4_Runs.txt", sep = "\t", header = TRUE) %>% data.frame()
ccre.nonG4.Runs <- fread("../data/G4/nonG4_Runs_incCRE.txt", sep = "\t", header = TRUE) %>% data.frame()
G4.Runs$id <- paste0("G4R", 1:dim(G4.Runs)[1])
ccre.nonG4.Runs$id <- paste0("nonR", 1:dim(ccre.nonG4.Runs)[1])

G4.Runs.value <- GetCSTR(G4.Runs[, c(1:3, 7)], lins.path)
G4.Runs.value <- G4.Runs.value %>% left_join(G4.Runs, by = "id")
ccre.nonG4.Runs.value <- GetCSTR(ccre.nonG4.Runs[, c(1:3, 7)], lins.path)
ccre.nonG4.Runs.value <- ccre.nonG4.Runs.value %>% left_join(ccre.nonG4.Runs, by = "id")

G4.Runs.value <- G4.Runs.value[, c("score", "encodeLabel")]
ccre.nonG4.Runs.value <- ccre.nonG4.Runs.value[, c("score", "encodeLabel")]

Runs.value <- bind_rows(data.frame(G4.Runs.value, group = "G4 G-Runs"),
                        data.frame(ccre.nonG4.Runs.value, group = "cCREs G-Runs"))
Runs.value$encodeLabel <- factor(Runs.value$encodeLabel, levels = c("PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3"))
#
wilcox.test(Runs.value %>% filter(encodeLabel == "DNase-H3K4me3", group == "G4 G-Runs") %>% select(score) %>% unlist(),
            Runs.value %>% filter(encodeLabel == "DNase-H3K4me3", group == "cCREs G-Runs") %>% select(score) %>% unlist())$p.value

fig <- ggerrorplot(Runs.value, x = "group", y = "score", add = "mean_sd", color = "group",  ylab = "LINSIGHT score",
                   palette = c("firebrick", "goldenrod1")) +
       rremove("xlab") + rremove("legend.title") + rotate_x_text(90) + 
       facet_wrap(vars(encodeLabel), ncol = 5)
ggsave("../figure/LINSIGHT/LINSIGHTscore_G4Runs_VS_cCRERuns.pdf", fig, width = 7, height = 5.86)

G4.Runs.filter <- bt.intersect(a = G4.Runs, b = G4 %>% filter(overlap_count > 1), wa = TRUE, v = TRUE) %>% unique() %>% data.frame()
colnames(G4.Runs.filter)[7] <- "id"
G4.use.ccre <- bt.intersect(a = G4 %>% filter(overlap_count == 1), b = ccre, wb = TRUE) %>% unique() %>% data.frame() 
G4.use.ccre <- G4.use.ccre %>% select(V1:V10, V23) %>% data.frame()
G4.remain <- bt.subtract(a = G4.use.ccre, b = G4.Runs.filter) %>% unique() %>% data.frame()
G4.remain$id <- paste0("G4L", 1:dim(G4.remain)[1])

G4.Runs.filter.value <- GetCSTR(G4.Runs.filter[, c(1:3, 7)], lins.path)
G4.Runs.filter.value <- G4.Runs.filter.value %>% left_join(G4.Runs.filter, by = "id")
G4.remain.value <- GetCSTR(G4.remain[, c(1:3, 12)], lins.path)
G4.remain.value <- G4.remain.value %>% left_join(G4.remain, by = "id")

G4.Runs.filter.value <- G4.Runs.filter.value[, c("score", "V5")]
G4.remain.value <- G4.remain.value[, c("score", "V11")]
colnames(G4.Runs.filter.value) <- colnames(G4.remain.value) <- c("score", "cCRE")

G4.RLvalue <- bind_rows(data.frame(G4.Runs.filter.value, group = "G4 G-Runs"),
                        data.frame(G4.remain.value, group = "Other nucleotides in G4s"))
G4.RLvalue$cCRE <- factor(G4.RLvalue$cCRE, levels = c("PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3"))

#
wilcox.test(G4.RLvalue %>% filter(cCRE == "DNase-H3K4me3", group == "G4 G-Runs") %>% select(score) %>% unlist(),
            G4.RLvalue %>% filter(cCRE == "DNase-H3K4me3", group == "Other nucleotides in G4s") %>% select(score) %>% unlist())$p.value

fig <- ggerrorplot(G4.RLvalue, x = "group", y = "score", add = "mean_sd", color = "group",  ylab = "LINSIGHT score",
                   palette = c("firebrick", "goldenrod1")) +
       rremove("xlab") + rremove("legend") + rotate_x_text(90) + 
       facet_wrap(vars(cCRE), ncol = 5)
ggsave("../figure/LINSIGHT/LINSIGHTscore_G4Runs_VS_G4Loop.pdf", fig, width = 7, height = 6)
