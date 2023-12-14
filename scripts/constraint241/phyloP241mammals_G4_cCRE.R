
set.seed(1)
options(scipen=200)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(bedtoolsr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(random))
suppressPackageStartupMessages(library(paletteer))
suppressPackageStartupMessages(library(EnrichedHeatmap))
suppressPackageStartupMessages(library(tidyr))

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

ccre <- fread("../data/ccre/V3/G4_cCRE_annotation.txt", sep = "\t", header = TRUE) %>% data.frame() %>% filter(chr %in% paste0("chr", c(1:22, "X")))
G4 <- fread("../data/G4/cCRE_G4_annotation.txt", sep = "\t", header = TRUE) %>% data.frame()
G4.Runs <- fread("../data/G4/G4_Runs.txt", sep = "\t", header = TRUE) %>% data.frame()
ccre.nonG4.Runs <- fread("../data/G4/nonG4_Runs_incCRE.txt", sep = "\t", header = TRUE) %>% data.frame()

G4.phylop <- GetCSTR(G4 %>% filter(overlap_count <= 1) %>% select(chr:end, G4_ID) %>% data.frame(), "../data/conscore/241-mammalian-2020v2.bigWig")
colnames(G4.phylop)[4] <- "G4_ID"
G4.phylop.group <- G4 %>% inner_join(G4.phylop, by = "G4_ID")

G4.phylop.group <- bind_rows(data.frame(score = G4.phylop.group %>% filter(PLS > 0) %>% select(score.y) %>% data.frame() %>% unlist(), group = "PLS G4s"), 
                             data.frame(score = G4.phylop.group %>% filter(pELS > 0) %>% select(score.y) %>% data.frame() %>% unlist(), group = "pELS G4s"), 
                             data.frame(score = G4.phylop.group %>% filter(dELS > 0) %>% select(score.y) %>% data.frame() %>% unlist(), group = "dELS G4s"), 
                             data.frame(score = G4.phylop.group %>% filter(CTCF.only > 0) %>% select(score.y) %>% data.frame() %>% unlist(), group = "CTCF-only G4s"), 
                             data.frame(score = G4.phylop.group %>% filter(DNase.H3K4me3 > 0) %>% select(score.y) %>% data.frame() %>% unlist(), group = "DNase-H3K4me3 G4s"), 
                             data.frame(score = G4.phylop.group %>% filter(overlap_count == 0) %>% select(score.y) %>% data.frame() %>% unlist(), group = "Other G4s"))

fig <- ggerrorplot(G4.phylop.group, x = "group", y = "score", add = "mean_sd", color = "group",  ylab = "phyloP score",
                   palette = c("#e64b35", "#4ebad5", "#049f87", "#3c5488", "#f49b7f", "#7d6148")) +
       rremove("xlab") + rremove("legend") + rotate_x_text(90)
ggsave("../figure/constrained/phyloP/phyloP_G4groups.pdf", fig, width = 5, height = 5)

ccre.phylop <- GetCSTR(ccre %>% select(chr:end, ID) %>% data.frame(), "../data/conscore/241-mammalian-2020v2.bigWig")
colnames(ccre.phylop)[4] <- "ID"
ccre.phylop <- ccre %>% inner_join(ccre.phylop, by = "ID")

G4.ccre.phylop <- bind_rows(data.frame(subgroup = "G4", G4.phylop.group %>% filter(group != "Other G4s"), label = str_split(G4.phylop.group %>% filter(group != "Other G4s") %>% select(group) %>% unlist(), " ", simplify = TRUE)[, 1]),
                            data.frame(subgroup = "cCRE", score = ccre.phylop$score, group = paste0(ccre.phylop$encodeLabel, " cCRE"), label = ccre.phylop$encodeLabel))
G4.ccre.phylop$subgroup <- factor(G4.ccre.phylop$subgroup, levels=c("G4", "cCRE"))

wilcox.test(G4.ccre.phylop %>% filter(label == "DNase-H3K4me3", subgroup == "G4") %>% select(score) %>% unlist(),
            G4.ccre.phylop %>% filter(label == "DNase-H3K4me3", subgroup == "cCRE") %>% select(score) %>% unlist())$p.value

fig <- ggerrorplot(G4.ccre.phylop, x = "label", y = "score", add = "mean_sd", color = "subgroup",  ylab = "phyloP score",
                   palette = c("#b22222", "#7aa6dd")) +
       rremove("xlab") + rremove("legend.title") + rotate_x_text(90)
ggsave("../figure/constrained/phyloP/phyloP_G4_cCRE.pdf", fig, width = 5, height = 5)

G4.Runs$id <- paste0("x_", 1:dim(G4.Runs)[1])
ccre.nonG4.Runs$id <- paste0("y_", 1:dim(ccre.nonG4.Runs)[1])
##
G4.Runsx <- GetCSTR(G4.Runs[, c(1:3, 7)], "../data/conscore/241-mammalian-2020v2.bigWig")
ccre.nonG4.Runsy <- GetCSTR(ccre.nonG4.Runs[, c(1:3, 7)], "../data/conscore/241-mammalian-2020v2.bigWig")
G4.Runsx <- G4.Runsx %>% inner_join(G4.Runs, by = "id")
ccre.nonG4.Runsy <- ccre.nonG4.Runsy %>% inner_join(ccre.nonG4.Runs, by = "id")

G4.ccre.GRuns <- bind_rows(data.frame(G4.Runsx[, c("encodeLabel", "score")], group = "G4 G-Runs"),
                           data.frame(ccre.nonG4.Runsy[, c("encodeLabel", "score")], group = "cCRE G-Runs"))
G4.ccre.GRuns$group <- factor(G4.ccre.GRuns$group, levels = c("G4 G-Runs", "cCRE G-Runs"))
G4.ccre.GRuns$encodeLabel <- factor(G4.ccre.GRuns$encodeLabel, levels = c("PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3"))

fig <- ggerrorplot(G4.ccre.GRuns, x = "encodeLabel", y = "score", add = "mean_sd", color = "group",  ylab = "phyloP score",
                   palette = c("#b22222", "#7aa6dd")) +
       rremove("xlab") + rremove("legend.title") + rotate_x_text(90)
ggsave("../figure/constrained/phyloP/phyloP_GRuns_G4_cCRE.pdf", fig, width = 5, height = 5)

wilcox.test(G4.Runsx[G4.Runsx[, "encodeLabel"] == "DNase-H3K4me3", "score"],
            ccre.nonG4.Runsy[ccre.nonG4.Runsy[, "encodeLabel"] == "DNase-H3K4me3", "score"])$p.value

G4.Runs.filter <- bt.intersect(a = G4.Runs, b = G4 %>% filter(overlap_count > 1), wa = TRUE, v = TRUE) %>% unique() %>% data.frame()
colnames(G4.Runs.filter)[7] <- "id"
G4.use.ccre <- bt.intersect(a = G4 %>% filter(overlap_count == 1), b = ccre, wb = TRUE) %>% unique() %>% data.frame() 
G4.use.ccre <- G4.use.ccre %>% select(V1:V10, V23) %>% data.frame()
G4.remain <- bt.subtract(a = G4.use.ccre, b = G4.Runs.filter) %>% unique() %>% data.frame()
G4.remain$id <- paste0("G4L", 1:dim(G4.remain)[1])

G4.Runs.filter.value <- GetCSTR(G4.Runs.filter[, c(1:3, 7)], "../data/conscore/241-mammalian-2020v2.bigWig")
G4.Runs.filter.value <- G4.Runs.filter.value %>% left_join(G4.Runs.filter, by = "id")
G4.remain.value <- GetCSTR(G4.remain[, c(1:3, 12)], "../data/conscore/241-mammalian-2020v2.bigWig")
G4.remain.value <- G4.remain.value %>% left_join(G4.remain, by = "id")

G4.Runs.filter.value <- G4.Runs.filter.value[, c("score", "V5")]
G4.remain.value <- G4.remain.value[, c("score", "V11")]
colnames(G4.Runs.filter.value) <- colnames(G4.remain.value) <- c("score", "cCRE")

G4.RLvalue <- bind_rows(data.frame(G4.Runs.filter.value, group = "G4 G-Runs"),
                        data.frame(G4.remain.value, group = "Other nucleotides in G4s"))
G4.RLvalue$cCRE <- factor(G4.RLvalue$cCRE, levels = c("PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3"))

wilcox.test(G4.RLvalue %>% filter(cCRE == "DNase-H3K4me3", group == "G4 G-Runs") %>% select(score) %>% unlist(),
            G4.RLvalue %>% filter(cCRE == "DNase-H3K4me3", group == "Other nucleotides in G4s") %>% select(score) %>% unlist())$p.value

fig <- ggerrorplot(G4.RLvalue, x = "cCRE", y = "score", add = "mean_sd", color = "group",  ylab = "phyloP score",
                   palette = c("#b22222", "#7aa6dd")) +
       rremove("xlab") + rremove("legend.title") + rotate_x_text(90)
ggsave("../figure/constrained/phyloP/phyloP_GRuns_G4_nonG4Runs.pdf", fig, width = 5, height = 5)

