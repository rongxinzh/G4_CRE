
set.seed(1)
options(scipen=200)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(bedtoolsr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(vcfR))
suppressPackageStartupMessages(library(random))
suppressPackageStartupMessages(library(paletteer))

source("./G4Hunter.R")

# load ccre
ccre <- fread("../data/ccre/V3/G4_cCRE_annotation.txt", sep = "\t", header = TRUE) %>% data.frame() %>% filter(chr %in% paste0("chr", c(1:22, "X")))

# load G4
G4 <- fread("../data/G4/cCRE_G4_annotation.txt", sep = "\t", header = TRUE) %>% data.frame()

regdb <- fread("../data/regulomedb/ENCFF250UJY.tsv", sep = "\t", header = TRUE) %>% data.frame()
#
regdb$regdbid <- paste0("v_", 1:dim(regdb)[1]) 

reg.G4 <- bt.intersect(a = regdb, b = G4 %>% filter(overlap_count <= 1), wa = TRUE, wb = TRUE) %>% unique() %>% data.frame()

reg.G4.score <- bind_rows(data.frame(score = reg.G4 %>% filter(V27 == 1) %>% select(V6) %>% unlist(), group = "PLS-G4s"), 
                          data.frame(score = reg.G4 %>% filter(V28 == 1) %>% select(V6) %>% unlist(), group = "pELS-G4s"),
                          data.frame(score = reg.G4 %>% filter(V29 == 1) %>% select(V6) %>% unlist(), group = "dELS-G4s"),
                          data.frame(score = reg.G4 %>% filter(V30 == 1) %>% select(V6) %>% unlist(), group = "CTCF-only-G4s"),
                          data.frame(score = reg.G4 %>% filter(V31 == 1) %>% select(V6) %>% unlist(), group = "DNase-H3K4me3-G4s"),
                          data.frame(score = reg.G4 %>% filter(V32 == 0) %>% select(V6) %>% unlist(), group = "Other-G4s"))
reg.G4.score$group <- factor(reg.G4.score$group, level = unique(reg.G4.score$group))

kruskal.test(score ~ group, data = reg.G4.score)$p.value

pairwise.wilcox.test(reg.G4.score$score, reg.G4.score$group)$p.value

fig <- ggboxplot(data = reg.G4.score, x = 'group', y = 'score', fill = 'group', 
                 palette = c("#e64b35", "#4ebad5", "#049f87", "#3c5488", "#f49b7f", "#7d6148"), 
                 ylab = "RegulomeDB score", width = .5) + 
       rremove("xlab") + rremove("legend") + rotate_x_text(90)
ggsave("../figure/RegulomeDB/regscore_cCREG4s_VS_otherG4s.pdf", fig, width = 4, height = 5)

reg.G4.rank <- bind_rows(data.frame(rank = reg.G4 %>% filter(V27 == 1) %>% select(V5) %>% unlist(), group = "PLS-G4s"), 
                         data.frame(rank = reg.G4 %>% filter(V28 == 1) %>% select(V5) %>% unlist(), group = "pELS-G4s"),
                         data.frame(rank = reg.G4 %>% filter(V29 == 1) %>% select(V5) %>% unlist(), group = "dELS-G4s"),
                         data.frame(rank = reg.G4 %>% filter(V30 == 1) %>% select(V5) %>% unlist(), group = "CTCF-only-G4s"),
                         data.frame(rank = reg.G4 %>% filter(V31 == 1) %>% select(V5) %>% unlist(), group = "DNase-H3K4me3-G4s"),
                         data.frame(rank = reg.G4 %>% filter(V32 == 0) %>% select(V5) %>% unlist(), group = "Other-G4s"))
reg.G4.rank$group <- factor(reg.G4.rank$group, level = unique(reg.G4.rank$group))
reg.G4.rank <- reg.G4.rank %>% group_by(rank, group) %>% summarise(count = n()) %>% data.frame()
reg.G4.rank$rank <- factor(reg.G4.rank$rank, levels = unique(reg.G4.rank$rank))

fig <- ggplot(reg.G4.rank, aes(fill = rank, y = count, x = group)) + 
       geom_bar(position = "fill", stat = "identity") +
       theme_light() +
       ylab("RegulomeDB Ranking") + rremove("xlab") + rotate_x_text(90) +
       scale_fill_manual(values = paletteer_d("ggthemes::Tableau_20"))
ggsave("../figure/RegulomeDB/regrank_cCREG4s_VS_otherG4s.pdf", fig, width = 5, height = 5.5)

regdb.ccre <- bt.intersect(a = regdb, b = ccre, wa = TRUE, wb = TRUE) %>% unique() %>% data.frame()

regdb.ccre.rank <- regdb.ccre[, c("V5", "V23", "V24")]
colnames(regdb.ccre.rank) <- c("rank", "cCRE", "group")
regdb.ccre.rank$group <- ifelse(regdb.ccre.rank$group == 1, "G4-associated", "Others")
regdb.ccre.rank$cCRE <- factor(regdb.ccre.rank$cCRE, levels = c("PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3"))
regdb.ccre.rank <- regdb.ccre.rank %>% group_by(rank, group, cCRE) %>% summarise(count = n()) %>% data.frame()
regdb.ccre.rank$rank <- factor(regdb.ccre.rank$rank, levels = unique(regdb.ccre.rank$rank))
fig <- ggplot(regdb.ccre.rank, aes(fill = rank, y = count, x = group)) + 
       geom_bar(position="fill", stat="identity") +
       theme_light() +
       ylab("RegulomeDB Ranking") + rremove("xlab") + rotate_x_text(90) +
       scale_fill_manual(values = paletteer_d("ggthemes::Tableau_20")) +
       facet_wrap(vars(cCRE), ncol = 5)
ggsave("../figure/RegulomeDB/regrank_G4cCREs_VS_othercCREs.pdf", fig, width = 8, height = 5)

regdb.ccre <- regdb.ccre[, c("V6", "V23", "V24")]
colnames(regdb.ccre) <- c("score", "cCRE", "group")
regdb.ccre$group <- ifelse(regdb.ccre$group == 1, "G4-associated", "Others")
regdb.ccre$cCRE <- factor(regdb.ccre$cCRE, levels = c("PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3"))
#
wilcox.test(regdb.ccre %>% filter(cCRE == "DNase-H3K4me3", group == "G4-associated") %>% select(score) %>% unlist(),
            regdb.ccre %>% filter(cCRE == "DNase-H3K4me3", group == "Others") %>% select(score) %>% unlist())$p.value

fig <- ggboxplot(data = regdb.ccre, x = 'cCRE', y = 'score', fill = 'group', 
                 palette = c("#ce534c", "#7aa6dd"), 
                 ylab = "RegulomeDB score", width = .5) + 
       rremove("xlab") + rremove("legend.title") + rotate_x_text(90)
ggsave("../figure/RegulomeDB/regscore_G4cCREs_VS_othercCREs.pdf", fig, width = 5, height = 5)

G4.Runs <- fread("../data/G4/G4_Runs.txt", sep = "\t", header = TRUE) %>% data.frame()
G4.Runs.filter <- bt.intersect(a = G4.Runs, b = G4 %>% filter(overlap_count > 1), wa = TRUE, v = TRUE) %>% unique() %>% data.frame()
regdb.G4.Runs <- bt.intersect(a = regdb, b = G4.Runs.filter, wa = TRUE, wb = TRUE) %>% unique() %>% data.frame()
regdb.G4.Runs <- regdb.G4.Runs[, c("V16", "V5", "V6", "V21")]
regdb.G4 <- bt.intersect(a = regdb, b = G4 %>% filter(overlap_count == 1), wa = TRUE, wb = TRUE) %>% unique() %>% data.frame()
regdb.G4 <- regdb.G4 %>% filter(!(V16 %in% regdb.G4.Runs$V16)) %>% unique()

reg.ccreG4.score <- bind_rows(data.frame(score = regdb.G4.Runs[, "V6"], group = "G4 G-Runs", cCRE = regdb.G4.Runs[, "V21"]),
                              data.frame(score = regdb.G4 %>% filter(V27 == 1) %>% select(V6) %>% unlist(), group = "Other nucleotides in G4s", cCRE = "PLS"),
                              data.frame(score = regdb.G4 %>% filter(V28 == 1) %>% select(V6) %>% unlist(), group = "Other nucleotides in G4s", cCRE = "pELS"),
                              data.frame(score = regdb.G4 %>% filter(V29 == 1) %>% select(V6) %>% unlist(), group = "Other nucleotides in G4s", cCRE = "dELS"),
                              data.frame(score = regdb.G4 %>% filter(V30 == 1) %>% select(V6) %>% unlist(), group = "Other nucleotides in G4s", cCRE = "CTCF-only"),
                              data.frame(score = regdb.G4 %>% filter(V31 == 1) %>% select(V6) %>% unlist(), group = "Other nucleotides in G4s", cCRE = "DNase-H3K4me3"))
reg.ccreG4.score$cCRE <- factor(reg.ccreG4.score$cCRE, levels = c("PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3"))

#
wilcox.test(reg.ccreG4.score %>% filter(cCRE == "DNase-H3K4me3", group == "G4 G-Runs") %>% select(score) %>% unlist(),
            reg.ccreG4.score %>% filter(cCRE == "DNase-H3K4me3", group == "Other nucleotides in G4s") %>% select(score) %>% unlist())$p.value

fig <- ggboxplot(data = reg.ccreG4.score, x = 'cCRE', y = 'score', fill = 'group', 
                 palette = c("#ce534c", "#7aa6dd"), 
                 ylab = "RegulomeDB score", width = .5) + 
       rremove("xlab") + rremove("legend.title") + rotate_x_text(90)
ggsave("../figure/RegulomeDB/regscore_G4Runs_VS_nonG4Runs.pdf", fig, width = 5, height = 5)

reg.ccreG4.rank <- bind_rows(data.frame(rank = regdb.G4.Runs[, "V5"], group = "G4 G-Runs", cCRE = regdb.G4.Runs[, "V21"]),
                             data.frame(rank = regdb.G4 %>% filter(V27 == 1) %>% select(V5) %>% unlist(), group = "Other nucleotides in G4s", cCRE = "PLS"),
                             data.frame(rank = regdb.G4 %>% filter(V28 == 1) %>% select(V5) %>% unlist(), group = "Other nucleotides in G4s", cCRE = "pELS"),
                             data.frame(rank = regdb.G4 %>% filter(V29 == 1) %>% select(V5) %>% unlist(), group = "Other nucleotides in G4s", cCRE = "dELS"),
                             data.frame(rank = regdb.G4 %>% filter(V30 == 1) %>% select(V5) %>% unlist(), group = "Other nucleotides in G4s", cCRE = "CTCF-only"),
                             data.frame(rank = regdb.G4 %>% filter(V31 == 1) %>% select(V5) %>% unlist(), group = "Other nucleotides in G4s", cCRE = "DNase-H3K4me3"))
reg.ccreG4.rank$cCRE <- factor(reg.ccreG4.rank$cCRE, levels = c("PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3"))
reg.ccreG4.rank <- reg.ccreG4.rank %>% group_by(rank, group, cCRE) %>% summarise(count = n()) %>% data.frame()
reg.ccreG4.rank$rank <- factor(reg.ccreG4.rank$rank, levels = unique(reg.ccreG4.rank$rank))

fig <- ggplot(reg.ccreG4.rank, aes(fill = rank, y = count, x = group)) + 
       geom_bar(position="fill", stat="identity") +
       theme_light() +
       ylab("RegulomeDB Ranking") + rremove("xlab") + rotate_x_text(90) +
       scale_fill_manual(values = paletteer_d("ggthemes::Tableau_20")) +
       facet_wrap(vars(cCRE), ncol = 5)
ggsave("../figure/RegulomeDB/regrank_G4Runs_VS_nonG4Runs.pdf", fig, width = 8, height = 5)

ccre.nonG4.Runs <- fread("../data/G4/nonG4_Runs_incCRE.txt", sep = "\t", header = TRUE) %>% data.frame()

reg.ccre.Runs <- bt.intersect(a = regdb, b = ccre.nonG4.Runs, wa = TRUE, wb = TRUE) %>% unique() %>% data.frame()
reg.ccre.Runs <- reg.ccre.Runs[, c("V5", "V6", "V21")]
regdb.G4.Runs <- bt.intersect(a = regdb, b = G4.Runs, wa = TRUE, wb = TRUE) %>% unique() %>% data.frame()
regdb.G4.Runs <- regdb.G4.Runs[, c("V5", "V6", "V21")]

colnames(reg.ccre.Runs) <- colnames(regdb.G4.Runs) <- c("rank", "score", "cCRE")
regdb.G4.Runs$group <- "G4 G-Runs"
reg.ccre.Runs$group <- "cCRE G-Runs"
reg.ccreGRuns.score <- bind_rows(regdb.G4.Runs, reg.ccre.Runs)
reg.ccreGRuns.score$cCRE <- factor(reg.ccreGRuns.score$cCRE, levels = c("PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3"))
reg.ccreGRuns.score$group <- factor(reg.ccreGRuns.score$group, levels = c("G4 G-Runs", "cCRE G-Runs"))

wilcox.test(reg.ccreGRuns.score %>% filter(cCRE == "DNase-H3K4me3", group == "G4 G-Runs") %>% select(score) %>% unlist(),
            reg.ccreGRuns.score %>% filter(cCRE == "DNase-H3K4me3", group == "cCRE G-Runs") %>% select(score) %>% unlist())$p.value

fig <- ggboxplot(data = reg.ccreGRuns.score, x = 'cCRE', y = 'score', fill = 'group', 
                 palette = c("#ce534c", "#7aa6dd"), 
                 ylab = "RegulomeDB score", width = .5) + 
       rremove("xlab") + rremove("legend.title") + rotate_x_text(90)
ggsave("../figure/RegulomeDB/regscore_G4Runs_VS_cCRERuns.pdf", fig, width = 5, height = 5)

reg.ccreGRuns.rank <- reg.ccreGRuns.score %>% group_by(rank, group, cCRE) %>% summarise(count = n()) %>% data.frame()
reg.ccreGRuns.rank$rank <- factor(reg.ccreGRuns.rank$rank, levels = unique(reg.ccreGRuns.rank$rank))

fig <- ggplot(reg.ccreGRuns.rank, aes(fill = rank, y = count, x = group)) + 
       geom_bar(position="fill", stat="identity") +
       theme_light() +
       ylab("RegulomeDB Ranking") + rremove("xlab") + rotate_x_text(90) +
       scale_fill_manual(values = paletteer_d("ggthemes::Tableau_20")) +
       facet_wrap(vars(cCRE), ncol = 5)
ggsave("../figure/RegulomeDB/regrank_G4Runs_VS_cCRERuns.pdf", fig, width = 8, height = 5)
