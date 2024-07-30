
set.seed(1)
options(scipen=200)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(paletteer))
suppressPackageStartupMessages(library(gridExtra))

load("./TSCNE-PLS.RData")
load("./TSCNE-ELS.RData")

ct.map <- data.frame(ct = unique(pls.mat$ct), abbr = c("STAD", "COAD/READ", "PRAD", "KIRC", "UCS", "BRCA", "PAAD", "HNSC"))

pls.mat <- pls.mat %>% left_join(ct.map, by = "ct")
els.mat <- els.mat %>% left_join(ct.map, by = "ct")

pls.mat$group <- ifelse(pls.mat$group == "G4 cCREs", "G4-associated cCREs", pls.mat$group)
els.mat$group <- ifelse(els.mat$group == "G4 cCREs", "G4-associated cCREs", els.mat$group)

p1 <- ggplot(pls.mat %>% filter(cCRE == "PLS"), aes(x = range, y = value, group = group)) +
        geom_line(aes(color = group), linewidth = 0.8) +
        theme_classic() + rremove("xlab") + ylab("Gained CRE in cancer\nwith H3K4me3 marker") + rremove("legend.title") +
        theme(legend.position="top") + 
        scale_x_continuous(breaks = c(-250, 0, 250), labels = c("-5kb", "PLS", "5kb")) + 
        scale_color_manual(values = paletteer_d("ggthemes::excel_Main_Event")) +
        facet_wrap(vars(abbr), scale = "free_y", ncol = 8)
p2 <- ggplot(els.mat %>% filter(cCRE == "pELS"), aes(x = range, y = value, group = group)) +
        geom_line(aes(color = group), linewidth = 0.8) +
        theme_classic() + rremove("xlab") + ylab("Gained CRE in cancer\nwith H3K27ac marker") + rremove("legend.title") +
        theme(legend.position="top") + 
        scale_x_continuous(breaks = c(-250, 0, 250), labels = c("-5kb", "pELS", "5kb")) + 
        scale_color_manual(values = paletteer_d("ggthemes::excel_Main_Event")) +
        facet_wrap(vars(abbr), scale = "free_y", ncol = 8) 
p3 <- ggplot(els.mat %>% filter(cCRE == "dELS"), aes(x = range, y = value, group = group)) +
        geom_line(aes(color = group), linewidth = 0.8) +
        theme_classic() + rremove("xlab") + ylab("Gained CRE in cancer\nwith H3K27ac marker") + rremove("legend.title") +
        theme(legend.position="top") + 
        scale_x_continuous(breaks = c(-250, 0, 250), labels = c("-5kb", "dELS", "5kb")) + 
        scale_color_manual(values = paletteer_d("ggthemes::excel_Main_Event")) +
        facet_wrap(vars(abbr), scale = "free_y", ncol = 8) 

ggsave("../figure/gained_CRE.pdf", plot = grid.arrange(p1, p2, p3, ncol = 1), width = 18, height = 7)
