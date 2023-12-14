
set.seed(1)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(bedtoolsr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyr))
#suppressPackageStartupMessages(library(ggstatsplot))
suppressPackageStartupMessages(library(EnrichedHeatmap))
suppressPackageStartupMessages(library(ComplexHeatmap))

# load ccre
ccre <- fread("../data/ccre/V3/G4_cCRE_annotation.txt", sep = "\t", header = TRUE) %>% data.frame() %>% filter(chr %in% paste0("chr", c(1:22, "X")))
colnames(ccre)[4] <- "rDHS"

# load G4
G4 <- fread("../data/G4/cCRE_G4_annotation.txt", sep = "\t", header = TRUE) %>% data.frame()

zscore.path <- paste0("../data/ccre/V3/", 
                      c("GRCh38.CTCF-zscore.rDHS-V3.txt", "GRCh38.DNase-zscore.rDHS-V3.txt", 
                        "GRCh38.H3K4me3-zscore.rDHS-V3.txt", "GRCh38.H3K27ac-zscore.rDHS-V3.txt"))

CalP <- function(mat = NULL) {
  
  mat.p <- NULL
  for (i in 9:dim(mat)[2]) {
    p1 <- wilcox.test(mat[mat[, 8] > 0, i] %>% as.numeric, mat[mat[, 8] == 0, i] %>% as.numeric, alternative = "greater")$p.value
    p2 <- wilcox.test(mat[mat[, 8] > 0, i] %>% as.numeric, mat[mat[, 8] == 0, i] %>% as.numeric, alternative = "less")$p.value
    #tmp.p <- ifelse(p1 < 0.05, "Greater", ifelse(p2 < 0.05, "Lesser", "ns"))
    tmp.p <- ifelse(p1 < 0.05, 1, ifelse(p2 < 0.05, -1, 0))
    mat.p <- c(mat.p, tmp.p)
  }

  return(mat.p)
}

res.list <- NULL

for (tmp.f in zscore.path) {
  message(tmp.f)
  tmp.signal <- fread(tmp.f, sep = "\t", header = TRUE) %>% data.frame
  tmp.signal <- tmp.signal[-1, ]

  pls.ccre.signal <- ccre %>% filter(encodeLabel == "PLS") %>% left_join(tmp.signal, by = "rDHS")
  pels.ccre.signal <- ccre %>% filter(encodeLabel == "pELS") %>% left_join(tmp.signal, by = "rDHS")
  dels.ccre.signal <- ccre %>% filter(encodeLabel == "dELS") %>% left_join(tmp.signal, by = "rDHS")
  ctcf.ccre.signal <- ccre %>% filter(encodeLabel == "CTCF-only") %>% left_join(tmp.signal, by = "rDHS")
  dh.ccre.signal <- ccre %>% filter(encodeLabel == "DNase-H3K4me3") %>% left_join(tmp.signal, by = "rDHS")

  pls.mat.p <- CalP(pls.ccre.signal)
  pels.mat.p <- CalP(pels.ccre.signal)
  dels.mat.p <- CalP(dels.ccre.signal)
  ctcf.mat.p <- CalP(ctcf.ccre.signal)
  dh.mat.p <- CalP(dh.ccre.signal)

  mat.p <- bind_cols(pls.mat.p, pels.mat.p, dels.mat.p, ctcf.mat.p, dh.mat.p) %>% as.matrix
  colnames(mat.p) <- c("PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3")
  mat.p <- t(mat.p)
  colnames(mat.p) <- colnames(tmp.signal)[2:dim(tmp.signal)[2]]

  res.list <- c(res.list, list(mat.p))

  tmp.signal <- 1 
  pls.ccre.signal <- 1 
  pels.ccre.signal <- 1 
  dels.ccre.signal <- 1 
  ctcf.ccre.signal <- 1 
  dh.ccre.signal <- 1 
  pls.mat.p <- 1 
  pels.mat.p <- 1 
  dels.mat.p <- 1 
  ctcf.mat.p <- 1 
  dh.mat.p <- 1 
  mat.p <- 1 

  rm(tmp.signal)
  rm(pls.ccre.signal)
  rm(pels.ccre.signal)
  rm(dels.ccre.signal)
  rm(ctcf.ccre.signal)
  rm(dh.ccre.signal)
  rm(pls.mat.p)
  rm(pels.mat.p)
  rm(dels.mat.p)
  rm(ctcf.mat.p)
  rm(dh.mat.p)
  rm(mat.p)

  gc()

}

save(list="res.list", file = "./z-score.RData")

message("DONE!")

all.data <- NULL
for (i in 1:length(zscore.path)) {
  tmp.f <- zscore.path[i]
  tmp.data <- gather(data.frame(cCRE = rownames(res.list[[i]]), res.list[[i]]), key = "file", value = "group", -cCRE) %>% 
              group_by(cCRE, group) %>% summarise(count = n()) %>% data.frame()
  tmp.data$group <- ifelse(tmp.data$group == 1, "Greater", ifelse(tmp.data$group == 0, "NS", "Less"))

  tmp.data$property <- str_split(tmp.f, "(\\.|-)", simplify = TRUE)[, 4]

  all.data <- bind_rows(all.data, tmp.data)

}

all.data <- all.data %>% filter(cCRE %in% c("PLS", "pELS", "dELS")) 
all.data$cCRE <- factor(all.data$cCRE, level = c("PLS", "pELS", "dELS"))
all.data$property <- factor(all.data$property, level = c("DNase", "H3K4me3", "H3K27ac", "CTCF"))

pdf(paste0("../figure/epigenome/cCREs_property_diff.pdf"), width = 8, height = 2.5)
p1 <- ggplot(all.data, aes(x = cCRE, y = count, fill = group)) +
      geom_col(colour = "black", position = "fill") +
      scale_y_continuous(labels = scales::percent) +
      scale_fill_brewer(palette = "Pastel1") + 
      ylab("Proportion") + rremove("xlab") + rotate_x_text(90) + labs(fill = "G4 cCREs VS other cCREs") +
      facet_wrap(vars(property), ncol = 4)
print(p1)
dev.off()

