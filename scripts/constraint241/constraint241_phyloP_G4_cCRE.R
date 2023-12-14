
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

phylop241 <- fread("../data/conscore/241-mammalian-2020v2_threshold_2.27.wig", sep = "\t", header = FALSE) %>% data.frame()

GrData <- function(data = NULL, center = TRUE, score = 1) {
  data <- data[, 1:3]
  if (center == TRUE) {
    data[, 2] <- ceiling((data[, 2] + data[, 3]) /  2)
    data[, 3] <- data[, 2]
  }

  data$score <- score
  colnames(data)[1:3] <- c("chr", "start", "end")
  data <- makeGRangesFromDataFrame(data,
                                   keep.extra.columns = TRUE,
                                   ignore.strand = TRUE)

  return(data)
}

phylop241.gr <- GrData(phylop241, center = FALSE)

G4.Runs <- fread("../data/G4/G4_Runs.txt", sep = "\t", header = TRUE) %>% data.frame()
ccre.nonG4.Runs <- fread("../data/G4/nonG4_Runs_incCRE.txt", sep = "\t", header = TRUE) %>% data.frame()
ccre <- fread("../data/ccre/V3/G4_cCRE_annotation.txt", sep = "\t", header = TRUE) %>% data.frame() %>% filter(chr %in% paste0("chr", c(1:22, "X")))
G4 <- fread("../data/G4/cCRE_G4_annotation.txt", sep = "\t", header = TRUE) %>% data.frame()

extend.width <- 250
width = 1

for (tmp.tp in c("PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3")) {

  message(tmp.tp)
  tmp.G4 <- G4 %>% filter(get(str_replace(tmp.tp, "-", ".")) == 1)
  tmp.G4.gr <- GrData(tmp.G4)
  #tmp.ccre.G4 <- ccre %>% filter(encodeLabel == tmp.tp, contain_G4 == 1)
  tmp.ccre <- ccre %>% filter(encodeLabel == tmp.tp)
  tmp.ccre.gr <- GrData(tmp.ccre)  

  x <- normalizeToMatrix(phylop241.gr, tmp.G4.gr, extend = extend.width, mean_mode = "coverage", w = width)
  gc()

  z <- normalizeToMatrix(phylop241.gr, tmp.ccre.gr, extend = extend.width, mean_mode = "coverage", w = width)
  gc()

  mat <- bind_rows(data.frame(value = colMeans(x), group = paste0("G4s overlapped with ", tmp.tp), range = c(-(extend.width/width):(extend.width/width))),
                   data.frame(value = colMeans(z), group = tmp.tp, range = c(-(extend.width/width):(extend.width/width)))
                   )
  mat$group <- factor(mat$group, levels = c(paste0("G4s overlapped with ", tmp.tp), tmp.tp))
  
  p1 <- ggplot(mat, aes(x = range, y = value, group = group)) +
        geom_line(aes(color = group), linewidth = 0.5) +
        geom_vline(xintercept = 0, linetype = "dotted", color = "grey", linewidth = 0.35) + 
        theme_classic() + rremove("xlab") + theme(legend.position = "top", legend.direction = "vertical") + 
        ylab("Constraint base density") + rremove("legend.title") +
        scale_color_manual(values = paletteer_d("ggsci::default_aaas"))
  ggsave(paste0("../figure/constrained/phyloP/phyloP_G4_cCRE_", tmp.tp, ".pdf"), p1, width = 4.5, height = 4)

}


round <- 100
sample.n <- 1000

for (tmp.tp in c("PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3")) {

  message(tmp.tp)
  
  mat <- NULL

  tmp.G4.Runs <- G4.Runs %>% filter(encodeLabel == tmp.tp)
  tmp.G4.Runs.gr <- GrData(tmp.G4.Runs)
  tmp.ccre.nonG4.Runs <- ccre.nonG4.Runs %>% filter(encodeLabel == tmp.tp)
  tmp.ccre.nonG4.Runs.gr <- GrData(tmp.ccre.nonG4.Runs)

  for (tmp.round in 1:round) {
    message(tmp.round)
    x <- normalizeToMatrix(phylop241.gr, tmp.G4.Runs.gr[sample(1:length(tmp.G4.Runs.gr), sample.n)], extend = extend.width, mean_mode = "coverage", w = width)
    gc()

    y <- normalizeToMatrix(phylop241.gr, tmp.ccre.nonG4.Runs.gr[sample(1:length(tmp.ccre.nonG4.Runs.gr), sample.n)], extend = extend.width, mean_mode = "coverage", w = width)
    gc()

    mat <- bind_rows(mat,
                     data.frame(value = colMeans(x), group = paste0("G-Runs in ", tmp.tp," G4s"), range = c(-(extend.width/width):(extend.width/width)), round = tmp.round),
                     data.frame(value = colMeans(y), group = paste0("G-Runs in ", tmp.tp), range = c(-(extend.width/width):(extend.width/width)), round = tmp.round))

  }

  mat$group <- factor(mat$group, levels = c(paste0("G-Runs in ", tmp.tp," G4s"), paste0("G-Runs in ", tmp.tp)))

  p1 <- ggplot(mat, aes(x = range, y = value, group = group, color = group)) +
        #geom_line(aes(color = group), linewidth = 0.5) +
        stat_summary(fun.data = 'mean_sdl',
               fun.args = list(mult = 1), linewidth = 0.5, 
               geom = 'smooth', se = TRUE) +
        geom_vline(xintercept = 0, linetype = "dotted", color = "grey", linewidth = 0.35) + 
        theme_classic() + rremove("xlab") + theme(legend.position="top", legend.direction = "vertical") + 
        ylab("Constraint base density") + rremove("legend.title") +
        scale_color_manual(values = paletteer_d("ggsci::default_aaas"))
  ggsave(paste0("../figure/constrained/phyloP/phyloP_G4_cCRE_GRuns_", tmp.tp, ".pdf"), p1, width = 4.5, height = 4)

}


round <- 100
sample.n <- 1000
mat <- NULL

for (tmp.tp in c("PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3")) {

  message(tmp.tp)
  
  tmp.G4.Runs <- G4.Runs %>% filter(encodeLabel == tmp.tp)
  tmp.ccre.nonG4.Runs <- ccre.nonG4.Runs %>% filter(encodeLabel == tmp.tp)

  for (tmp.round in 1:round) {
    message(tmp.round)
    
    tmp.x <- bt.intersect(a = tmp.G4.Runs[sample(1:nrow(tmp.G4.Runs), sample.n),], b = phylop241, c = TRUE)
    tmp.y <- bt.intersect(a = tmp.ccre.nonG4.Runs[sample(1:nrow(tmp.ccre.nonG4.Runs), sample.n),], b = phylop241, c = TRUE)
    tmp.x <- sum(tmp.x[, 7]) / sum(tmp.x[, 3] - tmp.x[, 2] + 1)
    tmp.y <- sum(tmp.y[, 7]) / sum(tmp.y[, 3] - tmp.y[, 2] + 1)

    gc()

    mat <- bind_rows(mat,
                     data.frame(value = tmp.x, group = "G4 G-Runs", round = tmp.round, subgroup = tmp.tp),
                     data.frame(value = tmp.y, group = "cCRE G-Runs", round = tmp.round, subgroup = tmp.tp))

  }

}
mat$group <- factor(mat$group, levels = c("G4 G-Runs", "cCRE G-Runs"))
mat$subgroup <- factor(mat$subgroup, levels = c("PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3"))

save(list=ls(), file = "constraint241_phyloP_G4_cCRE_2.RData")

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
 return(data_sum)
}

mat2 <- data_summary(mat, varname = "value", 
                     groupnames = c("group", "subgroup"))

fig <- ggplot(mat2, aes(x = subgroup, y = value, fill = group)) + 
       geom_bar(stat = "identity", color = "black", position = position_dodge()) +
       geom_errorbar(aes(ymin = value-sd, ymax = value+sd), width = .2, position = position_dodge(.9)) +
       theme_classic() +
       rremove("xlab") + rremove("legend.title") + ylab("Proportion of constrained bases") + rotate_x_text(90) +
       scale_fill_manual(values = c("#B22C2C", "#2C85B2"))
ggsave(paste0("../figure/constrained/phyloP/phyloP241_G4_cCRE_GRuns_Proportion.pdf"), fig, width = 5, height = 3.8)


