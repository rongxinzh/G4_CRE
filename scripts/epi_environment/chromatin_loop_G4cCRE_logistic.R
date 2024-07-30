
set.seed(1)
options(scipen=200)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(bedtoolsr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(paletteer))
suppressPackageStartupMessages(library(broom))
suppressPackageStartupMessages(library(RColorBrewer))

hg19tohg38 <- function(data = NULL) {
  # liftOver oldFile map.chain newFile unMapped
  fwrite(data, "./tmp_hg19_ELS.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  system(paste0("liftOver ./tmp_hg19_ELS.txt ../data/ref/hg19ToHg38.over.chain.gz ./tmp_hg38_ELS.txt unMapped"))
  map.data <- fread("./tmp_hg38_ELS.txt", sep = "\t", header = FALSE) %>% data.frame()
  unlink("./tmp_hg19_ELS.txt")
  unlink("./tmp_hg38_ELS.txt")
  unlink("unMapped")
  return(map.data)
}

# load signals
enh.sig <- fread("../data/ccre/V3/H3K27ac.avg.txt", sep = " ", header = FALSE) %>% data.frame()
colnames(enh.sig) <- c("IDD", "value")

# load ccre
ccre <- fread("../data/ccre/V3/G4_cCRE_annotation.txt", sep = "\t", header = TRUE) %>% data.frame() %>% filter(chr %in% paste0("chr", c(1:22, "X")))
ccre$contain_G4 <- factor(ccre$contain_G4)
ccre <- ccre %>% left_join(enh.sig, by = "IDD") %>% data.frame()

ccre.pls <- ccre %>% filter(encodeLabel %in% c("PLS"))
ccre.pels <- ccre %>% filter(encodeLabel %in% c("pELS"))
ccre.dels <- ccre %>% filter(encodeLabel %in% c("dELS"))

loop <- fread("../data/TAD/hg19/41586_2020_2151_MOESM5_ESM.txt", sep = "\t", header = TRUE) %>% data.frame()
loop1 <- loop[, 1:3]
loop2 <- loop[, 4:6]
colnames(loop1) <- colnames(loop2) <- c("chr", "start", "end")
loop.anchor <- bind_rows(loop1, loop2)
loop.anchor <- hg19tohg38(loop.anchor)

ccre.pls.loop <- bt.intersect(a = ccre.pls, b = loop.anchor, wa = TRUE, c = TRUE)
ccre.pls.loop[, 10] <- ifelse(ccre.pls.loop[, 10] == 0, 0, 1)

ccre.pels.loop <- bt.intersect(a = ccre.pels, b = loop.anchor, wa = TRUE, c = TRUE)
ccre.pels.loop[, 10] <- ifelse(ccre.pels.loop[, 10] == 0, 0, 1)

ccre.dels.loop <- bt.intersect(a = ccre.dels, b = loop.anchor, wa = TRUE, c = TRUE)
ccre.dels.loop[, 10] <- ifelse(ccre.dels.loop[, 10] == 0, 0, 1)

model1 <- glm(V10 ~ V8+V9, data = ccre.pls.loop, family = binomial)
model2 <- glm(V10 ~ V8+V9, data = ccre.pels.loop, family = binomial)
model3 <- glm(V10 ~ V8+V9, data = ccre.dels.loop, family = binomial)

model1.info <- tidy(model1) %>% data.frame() 
model2.info <- tidy(model2) %>% data.frame() 
model3.info <- tidy(model3) %>% data.frame()

model.res <- bind_rows(data.frame(cCRE = "PLS", Estimate = model1.info[2:3, 2], PValue = model1.info[2:3, 5], Feature = c("G4", "H3K27ac")),
                       data.frame(cCRE = "pELS", Estimate = model2.info[2:3, 2], PValue = model2.info[2:3, 5], Feature = c("G4", "H3K27ac")),
                       data.frame(cCRE = "dELS", Estimate = model3.info[2:3, 2], PValue = model3.info[2:3, 5], Feature = c("G4", "H3K27ac")))
model.res$cCRE <- factor(model.res$cCRE, levels = c("PLS", "pELS", "dELS"))

model.res <- model.res %>%
  mutate(Significance = case_when(
    PValue < 0.0001 ~ "****",
    PValue < 0.001  ~ "***",
    PValue < 0.01   ~ "**",
    PValue < 0.05   ~ "*",
    TRUE            ~ "ns"
  ))

p1 <- ggplot(model.res, aes(x = cCRE, y = Estimate, fill = Feature)) +
        geom_bar(stat = "identity", position = "dodge") +
        geom_text(aes(label = Significance), 
                  position = position_dodge(width = 0.95), vjust = -0.1, size = 4) +
        labs(x = "", 
             y = "Estimate",
             fill = "") +
        scale_fill_manual(values = c("G4" = "#cb5815", "H3K27ac" = "#929293")) +
        theme_minimal()

ggsave("../figure/loop_model.pdf", p1, width = 4.5, height = 3.5)
