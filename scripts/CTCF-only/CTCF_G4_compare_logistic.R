
set.seed(1)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(bedtoolsr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(EnrichedHeatmap))
suppressPackageStartupMessages(library(paletteer))
suppressPackageStartupMessages(library(RColorBrewer))

# load ccre
ccre <- fread("../data/ccre/V3/G4_cCRE_annotation.txt", sep = "\t", header = TRUE) %>% data.frame() %>% filter(chr %in% paste0("chr", c(1:22, "X")))

enh.sig <- fread("../data/ccre/V3/H3K27ac.avg.txt", sep = " ", header = FALSE) %>% data.frame()
colnames(enh.sig) <- c("IDD", "H3K27ac")
prm.sig <- fread("../data/ccre/V3/H3K4me3.avg.txt", sep = " ", header = FALSE) %>% data.frame()
colnames(prm.sig) <- c("IDD", "H3K4me3")
ctc.sig <- fread("../data/ccre/V3/CTCF.avg.txt", sep = " ", header = FALSE) %>% data.frame()
colnames(ctc.sig) <- c("IDD", "CTCF")

ccre <- ccre %>% left_join(enh.sig, by = "IDD")
ccre <- ccre %>% left_join(prm.sig, by = "IDD")
ccre <- ccre %>% left_join(ctc.sig, by = "IDD")

ccre.pls <- ccre %>% filter(encodeLabel %in% c("PLS"))
ccre.pels <- ccre %>% filter(encodeLabel %in% c("pELS"))
ccre.dels <- ccre %>% filter(encodeLabel %in% c("dELS"))
ccre.dh <- ccre %>% filter(encodeLabel %in% c("DNase-H3K4me3"))

model1 <- summary(glm(contain_G4 ~ H3K4me3*CTCF, data = ccre.pls, family = binomial))$coefficients
model2 <- summary(glm(contain_G4 ~ H3K27ac*CTCF, data = ccre.pels, family = binomial))$coefficients
model3 <- summary(glm(contain_G4 ~ H3K27ac*CTCF, data = ccre.dels, family = binomial))$coefficients
model4 <- summary(glm(contain_G4 ~ H3K4me3*CTCF, data = ccre.dh, family = binomial))$coefficients

model.res <- bind_rows(data.frame(cCRE = "PLS", Feature = c("H3K4me3", "CTCF", "H3K4me3×CTCF"), Estimate = model1[2:4, 1], PValue = model1[2:4, 4]),
                       data.frame(cCRE = "pELS", Feature = c("H3K27ac", "CTCF", "H3K27ac×CTCF"), Estimate = model2[2:4, 1], PValue = model2[2:4, 4]),
                       data.frame(cCRE = "dELS", Feature = c("H3K27ac", "CTCF", "H3K27ac×CTCF"), Estimate = model3[2:4, 1], PValue = model3[2:4, 4]),
                       data.frame(cCRE = "DNase-H3K4me3", Feature = c("H3K4me3", "CTCF", "H3K4me3×CTCF"), Estimate = model4[2:4, 1], PValue = model4[2:4, 4]))

model.res$PValue[model.res$PValue <= 1e-100] <- 1e-100
model.res$Significant <- model.res$PValue > 0.05
model.res$Shape <- ifelse(model.res$Estimate > 0, 24, 25) 
model.res$AbsEstimate <- abs(model.res$Estimate)

model.res$cCRE <- factor(model.res$cCRE, levels = c("dELS", "pELS", "DNase-H3K4me3", "PLS"))
model.res$Feature <- factor(model.res$Feature, levels = c("H3K27ac", "H3K27ac×CTCF", "CTCF", "H3K4me3", "H3K4me3×CTCF"))

p <- ggplot(model.res, aes(x = Feature, y = cCRE)) +
  geom_point(aes(color = -log10(PValue), fill = -log10(PValue), shape = as.factor(Shape), size = AbsEstimate), stroke = 1) +
  scale_color_gradient(low = "#fbb2bc", high = "#981320", name = "Significance\n(-log10 P-value)") +
  scale_fill_gradient(low = "#fbb2bc", high = "#981320", name = "Significance\n(-log10 P-value)") +
  scale_size_continuous(name = "Absolute estimate", range = c(1, 5), 
                        guide = guide_legend(override.aes = list(shape = 24))) +
  scale_shape_manual(values = c("24" = 24, "25" = 25), name = "Estimate direction", labels = c("Positive", "Negative")) +
  theme_minimal() +
  labs(x = "", y = "") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  geom_point(data = model.res[model.res$Significant, ], aes(x = Feature, y = cCRE), shape = 4, size = 4, color = "grey", stroke = 0.8)
ggsave("../figure/ctcf_model.pdf", p, width = 5, height = 6)

