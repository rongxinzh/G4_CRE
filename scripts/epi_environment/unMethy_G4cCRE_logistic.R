
set.seed(1)
options(scipen=200)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(bedtoolsr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(paletteer))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(gridExtra))

extend.width <- 2500
width = 10

# load ccre
ccre <- fread("../data/ccre/V3/G4_cCRE_annotation.txt", sep = "\t", header = TRUE) %>% data.frame() %>% filter(chr %in% paste0("chr", c(1:22, "X")))

enh.sig <- fread("../data/ccre/V3/H3K27ac.avg.txt", sep = " ", header = FALSE) %>% data.frame()
colnames(enh.sig) <- c("IDD", "value")
prm.sig <- fread("../data/ccre/V3/H3K4me3.avg.txt", sep = " ", header = FALSE) %>% data.frame()
colnames(prm.sig) <- c("IDD", "value")
ctc.sig <- fread("../data/ccre/V3/CTCF.avg.txt", sep = " ", header = FALSE) %>% data.frame()
colnames(ctc.sig) <- c("IDD", "value")
dna.sig <- fread("../data/ccre/V3/DNase.avg.txt", sep = " ", header = FALSE) %>% data.frame()
colnames(dna.sig) <- c("IDD", "value")

ccre <- ccre %>% left_join(enh.sig, by = "IDD")
ccre <- ccre %>% left_join(prm.sig, by = "IDD")
ccre <- ccre %>% left_join(ctc.sig, by = "IDD")
ccre <- ccre %>% left_join(dna.sig, by = "IDD")

ccre.pls <- ccre %>% filter(encodeLabel %in% c("PLS"))
ccre.pels <- ccre %>% filter(encodeLabel %in% c("pELS"))
ccre.dels <- ccre %>% filter(encodeLabel %in% c("dELS"))

model.res <- NULL
methy.files <- list.files("../data/methy/41586_2022_5580_MOESM5_ESM/") 
for (tmp.methy.f in methy.files) {
  message(tmp.methy.f)
  tmp.methy <- fread(paste0("../data/methy/41586_2022_5580_MOESM5_ESM/", tmp.methy.f), sep = "\t", header = TRUE) %>% data.frame()

  tmp.ccre.pls <- bt.coverage(a = ccre.pls, b = tmp.methy)
  tmp.ccre.pls[, ncol(tmp.ccre.pls)] <- ifelse(tmp.ccre.pls[, ncol(tmp.ccre.pls)] > 0, 1, 0)

  tmp.ccre.pels <- bt.coverage(a = ccre.pels, b = tmp.methy)
  tmp.ccre.pels[, ncol(tmp.ccre.pels)] <- ifelse(tmp.ccre.pels[, ncol(tmp.ccre.pels)] > 0, 1, 0)

  tmp.ccre.dels <- bt.coverage(a = ccre.dels, b = tmp.methy)
  tmp.ccre.dels[, ncol(tmp.ccre.dels)] <- ifelse(tmp.ccre.dels[, ncol(tmp.ccre.dels)] > 0, 1, 0)

  tmp.pls.model <- glm(V16 ~ V8 + V9 + V10 + V11 + V12, data = tmp.ccre.pls, family = binomial)
  tmp.pels.model <- glm(V16 ~ V8 + V9 + V10 + V11 + V12, data = tmp.ccre.pels, family = binomial)
  tmp.dels.model <- glm(V16 ~ V8 + V9 + V10 + V11 + V12, data = tmp.ccre.dels, family = binomial)

  model.res <- bind_rows(model.res, 
                         data.frame(cCRE = "PLS", tissue = str_split(tmp.methy.f, "\\.", simplify = TRUE)[, 1],
                                    Estimate = as.numeric(summary(tmp.pls.model)$coefficients[2:6, 1]), 
                                    Feature = c("G4", "H3K27ac", "H3K4me3", "CTCF", "DNase"),
                                    PValue = as.numeric(summary(tmp.pls.model)$coefficients[2:6, 4])),
                         data.frame(cCRE = "pELS", tissue = str_split(tmp.methy.f, "\\.", simplify = TRUE)[, 1],
                                    Estimate = as.numeric(summary(tmp.pels.model)$coefficients[2:6, 1]), 
                                    Feature = c("G4", "H3K27ac", "H3K4me3", "CTCF", "DNase"), 
                                    PValue = as.numeric(summary(tmp.pels.model)$coefficients[2:6, 4])),
                         data.frame(cCRE = "dELS", tissue = str_split(tmp.methy.f, "\\.", simplify = TRUE)[, 1],
                                    Estimate = as.numeric(summary(tmp.dels.model)$coefficients[2:6, 1]), 
                                    Feature = c("G4", "H3K27ac", "H3K4me3", "CTCF", "DNase"), 
                                    PValue = as.numeric(summary(tmp.dels.model)$coefficients[2:6, 4])))
  print("*******************")
}

model.res$Significant <- ifelse(model.res$PValue > 0.05, TRUE, FALSE)
model.res$Feature <- factor(model.res$Feature, levels = c("G4", "H3K27ac", "H3K4me3", "CTCF", "DNase"))

pls.p <- ggplot(model.res %>% filter(cCRE == "PLS"), aes(x = Feature, y = tissue)) +
  geom_point(aes(color = Estimate, size = -log10(PValue)), shape = 19) +
  scale_color_gradientn(colors = brewer.pal(11, "RdYlBu")[11:1], name = "Estimate") +
  scale_size_continuous(name = "Significance\n(-log10 P-value)", range = c(1, 10)) +
  theme_minimal() +
  labs(title = "PLS", x = "", y = "") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_point(data = model.res[model.res$Significant, ], aes(x = Feature, y = tissue), shape = 4, size = 3, color = "grey")

pels.p <- ggplot(model.res %>% filter(cCRE == "pELS"), aes(x = Feature, y = tissue)) +
  geom_point(aes(color = Estimate, size = -log10(PValue)), shape = 19) +
  scale_color_gradientn(colors = brewer.pal(11, "RdYlBu")[11:1], name = "Estimate") +
  scale_size_continuous(name = "Significance\n(-log10 P-value)", range = c(1, 10)) +
  theme_minimal() +
  labs(title = "pELS", x = "", y = "") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_point(data = model.res[model.res$Significant, ], aes(x = Feature, y = tissue), shape = 4, size = 4, color = "grey")

dels.p <- ggplot(model.res %>% filter(cCRE == "dELS"), aes(x = Feature, y = tissue)) +
  geom_point(aes(color = Estimate, size = -log10(PValue)), shape = 19) +
  scale_color_gradientn(colors = brewer.pal(11, "RdYlBu")[11:1], name = "Estimate") +
  scale_size_continuous(name = "Significance\n(-log10 P-value)", range = c(1, 10)) +
  theme_minimal() +
  labs(title = "dELS", x = "", y = "") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_point(data = model.res[model.res$Significant, ], aes(x = Feature, y = tissue), shape = 4, size = 4, color = "grey")

ggsave("../figure/unmethy_model.pdf", grid.arrange(pls.p, pels.p, dels.p, ncol = 3), width = 15, height = 8)
