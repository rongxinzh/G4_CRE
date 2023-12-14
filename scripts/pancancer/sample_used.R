
set.seed(1)
options(scipen=200)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))

pan.atac.files <- list.files("../data/pancATAC/TCGA-ATAC_Cancer_Type-specific_PeakCalls/")
pan.atac.ct <- str_split(pan.atac.files, "_", simplify = TRUE)[, 1]

superenh <- fread("../data/pancEnh/super_enhancer.txt", sep = "\t", header = TRUE) %>% data.frame() %>% filter(chr %in% paste0("chr", c(1:22, "X")))
superenh.flt <- superenh[ , c("chr", "start", "end", "enhancer_id", "peak_score")]
pan.enh.ct <- str_split(superenh.flt$enhancer_id, "_", simplify = TRUE)[, 1] %>% unique()

raw.prom <- fread("../data/pancExp/alternativePromoters_rawPromoterActivity.tsv", sep = "\t", header = FALSE) %>% data.frame()
sample <- fread("../data/pancExp/pcawg_sample_sheet.tsv", sep = "\t", header = TRUE) %>% data.frame()
raw.prom.sample <- raw.prom[1, 3:ncol(raw.prom)] %>% t()
raw.prom.sample <- data.frame(aliquot_id = c(raw.prom.sample)) %>% inner_join(sample, by = "aliquot_id")
raw.prom.sample <- raw.prom.sample[!(str_starts(raw.prom.sample$dcc_specimen_type, "Normal")),]
pan.prom.ct <- str_split(raw.prom.sample$dcc_project_code, "-", simplify = TRUE)[, 1] %>% unique()

fwrite(intersect(intersect(pan.atac.ct, pan.enh.ct), pan.prom.ct) %>% data.frame(), "../data/PANC_CT_used.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)