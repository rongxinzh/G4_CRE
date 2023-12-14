
library(dplyr)
library(data.table)
library(stringr)
library(Biostrings)
source("./G4Hunter.R")

chr.seq.path <- "../data/ref/hg38.fa.gz"
chr.used <- paste0("chr", c(1:22, "X"))
window <- 25
score <- 1.5

seq <- readBStringSet(chr.seq.path, format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)
seq <- data.frame(seq = as.vector(as.character(seq)), chr = str_split(names(seq), (" "), simplify = TRUE)[, 1])
seq <- seq %>% filter(chr %in% chr.used)

#unlink(paste0("../data/G4/G4Hunter_w", window, "_s", score,"_hg38.txt"))

for (i in 1:dim(seq)[1]) {
  
  message("Loading sequence")
  tmp.seq <- seq[i, 1]
  tmp.chr <- seq[i ,2]
  message(paste0("Start predict: ", tmp.chr))

  tmp.predicted <- modG4huntref(k = window, hl = score, DNAStringSet(tmp.seq)[[1]], seqname = tmp.chr, with.seq = T, Gseq.only = F)
  tmp.predicted <- as.data.frame(tmp.predicted) %>% select(seqnames:k)

  fwrite(tmp.predicted, paste0("../data/G4/G4Hunter_w", window, "_s", score,"_hg38.txt"), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)

  message(paste0("Finish: ", tmp.chr))
  
}
