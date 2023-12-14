
set.seed(1)
options(scipen=200)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(bedtoolsr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(random))

GetRuns <- function(loc = NULL, mode = "G") {
  
  rand.f <- randomStrings(n = 1, len = 20) %>% as.character
  unlink(rand.f)
  loc[, 2] <- loc[, 2] - 1
  loc.seq <- bt.getfasta(fi = "../data/ref/hg38.fa", bed = loc, tab = TRUE)
  for (i in 1:dim(loc.seq)[1]) {
    tmp.fa <- loc.seq[i, 2]
    #tmp.runs.loc <- str_locate_all(toupper(tmp.fa), "(C{3,}|G{3,})") %>% data.frame
    if (mode == "G") {
      tmp.runs.loc <- str_locate_all(toupper(tmp.fa), "G{3,}") %>% data.frame
    } else {
      tmp.runs.loc <- str_locate_all(toupper(tmp.fa), "C{3,}") %>% data.frame
    }
    
    if (dim(tmp.runs.loc)[1] != 0) {
      tmp.runs.loc[, 1] <- tmp.runs.loc[, 1] + loc[i, 2]
      tmp.runs.loc[, 2] <- tmp.runs.loc[, 2] + loc[i, 2]
      fwrite(data.frame(chr = loc[i, 1], tmp.runs.loc), rand.f, 
        sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
    }

    tmp.runs.loc <- NULL
  }

  G.runs <- fread(rand.f, sep = "\t", header = FALSE) %>% data.frame
  unlink(rand.f)
  return(G.runs)
}

G4 <- fread("../data/G4/cCRE_G4_annotation.txt", sep = "\t", header = TRUE) %>% data.frame()
ccre <- fread("../data/ccre/V3/G4_cCRE_annotation.txt", sep = "\t", header = TRUE) %>% data.frame() %>% filter(chr %in% paste0("chr", c(1:22, "X")))

G4.ccre <- G4 %>% filter(overlap_count > 0)
G4.GRuns <- GetRuns(G4.ccre[G4.ccre[, "strand"] == "+", 1:3], "G") %>% unique()
G4.CRuns <- GetRuns(G4.ccre[G4.ccre[, "strand"] == "-", 1:3], "C") %>% unique()

ccre.GRuns <- GetRuns(ccre[, 1:3], "G")
ccre.CRuns <- GetRuns(ccre[, 1:3], "C")

# G4 Runs in cCREs
G4.GRuns.ccre <- bt.intersect(a = G4.GRuns, b = ccre, wa = TRUE, wb = TRUE) %>% unique() %>% data.frame()
G4.GRuns.ccre <- G4.GRuns.ccre[, c("V1", "V2", "V3", "V8", "V10")] %>% unique()
G4.CRuns.ccre <- bt.intersect(a = G4.CRuns, b = ccre, wa = TRUE, wb = TRUE) %>% unique() %>% data.frame()
G4.CRuns.ccre <- G4.CRuns.ccre[, c("V1", "V2", "V3", "V8", "V10")] %>% unique()

# Other Runs in cCREs
ccre.GRuns <- bt.intersect(a = ccre.GRuns, b = ccre, wa = TRUE, wb = TRUE) %>% unique() %>% data.frame()
ccre.GRuns <- bt.intersect(a = ccre.GRuns, b = G4.GRuns.ccre, wa = TRUE, v = TRUE) %>% unique() %>% data.frame()
ccre.GRuns <- ccre.GRuns[, c("V1", "V2", "V3", "V8", "V10")]
ccre.CRuns <- bt.intersect(a = ccre.CRuns, b = ccre, wa = TRUE, wb = TRUE) %>% unique() %>% data.frame()
ccre.CRuns <- bt.intersect(a = ccre.CRuns, b = G4.CRuns.ccre, wa = TRUE, v = TRUE) %>% unique() %>% data.frame()
ccre.CRuns <- ccre.CRuns[, c("V1", "V2", "V3", "V8", "V10")]

colnames(G4.GRuns.ccre) <- colnames(G4.CRuns.ccre) <- colnames(ccre.GRuns) <- colnames(ccre.CRuns) <- c("chr", "start", "end", "ID", "encodeLabel")

G4.Runs <- bind_rows(data.frame(G4.GRuns.ccre, strand = "+"), data.frame(G4.CRuns.ccre, strand = "-"))
ccre.nonG4.Runs <- bind_rows(data.frame(ccre.GRuns, strand = "+"), data.frame(ccre.CRuns, strand = "-"))
fwrite(G4.Runs, "../data/G4/G4_Runs.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
fwrite(ccre.nonG4.Runs, "../data/G4/nonG4_Runs_incCRE.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

