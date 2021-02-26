library(tidyverse)
library(ape)
library(phangorn)
library(ggtree)
library(plotly)
KnownSeqs <- read.dna("USA_WI_MinK.aligned.fasta", format = "fasta")

msa <- phyDat(KnownSeqs, type = "DNA")
dm  <- dist.ml(msa)
treeNJ <- NJ(dm)
phylocanvas(treeNJ, treetype = "radial", alignlabels = T)

t <- ggtree(treeNJ)
ggplotly(t)

