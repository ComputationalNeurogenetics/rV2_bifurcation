## ----Libraries----------------------------------------------------------------
library(Seurat)
library(Signac)
library(tidyverse)
library(GenomicRanges)
library(qs)
source("~/Workspace/generic code/AuxFunctions.R")


## ----Loading data-------------------------------------------------------------
s.data <- qread("~/Workspace/NeuronalFeatureSpace/scATAC_data/E14_DownstreamReady.300922.RNA.pos.idents.qs", nthreads = 6)
DefaultAssay(s.data) <- "peaks"


## -----------------------------------------------------------------------------
# Picking P2R and P3 precursor clusters based on the cluster annotation made from marker data
#ioi.di <- c("P3 GABAergic neuron precursors","P2 (P2r GABA, P2c Glut)")
# These correspond to gaba_8, gaba_10
# ioi.di <- c("P3 GABAergic neuron precursors","P2 (P2r GABA, P2c Glut)")

#ioi.mb <- c("GABA, Glut precursors (mixed Lhx4 or Otp)", "GABAergic precurors (Tal2, Nr2f1, Sst)")

#ioi.r1 <- c("rV2 GABA precursors","DL GABAergic precursors (Skor1, Pax2)")


## ----Writing E14DI barcodes out for TOBIAS calculations-----------------------
all.barcodes <- tibble(barcodes=Cells(s.data), seurat_clusters=Idents(s.data))
#di.pre.gaba <- all.barcodes %>%  filter(seurat_clusters %in% ioi.di)
#di.pre.gaba$barcodes <- str_remove(di.pre.gaba$barcodes,pattern = "_._")
write_tsv(all.barcodes, file = "/Volumes/ExtSSD/data/E14_BAM_merge/E14.barcodes.clusters", col_names = FALSE)

