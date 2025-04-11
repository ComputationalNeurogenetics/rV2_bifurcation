# Standalone script to export tracks to IGV-compatible formats

# Load required libraries
library(tidyverse)
library(Seurat)
library(Signac)
library(qs)
library(rtracklayer)
library(gUtils)
library(GenomicRanges)

# Set parameters
cores <- 6


# Load the rV2.data object
rV2.data <- qread("../scATAC_data/nmm_rV2_subset_relabeled_031023_links.qs", nthreads = cores)
con <- dbConnect(RSQLite::SQLite(), "../../TOBIAS.dr.h12_191124.sqlite") 

# Load gene annotations
annotations <- readRDS("../metadata/Ensembel_Mmusculus79_annotation.Rds")
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"
Signac::Annotation(rV2.data) <- annotations

# Define output directory for IGV-compatible files
igv_output_dir <- "~/OneDrive - University of Helsinki/E12R1 project/Manuscript I - Regulation of Tal1 dependent rV2 lineage bifurcation/Figures/IGV_tracks/"
if (!dir.exists(igv_output_dir)) {
  dir.create(igv_output_dir, recursive = TRUE)
}

# Load conservation data
cons.gr <- qread("../metadata/cons.filt.gr.qs", nthreads = cores)
cons.gr.reduced <- GenomicRanges::reduce(cons.gr)

# Export Conservation Track
export(cons.gr.reduced, con = file.path(igv_output_dir, "Conservation_track.bed"), format = "BED")

# Fetch links Lmo1
region.of.interest <- "chr7-109088977-109219982"

links <- fetchGRangesLinks(
  "../../TOBIAS.dr.h12_191124.sqlite",
  zscore_threshold = 2,
  pvalue_threshold = 0.05,
  coordinate_filter = region.of.interest,
  target.gene_name ="Lmo1"
)

# Export Link Track
export(StringToGRanges(links$peak), con = file.path(igv_output_dir, "Link_Lmo1_track.bed"), format = "BED")
write.table(as_tibble(links), file = file.path(igv_output_dir, "Link_Lmo1_track_num.csv"), sep = ",")

# Fetch links Ebf2
region.of.interest <- "chr14-67183957-67478090"
links <- fetchGRangesLinks(
  "../../TOBIAS.dr.h12_191124.sqlite",
  zscore_threshold = 2,
  pvalue_threshold = 0.05,
  coordinate_filter = region.of.interest,
  target.gene_name ="Ebf2"
)

# Export Link Track
export(StringToGRanges(links$peak), con = file.path(igv_output_dir, "Link_Ebf2_track.bed"), format = "BED")
write.table(as_tibble(links), file = file.path(igv_output_dir, "Link_Ebf2_track_num.csv"), sep = ",")

# Fetch links Plekha1
region.of.interest <- "chr7-130815984-130960831"
links <- fetchGRangesLinks(
  "../../TOBIAS.dr.h12_191124.sqlite",
  zscore_threshold = 2,
  pvalue_threshold = 0.05,
  coordinate_filter = region.of.interest,
  target.gene_name ="Plekha1"
)

# Export Link Track
export(StringToGRanges(links$peak), con = file.path(igv_output_dir, "Link_Plekha1_track.bed"), format = "BED")
write.table(as_tibble(links), file = file.path(igv_output_dir, "Link_Plekha1_track_num.csv"), sep = ",")

# Export Features Track
features <- rownames(rV2.data[["peaks"]])
export(StringToGRanges(features), con = file.path(igv_output_dir, "Features_track.bed"), format = "BED")

# Export footprints Tal1
query <- "SELECT t.seqnames,t.start,t.end, t.GA1_2_score FROM tobias as t WHERE t.TF_gene_name='TAL1' AND t.GA1_2_bound = 1"
Tal1.fp <- dbGetQuery(con, query)
Tal1.fp.gr <- makeGRangesFromDataFrame(Tal1.fp, keep.extra.columns = T)
export(Tal1.fp.gr, con = file.path(igv_output_dir, "Tal1.fp.bed"), format = "BED")
write.table(as_tibble(Tal1.fp.gr), file = file.path(igv_output_dir, "FP_Tal1_track_num.csv"), sep = ",")

# Export footprints Gata2
query <- "SELECT t.seqnames,t.start,t.end, t.GA1_2_score FROM tobias as t WHERE t.TF_gene_name='GATA2' AND t.GA1_2_bound = 1"
Gata2.fp <- dbGetQuery(con, query)
Gata2.fp.gr <- makeGRangesFromDataFrame(Gata2.fp, keep.extra.columns = T)
export(Gata2.fp.gr, con = file.path(igv_output_dir, "Gata2.fp.bed"), format = "BED")
write.table(as_tibble(Gata2.fp.gr), file = file.path(igv_output_dir, "FP_Gata2_track_num.csv"), sep = ",")

# Export footprints Gata3
query <- "SELECT t.seqnames,t.start,t.end, t.GA1_2_score FROM tobias as t WHERE t.TF_gene_name='GATA3' AND t.GA1_2_bound = 1"
Gata3.fp <- dbGetQuery(con, query)
Gata3.fp.gr <- makeGRangesFromDataFrame(Gata3.fp, keep.extra.columns = T)
export(Gata3.fp.gr, con = file.path(igv_output_dir, "Gata3.fp.bed"), format = "BED")
write.table(as_tibble(Gata3.fp.gr), file = file.path(igv_output_dir, "FP_Gata3_track_num.csv"), sep = ",")

# Print success message
cat("All tracks have been exported to IGV-compatible formats in:", igv_output_dir, "\n")
