#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1]]
}

sample_id <- get_arg("--sample")
h5_file   <- get_arg("--h5")
outdir    <- get_arg("--outdir", "results")
use_doubletprop <- tolower(get_arg("--doubletprop", "false")) %in% c("true", "1", "yes")

mat <- Read10X_h5(h5_file)

seu <- CreateSeuratObject(mat[["Gene Expression"]])
seu[["ADT"]] <- CreateAssay5Object(mat[["Antibody Capture"]])

DefaultAssay(seu) <- "RNA"
seu$orig.barcode <- colnames(seu)

seu <- RenameCells(seu, add.cell.id = sample_id)
seu <- subset(seu, subset = nFeature_RNA > 10)

# RNA / MT QC
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
mt_thresholds <- median(seu$percent.mt) + c(2, 3, 4) * mad(seu$percent.mt)
seu$qc_mt_high_2mad <- seu$percent.mt > mt_thresholds[1]
seu$qc_mt_high_3mad <- seu$percent.mt > mt_thresholds[2]
seu$qc_mt_high_4mad <- seu$percent.mt > mt_thresholds[3]

seu$qc_nFeature_RNA_low_300 <- seu$nFeature_RNA < 300
seu$qc_nFeature_RNA_high_4000 <- seu$nFeature_RNA > 4000

# ADT QC
adt_counts <- GetAssayData(seu, assay = "ADT", layer = "counts")
adt_total <- Matrix::colSums(adt_counts)
seu$nCount_ADT_portfolio <- as.numeric(adt_total)

log10_adt <- log10(seu$nCount_ADT_portfolio + 1)
adt_med <- median(log10_adt)
adt_mad <- mad(log10_adt)

adt_low  <- adt_med - 2 * adt_mad
adt_high <- adt_med + 4 * adt_mad

seu$qc_nCount_ADT_low_2mad  <- log10_adt < adt_low
seu$qc_nCount_ADT_high_4mad <- log10_adt > adt_high

keep <- !seu$qc_mt_high_3mad &
        !seu$qc_nFeature_RNA_high_4000 &
        !seu$qc_nFeature_RNA_low_300 &
        !seu$qc_nCount_ADT_low_2mad &
        !seu$qc_nCount_ADT_high_4mad

seu <- subset(seu, cells = colnames(seu)[keep])

# Dimension reduction
seu <- SCTransform(seu, verbose = FALSE)
seu <- RunPCA(
  seu,
  assay = "SCT",
  reduction.name = "pca.sct",
  verbose = FALSE
)
seu <- FindNeighbors(
  seu,
  reduction = "pca.sct",
  dims = 1:30,
  graph.name = c("SCT_nn", "SCT_snn"),
  verbose = FALSE
)
seu <- RunUMAP(
  seu,
  reduction = "pca.sct",
  dims = 1:30,
  reduction.name = "umap",
  verbose = FALSE
)

#  DoubletProp
seu <- DoubletProp::doublet_prop(
  seu,
  p = 0.05,
  doublet_col = "freemuxlet_DROPLET.TYPE",
  donor_col = "freemuxlet_SNG.BEST.GUESS",
  graph.name = "SCT_snn"
)

seu <- subset(seu, subset = doublet_prop_prediction == "SNG")

rds_file <- file.path(outdir, sprintf("%s.qced_seurat.rds", sample_id))
barcode_file <- file.path(outdir, sprintf("%s.qced_barcodes.txt", sample_id))

saveRDS(seu, file = rds_file)

barcodes <- data.frame(barcode = Cells(seu))
fwrite(barcodes, file = barcode_file, col.names = FALSE)
