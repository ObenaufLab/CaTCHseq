.libPaths("")

library(optparse)

option_list <- list(
    make_option(c("--sample"),
        type = "character", default = NULL,
        help = "sample name"
    ),
    make_option(c("--data10X"),
        type = "character", default = NULL,
        help = "path to the 10X data"
    ),
    make_option(c("--catchBC"),
        type = "character", default = NULL,
        help = "path to the file with CaTCH barcodes"
    ),
    make_option(c("--max_mt"),
        type = "numeric", default = 10,
        help = "maximum percent of mitochondrial reads in a valid cell"
    ),
    make_option(c("--min_features"),
        type = "numeric", default = 500,
        help = "minimum number of detected features in a valid cell"
    ),
    make_option(c("--hvg_cutoff"),
        type = "numeric", default = 0.1,
        help = "cutoff value to call percentage of high variable genes (must be between 0 and 1)"
    ),
    make_option(c("--out"),
        type = "character", default = NULL,
        help = "path to the output file"
    ),
    make_option(c("--marker"),
        type = "character", default = "/tools/data/R/stagemarkers_xue2020.rds",
        help = "path to cellcycle marker file in RDS format"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

print(paste0("CLI: ", opt))

if (is.null(opt$sample) || is.null(opt$out)) {
    print_help(opt_parser)
    stop("Not enough parameters")
}


#### Source Functions ####
source("../R_collection/singlecell_utils.R")

############################

library(tidyverse)
library(scater)
library(scran)
library(SingleCellExperiment)
library(Seurat)
# options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")
library(SeuratWrappers)
library(sctransform)
library(R.filesets)

### Load all experiments, add CaTCH barcodes as layer and converting to Seuratv5 Object
sce <- create_SCEs(opt$sample, opt$data10X, opt$catchBC)

### Run QC
# calculate percentage of MT reads
sce <- PercentageFeatureSet(sce, pattern = "^MT-|^mt-", col.name = "percent.mt")
sce[["is.low_yield"]] <- sce@meta.data %>%
    pull(nFeature_RNA) %>%
    {
        case_when(. < opt$min_features ~ TRUE, .default = FALSE)
    }
sce[["is.damaged"]] <- sce@meta.data %>%
    pull(percent.mt) %>%
    {
        case_when(. > opt$max_mt ~ TRUE, .default = FALSE)
    }

pdf(file = paste0(opt$out, "_QC_Violin_MT_content.pdf"), width = 30, height = 10)
VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "Sample")
dev.off()
#
pdf(file = paste0(opt$out, "_QC_Scatter_MT_content.pdf"), width = 30, height = 10)
FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Sample")
dev.off()
#
pdf(file = paste0(opt$out, "_QC_Scatter_Feature_content.pdf"), width = 30, height = 10)
FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Sample")
dev.off()

### Filter for MT content and min reads
sce <- subset(sce, subset = nFeature_RNA > opt$min_features & percent.mt < opt$max_mt)

### Split SCE for integrative analysis
### ALREADY DONE BY MERGE
# sce <- split_SCE(sce)

### Normalize and scale counts
sce <- normalize_SCE(sce)

### SCtransform counts
sce <- sctransform_SCE(sce)

### Run initial dim reduction for norm
sce <- reduceDims_SCE(sce)

### Run initial dim reduction for STC
sce <- reduceDims_SCE(sce, assay = "SCT", reduction.name = "pca_sct")

### Integrate the SCE layers for integrative analysis
sce <- integrate_SCE(sce, assay = "RNA", orig.reduction = "pca", new.reduction = "integrated.cca", normalization.method = "LogNormalize")

### Integrate the SCE layers after SCT for integrative analysis
sce <- integrate_SCE(sce, assay = "SCT", orig.reduction = "pca_sct", new.reduction = "sct_integrated.cca", normalization.method = "SCT")

### Cluster integrated layers for plotting

sce <- cluster_SCE(sce, assay = "RNA_integrated.cca", reduction = "integrated.cca", cluster.name = "integrated.cca_cluster")
sce <- cluster_SCE(sce, assay = "RNA", reduction = "pca", cluster.name = "pca_cluster")

sce <- cluster_SCE(sce, assay = "SCT", reduction = "pca_sct", cluster.name = "sct_cluster")
sce <- cluster_SCE(sce, assay = "SCT_integrated", reduction = "pca_sct", cluster.name = "sct_integrated_cluster")

## Run UMAPs
sce <- umap_SCE(sce, assay = "RNA_integrated.cca", reduction = "integrated.cca", dims = 1:30, reduction.name = "umap_integrated.cca", n.neighbors = 30L, min.dist = 0.1, spread = 5)
sce <- umap_SCE(sce, assay = "RNA", reduction = "pca", dims = 1:30, reduction.name = "umap_pca", n.neighbors = 30L, min.dist = 0.1, spread = 5)

sce <- umap_SCE(sce, assay = "SCT_integrated", reduction = "sct_integrated.cca", dims = 1:10, reduction.name = "umap_integrated.cca", n.neighbors = 30L, min.dist = 0.01, spread = 5)
sce <- umap_SCE(sce, assay = "SCT", reduction = "pca_sct", dims = 1:30, reduction.name = "umap_pca_sct", n.neighbors = 30L, min.dist = 0.1, spread = 5)
### Assign cell stage and categories

print(paste0("Loading cell stage markers from ", opt$marker, sep = ""))
markerfile <- loadRDS(opt$marker)

print("Attempting to assign cell stage...")
s.genes <- markerfile$S
g2m.genes <- markerfile$G2M
tryCatch(
    {
        sce <- CellCycleScoring(sce, assay = "RNA_integrated.cca", s.features = toupper(s.genes), g2m.features = toupper(g2m.genes), set.ident = FALSE)
    },
    error = function(e) {
        print(paste("ERROR COUGHT:  ", e, " WILL SKIP ASSIGNMENT OF CELL STAGE AND SET TO DEFAULT G0"))
        sce@meta.data$Phase <- "G0"
        sce@meta.data$S.Score <- 0
        sce@meta.data$G2M.Score <- 0
    }
)

print("Categorize the cells ...")
cell.categories <- c("Good", "Damaged", "Few features", "Damaged AND few features")

sce@meta.data$Category <- sce@meta.data %>%
    dplyr::mutate(tmp = as.integer(is.damaged) * 1 + as.integer(is.low_yield) * 2) %>%
    dplyr::select("tmp") %>%
    purrr::map(.f = ~ cell.categories[.x + 1]) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(Class = factor(tmp, levels = cell.categories)) %>%
    dplyr::select(Class)

saveRDS(sce, file = paste0(opt$out, "_SCE.rds.gz"), compress = "gzip")