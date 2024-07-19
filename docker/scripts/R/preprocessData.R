.libPaths("")

library(optparse)

option_list <- list(
    make_option(c("--sample"),
        type = "character", default = NULL,
        help = "sample name"
    ),
    make_option(c("--baseCond"),
        type = "character", default = NULL,
        help = "Base condition for DE Analysis"
    ),
    make_option(c("--data10X"),
        type = "character", default = NULL,
        help = "path to the 10X data"
    ),
    make_option(c("--catchBC"),
        type = "character", default = NULL,
        help = "path to the file with CaTCH barcodes"
    ),
    make_option(c("--annotation"),
        type = "character", default = NULL,
        help = "path to the matching annotation file in GTF format"
    ),
    make_option(c("--minBC"),
        type = "numeric", default = 10,
        help = "minimum number of barcode reads per cell"
    ),
    make_option(c("--singletCut"),
        type = "numeric", default = 0.9,
        help = "minimum ratio of barcode 1 per cell for category 'Singlet'"
    ),
    make_option(c("--bc1Cut"),
        type = "numeric", default = 0.4,
        help = "minimum ratio of barcode 1 reads per cell for category 'Double_Integration'"
    ),
    make_option(c("--bc2Cut"),
        type = "numeric", default = 0.3,
        help = "minimum ratio of barcode 2 reads per cell for category 'Double_Integration'"
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
    ),
    make_option(c("--libpath"),
        type = "character", default = "/tools/scripts/R/",
        help = "path to R libs, trailing slash is required"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

print(paste0("CLI: ", opt))

if (is.null(opt$sample) || is.null(opt$out) || is.null(opt$annotation)) {
    print_help(opt_parser)
    stop("Not enough parameters")
}

#### Source Functions ####
source(paste0(opt$libpath, "singlecell_utils.R"))

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
sl <- create_SCEs(opt$sample, opt$data10X, opt$catchBC, opt$annotation, opt$minBC, opt$singletCut, opt$bc1Cut, opt$bc2Cut)

sce <- sl$sce
seurat_sce <- sl$seurat_sce

print(paste("SCE object contains ", ncol(sce), " cells and ", nrow(sce), " genes.\n Seurat object contains ", ncol(seurat_sce), " cells and ", nrow(seurat_sce), " genes."))
# seurat_sce <- DietSeurat(seurat_sce, layers = "counts")  # Get rid of artificial data slot

rm(sl) # clean up

#### Annotate Samples, Conditions and Replicates ####
sce$Condition <- as.factor(unlist(lapply(sce$Sample, function(x) unlist(str_split(x, "_"))[2])))
sce$Replicate <- as.factor(unlist(lapply(sce$Sample, function(x) paste(unlist(str_split(x, "_"))[2], unlist(str_split(x, "_"))[3], sep = "_"))))
sce$Sample <- as.factor(unlist(lapply(sce$Sample, function(x) unlist(str_split(x, "_"))[1])))

seurat_sce$Sample.orig <- as.factor(seurat_sce$Sample)
seurat_sce$Condition <- as.factor(unlist(lapply(seurat_sce$Sample, function(x) unlist(str_split(x, "_"))[2])))
seurat_sce$Replicate <- as.factor(unlist(lapply(seurat_sce$Sample, function(x) paste(unlist(str_split(x, "_"))[2], unlist(str_split(x, "_"))[3], sep = "_"))))
seurat_sce$Sample <- as.factor(unlist(lapply(seurat_sce$Sample, function(x) unlist(str_split(x, "_"))[1])))

#### calculate percentage of MT reads SEURAT ####
seurat_sce <- PercentageFeatureSet(seurat_sce, pattern = "^MT-|^mt-", col.name = "percent.mt")
pdf(file = paste0(opt$out, "_QC_Violin_MT_content.pdf"), width = 30, height = 10)
VlnPlot(seurat_sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "Sample")
dev.off()
#
pdf(file = paste0(opt$out, "_QC_Scatter_MT_content.pdf"), width = 30, height = 10)
FeatureScatter(seurat_sce, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Sample")
dev.off()
#
pdf(file = paste0(opt$out, "_QC_Scatter_Feature_content.pdf"), width = 30, height = 10)
FeatureScatter(seurat_sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Sample")
dev.off()

#### normalize counts and calculate percentage of MT reads SCE ####
print("Normalize sce ...")

is.mitochondrial <- grepl(
    pattern = "^MT-|^mt-", # Needed e.g. for mouse 10x data with cellranger prebuilt index
    x = rownames(sce),
    ignore.case = FALSE,
    perl = TRUE
)

sce <- sce %>%
    scater::logNormCounts() %>%
    scater::addPerCellQC(subsets = list(MT = is.mitochondrial)) %>%
    scater::addPerFeatureQC()

### Normalize and scale counts Seurat ####
print("Normalize Seurat ...")

seurat_sce <- normalize_Seurat(seurat_sce)

#### annotate low yield and damaged cells ####
print("Annotate sce ...")
rowData(sce)["is.mitochondrial"] <- is.mitochondrial
colData(sce)["is.damaged"] <- colData(sce)[, "subsets_MT_percent"] > opt$max_mt
colData(sce)["is.low_yield"] <- colData(sce)[, "detected"] < opt$min_features

print("Annotate seurat ...")
seurat_sce[["is.damaged"]] <- seurat_sce@meta.data %>%
    pull(percent.mt) %>%
    {
        case_when(. > opt$max_mt ~ TRUE, .default = FALSE)
    }
seurat_sce[["is.low_yield"]] <- seurat_sce@meta.data %>%
    pull(nFeature_RNA) %>%
    {
        case_when(. < opt$min_features ~ TRUE, .default = FALSE)
    }

#### Categorize the cells ####
cell.categories <- c("Good", "Damaged", "Few features", "Damaged AND few features")

print("Categorize sce ...")
colData(sce)["Category"] <- colData(sce) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(tmp = as.integer(is.damaged) * 1 + as.integer(is.low_yield) * 2) %>%
    dplyr::select("tmp") %>%
    purrr::map(.f = ~ cell.categories[.x + 1]) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(Class = factor(tmp, levels = cell.categories)) %>%
    dplyr::select(Class)

print("Categorize seurat ...")

seurat_sce@meta.data$Category <- seurat_sce@meta.data %>%
    dplyr::mutate(tmp = as.integer(is.damaged) * 1 + as.integer(is.low_yield) * 2) %>%
    dplyr::select("tmp") %>%
    purrr::map(.f = ~ cell.categories[.x + 1]) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(Class = factor(tmp, levels = cell.categories)) %>%
    dplyr::pull(Class)

### Assign cell stage and categories
print(paste0("Loading cell stage markers from ", opt$marker, sep = ""))
markerfile <- loadRDS(opt$marker)

#### Joining seurat layers for downstream annotation ####
print("Joining seurat layers ...")
seurat_sce <- join_Seurat(seurat_sce, assay = "RNA", layers = "data", new = "data")
seurat_sce <- join_Seurat(seurat_sce, assay = "RNA", layers = "counts", new = "counts")

#### Attempting to assign cell phase via Seurat ####
print("Phase annotation seurat ...")

s.genes <- markerfile$S
g2m.genes <- markerfile$G2M

withCallingHandlers(
    {
        seurat_sce <- CellCycleScoring(seurat_sce, assay = "RNA", slot = "data", s.features = str_to_upper(s.genes), g2m.features = str_to_upper(g2m.genes), set.ident = FALSE)
    },
    warning = function(w) {
        if (grepl("Could not find enough features in the object", w$message)) {
            print("WARNING COUGHT:  Could not find enough features in the object. Will try to match case")
            seurat_sce <- CellCycleScoring(seurat_sce, assay = "RNA", slot = "data", s.features = str_to_title(s.genes), g2m.features = str_to_title(g2m.genes), set.ident = FALSE)
        } else {
            message(w$message)
        }
    },
    error = function(e) {
        print(paste("ERROR COUGHT:  ", e, " WILL SKIP ASSIGNMENT OF SEURAT CELL PHASE AND SET TO DEFAULT G0"))
        seurat_sce@meta.data$Phase <- "G0"
        seurat_sce@meta.data$S.Score <- 0
        seurat_sce@meta.data$G2M.Score <- 0
    }
)

#### Attempting to assign cell stage to SCEs ####
print("Cellstage sce ...")
sce <- sce %>% assignCategoryByMarker(markers = markerfile, col.name = "CellStage")
print("Cellstage seurat ...")
seurat_sce <- seurat_sce %>% assignCategoryByMarker(markers = markerfile, col.name = "CellStage")

### Split SCE again
seurat_sce <- split_Seurat(seurat_sce, by = seurat_sce$Sample.orig)

#### Save unfiltered sce and seurat_sce objects ####
print("Save unfiltered ...")
saveRDS(sce, file = paste0(opt$out, "_unfiltered_sce.rds.gz"), compress = "gzip")
saveRDS(seurat_sce, file = paste0(opt$out, "_unfiltered_seurat_sce.rds.gz"), compress = "gzip")

### Filter for MT content and min reads
print("Filter seurat ...")
# seurat_sce <- subset(seurat_sce, subset = nFeature_RNA > opt$min_features & percent.mt < opt$max_mt)
seurat_sce <- subset(seurat_sce, subset = is.low_yield == FALSE & is.damaged == FALSE)
print("Filter sce ...")
sce <- sce[, sce$is.low_yield == FALSE & sce$is.damaged == FALSE]
print(paste0("Keeping ", ncol(sce), " Cells from SCE object and ", ncol(seurat_sce), " Cells from Seurat object."))


# Plot distribution of CaTCH barcode ratios
toplot <- seurat_sce[[]] %>%
    group_by(CaTCH.BCs) %>%
    filter(n() > 1)

a <- ggplot(toplot, aes(x = Condition, y = CaTCH.BC1 / CaTCH.Sum)) +
    geom_violin() +
    scale_y_log10(guide = "axis_logticks") +
    theme_bw()

b <- ggplot(toplot, aes(x = Condition, y = CaTCH.BC2 / CaTCH.Sum)) +
    geom_violin() +
    scale_y_log10(guide = "axis_logticks") +
    theme_bw()

c <- ggplot(toplot, aes(x = Condition, y = (CaTCH.BC1 + CaTCH.BC2) / CaTCH.Sum)) +
    geom_violin() +
    scale_y_log10(guide = "axis_logticks") +
    theme_bw()

pdf(file = paste0(opt$out, "_QC_Barcode_Ratio_Distribution.pdf"), width = 10, height = 30)
print(a | b | c)
dev.off()

rm(a, b, c, toplot)

#### Run PCA and UMAP ####
print("Identify the top variable genes...")
gene.var <- modelGeneVar(sce)
hvg <- getTopHVGs(stats = gene.var, prop = opt$hvg_cutoff)

print("Run the PCA ...")
sce <- runPCA(sce, subset_row = hvg)
set.seed(42)

print("Run tSNE and UMAP analyses ...")
sce <- runTSNE(sce, dimred = "PCA")
sce <- runUMAP(sce, dimred = "PCA")

print("Clustering ...")
g <- buildSNNGraph(sce, use.dimred = "PCA")
cluster <- igraph::cluster_walktrap(g)$membership
colData(sce)["Cluster"] <- factor(cluster)


### SCtransform counts
print("SCtransform Seurat ...")
seurat_sce <- sctransform_Seurat(seurat_sce)
### Run initial dim reduction for norm
print("PCA Seurat ...")
seurat_sce <- reduceDims_Seurat(seurat_sce)
### Run initial dim reduction for STC
seurat_sce <- reduceDims_Seurat(seurat_sce, assay = "SCT", reduction.name = "pca_sct")
### Cluster
print("Clustering Seurat ...")
seurat_sce <- cluster_Seurat(seurat_sce, assay = "RNA", reduction = "pca", cluster.name = "pca_cluster")
seurat_sce <- cluster_Seurat(seurat_sce, assay = "SCT", reduction = "pca_sct", cluster.name = "sct_cluster")
## Run UMAPs
print("UMAP Seurat ...")
seurat_sce <- umap_Seurat(seurat_sce, assay = "RNA", reduction = "pca", dims = 1:30, reduction.name = "umap_pca", n.neighbors = 30L, min.dist = 0.1, spread = 5)
seurat_sce <- umap_Seurat(seurat_sce, assay = "SCT", reduction = "pca_sct", dims = 1:30, reduction.name = "umap_pca_sct", n.neighbors = 30L, min.dist = 0.1, spread = 5)

#### Assign CaTCH barcode indices based on their abundance in the reference samples ####
print("Assign CaTCH barcodes ...")

# Create unique BC for merge
colData(sce)["CaTCH.BC_unique"] <- colData(sce) %>%
  as_tibble() %>%
  select(CaTCH.Status, CaTCH.BCs) %>%
  rowwise() %>%
  mutate(CaTCH.BC_unique = ifelse(CaTCH.Status == "Singlet", str_split_1(CaTCH.BCs, ";")[1], ifelse(CaTCH.Status == "Double_Integration", paste0(str_split_1(CaTCH.BCs, ";")[1], str_split_1(CaTCH.BCs, ";")[2], sep="+"), paste0(CaTCH.BCs)))) %>%
  ungroup() %>%
  select(CaTCH.BC_unique)

seurat_sce@meta.data$CaTCH.BC_unique <- seurat_sce@meta.data %>%
  select(CaTCH.Status, CaTCH.BCs) %>%
  rowwise() %>%
  mutate(CaTCH.BC_unique = ifelse(CaTCH.Status == "Singlet", str_split_1(CaTCH.BCs, ";")[1], ifelse(CaTCH.Status == "Double_Integration", paste(str_split_1(CaTCH.BCs, ";")[1], str_split_1(CaTCH.BCs, ";")[2], sep="+"), CaTCH.BCs))) %>%
  ungroup() %>%
  pull(CaTCH.BC_unique)

tmp <- seurat_sce@meta.data %>%
    as_tibble() %>%
    dplyr::select(c(CaTCH.Status, Condition, CaTCH.BC_unique, Sample)) %>%
    filter(CaTCH.Status == "Singlet" | CaTCH.Status == "Double_Integration", Condition == opt$baseCond) %>%
    dplyr::select(CaTCH.BC_unique, Sample, CaTCH.Status) %>%
    group_by(CaTCH.BC_unique, Sample) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    distinct() %>%
    pivot_wider(names_from = Sample, values_from = n, values_fill = 0) %>%
    filter(rowSums(across(starts_with(opt$baseCond))) > 0) %>%
    mutate(.Means = rowMeans(across(starts_with(opt$baseCond)))) %>%
    arrange(by = desc(.Means)) %>%
    rowid_to_column(".ID") %>%
    mutate(CaTCH.BC_ID = ifelse(is.na(.ID), "BC_0", ifelse(CaTCH.Status == "Singlet", paste0("BC_", .ID), paste0("BC*_", .ID)))) %>%
    dplyr::select(-.Means, -.ID) %>%
    relocate(CaTCH.BC_ID, .after = CaTCH.BC_unique) %>%
    dplyr::select(CaTCH.BC_unique, CaTCH.BC_ID)


colData(sce)["CaTCH.BC_ID"] <- colData(sce) %>%
    as_tibble() %>%
    left_join(y = tmp, by = "CaTCH.BC_unique") %>%
    mutate(CaTCH.BC_ID = ifelse(is.na(CaTCH.BC_ID), "BC_0", CaTCH.BC_ID)) %>%
    mutate(CaTCH.BC_ID = factor(CaTCH.BC_ID, levels = str_sort(unique(CaTCH.BC_ID), numeric = TRUE))) %>%
    dplyr::select(CaTCH.BC_ID)

seurat_sce@meta.data$CaTCH.BC_ID <- seurat_sce@meta.data %>%
    left_join(y = tmp, by = "CaTCH.BC_unique") %>%
    mutate(CaTCH.BC_ID = ifelse(is.na(CaTCH.BC_ID), "BC_0", CaTCH.BC_ID)) %>%
    mutate(CaTCH.BC_ID = factor(CaTCH.BC_ID, levels = str_sort(unique(CaTCH.BC_ID), numeric = TRUE))) %>%
    pull(CaTCH.BC_ID)

#### Save final objects ####
print("Final Save ...")

saveRDS(sce, file = paste0(opt$out, "_filtered_sce.rds.gz"), compress = "gzip")
saveRDS(seurat_sce, file = paste0(opt$out, "_filtered_seurat_sce.rds.gz"), compress = "gzip")
