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
  make_option(c("--annotation"),
    type = "character", default = NULL,
    help = "path to the matching annotation file in GTF format"
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

if (is.null(opt$sample) || is.null(opt$out) || is.null(opt$annotation)) {
  print_help(opt_parser)
  stop("Not enough parameters")
}

#### Source Functions ####
source("singlecell_utils.R")

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
sl <- create_SCEs(opt$sample, opt$data10X, opt$catchBC, opt$annotation)

sce <- sl$sce
seurat_sce <- sl$seurat_sce
#seurat_sce <- DietSeurat(seurat_sce, layers = "counts")  # Get rid of artificial data slot

rm(sl) # clean up

#### Annotate Samples, Conditions and Replicates ####
sce$Condition <- as.factor(unlist(lapply(sce$Sample, function(x) unlist(str_split(x, "_"))[2])))
sce$Replicate <- as.factor(unlist(lapply(sce$Sample, function(x) paste(unlist(str_split(x, "_"))[2], unlist(str_split(x, "_"))[3], sep = "_"))))
sce$Sample <- as.factor(unlist(lapply(sce$Sample, function(x) unlist(str_split(x, "_"))[1])))

seurat_sce$Sample.orig <- as.factor(seurat_sce$Sample)
seurat_sce$Condition <- as.factor(unlist(lapply(seurat_sce$Sample, function(x) unlist(str_split(x, "_"))[2])))
seurat_sce$Replicate <- as.factor(unlist(lapply(seurat_sce$Sample, function(x) paste(unlist(str_split(x, "_"))[2], unlist(str_split(x, "_"))[3], sep = "_"))))
seurat_sce$Sample <- as.factor(unlist(lapply(seurat_sce$Sample, function(x) unlist(str_split(x, "_"))[1])))

ref.Condition <- opt$baseCond

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
print("Normalize the counts ...")
sce <- sce %>%
  scater::logNormCounts() %>%
  scater::addPerCellQC(subsets = list(MT = grepl(
    pattern = "^MT-",
    x = rownames(sce),
    ignore.case = TRUE # Needed e.g. for mouse 10x data with cellranger prebuilt index
  )
  )) %>%
  scater::addPerFeatureQC()

### Normalize and scale counts Seurat ####
seurat_sce <- normalize_Seurat(seurat_sce)

#### annotate low yield and damaged cells ####
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


is.mitochondrial <- grepl(
  pattern = "^MT-",
  x = rownames(sce),
  ignore.case = TRUE # Needed e.g. for mouse 10x data with cellranger prebuilt index
)

rowData(sce)["is.mitochondrial"] <- is.mitochondrial
colData(sce)["is.damaged"] <- colData(sce)[, "subsets_MT_percent"] > opt$max_mt
colData(sce)["is.low_yield"] <- colData(sce)[, "detected"] < opt$min_features

#### Categorize the cells ####
cell.categories <- c("Good", "Damaged", "Few features", "Damaged AND few features")

seurat_sce@meta.data$Category <- seurat_sce@meta.data %>%
  dplyr::mutate(tmp = as.integer(is.damaged) * 1 + as.integer(is.low_yield) * 2) %>%
  dplyr::select("tmp") %>%
  purrr::map(.f = ~ cell.categories[.x + 1]) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(Class = factor(tmp, levels = cell.categories)) %>%
  dplyr::pull(Class)


colData(sce)["Category"] <- colData(sce) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(tmp = as.integer(is.damaged) * 1 + as.integer(is.low_yield) * 2) %>%
  dplyr::select("tmp") %>%
  purrr::map(.f = ~ cell.categories[.x + 1]) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(Class = factor(tmp, levels = cell.categories)) %>%
  dplyr::select(Class)


### Assign cell stage and categories
print(paste0("Loading cell stage markers from ", opt$marker, sep = ""))
markerfile <- loadRDS(opt$marker)

#### Attempting to assign cell phase via Seurat ####
s.genes <- markerfile$S
g2m.genes <- markerfile$G2M

seurat_sce <- join_Seurat(seurat_sce, assay = "RNA", layers = "data", new = "data")
seurat_sce <- join_Seurat(seurat_sce, assay = "RNA", layers = "counts", new = "counts")

tryCatch(
  {
    seurat_sce <- CellCycleScoring(seurat_sce, assay = "RNA", slot = "data", s.features = toupper(s.genes), g2m.features = toupper(g2m.genes), set.ident = FALSE)
  },
  error = function(e) {
    print(paste("ERROR COUGHT:  ", e, " WILL SKIP ASSIGNMENT OF SEURAT CELL PHASE AND SET TO DEFAULT G0"))
    seurat_sce@meta.data$Phase <- "G0"
    seurat_sce@meta.data$S.Score <- 0
    seurat_sce@meta.data$G2M.Score <- 0
  }
)

#### Attempting to assign cell stage to SCEs ####
sce <- sce %>% assignCategoryByMarker(markers = markerfile, col.name = "CellStage")
seurat_sce <- seurat_sce %>% assignCategoryByMarker(markers = markerfile, col.name = "CellStage", obj.type = "SEURAT")

### Split SCE again
seurat_sce <- split_Seurat(seurat_sce, by = seurat_sce$Sample.orig)

#### Save unfiltered sce and seurat_sce objects ####
saveRDS(sce, file = paste0(opt$out, "_unfiltered_sce.rds"))
saveRDS(seurat_sce, file = paste0(opt$out, "_unfiltered_seurat_sce.rds"))

### Filter for MT content and min reads
#seurat_sce <- subset(seurat_sce, subset = nFeature_RNA > opt$min_features & percent.mt < opt$max_mt)
seurat_sce <- subset(seurat_sce, subset = is.low_yield == FALSE & is.damaged == FALSE)
sce <- sce[, sce$is.low_yield == FALSE & sce$is.damaged == FALSE]
print(paste0("Keeping ",ncol(sce)," Cells from SCE object and ", ncol(seurat_sce), " Cells from Seurat object."))


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
seurat_sce <- sctransform_Seurat(seurat_sce)
### Run initial dim reduction for norm
seurat_sce <- reduceDims_Seurat(seurat_sce)
### Run initial dim reduction for STC
seurat_sce <- reduceDims_Seurat(seurat_sce, assay = "SCT", reduction.name = "pca_sct")
### Cluster 
seurat_sce <- cluster_Seurat(seurat_sce, assay = "RNA", reduction = "pca", cluster.name = "pca_cluster")
seurat_sce <- cluster_Seurat(seurat_sce, assay = "SCT", reduction = "pca_sct", cluster.name = "sct_cluster")
## Run UMAPs
seurat_sce <- umap_Seurat(seurat_sce, assay = "RNA", reduction = "pca", dims = 1:30, reduction.name = "umap_pca", n.neighbors = 30L, min.dist = 0.1, spread = 5)
seurat_sce <- umap_Seurat(seurat_sce, assay = "SCT", reduction = "pca_sct", dims = 1:30, reduction.name = "umap_pca_sct", n.neighbors = 30L, min.dist = 0.1, spread = 5)

#### Assign CaTCH barcode indices based on their abundance in the reference samples ####
tmp <- colData(sce) %>%
  as_tibble() %>%
  filter(CaTCH.Status == "Singlet", Condition == opt$baseCond) %>%
  select(CaTCH.BCs, Sample) %>%
  group_by(CaTCH.BCs, Sample) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  distinct() %>%
  pivot_wider(names_from = Sample, values_from = n, values_fill = 0) %>%
  filter(rowSums(across(starts_with(opt$baseCond))) > 0) %>%
  mutate(.Means = rowMeans(across(starts_with(opt$baseCond)))) %>%
  arrange(by = desc(.Means)) %>%
  rowid_to_column(".ID") %>%
  mutate(CaTCH.BC_ID = paste0("BC_", .ID)) %>%
  select(-.Means, -.ID) %>%
  relocate(CaTCH.BC_ID, .after = CaTCH.BCs) %>%
  select(CaTCH.BCs, CaTCH.BC_ID)

colData(sce)["CaTCH.BC_ID"] <- colData(sce) %>%
  as_tibble() %>%
  left_join(y = tmp, by = "CaTCH.BCs") %>%
  mutate(CaTCH.BC_ID = factor(CaTCH.BC_ID, levels = str_sort(unique(CaTCH.BC_ID), numeric = TRUE))) %>%
  select(CaTCH.BC_ID)

seurat_sce@meta.data$CaTCH.BC_ID <- seurat_sce@meta.data %>%
  left_join(y = tmp, by = "CaTCH.BCs") %>%
  mutate(CaTCH.BC_ID = factor(CaTCH.BC_ID, levels = str_sort(unique(CaTCH.BC_ID), numeric = TRUE))) %>%
  pull(CaTCH.BC_ID)

#### Save final objects ####
saveRDS(sce, file = paste0(opt$out, "_filtered_sce.rds"))
saveRDS(seurat_sce, file = paste0(opt$out, "_filtered_seurat_sce.rds"))

## Prepare the data for the report plots
#tmp <- reducedDim(sce, "TSNE", withDimnames = FALSE)
#data.tsne <- tibble(
#  tSNE1 = tmp[, 1],
#  tSNE2 = tmp[, 2],
#  Sample = sce$Sample
#) %>%
#  write.table(x = ., file = paste0(opt$out, ".tsne"), quote = FALSE, row.names = FALSE)
#
#colData(sce) %>%
#  write.table(x = ., file = paste0(opt$out, ".metadata"), quote = FALSE, row.names = FALSE)