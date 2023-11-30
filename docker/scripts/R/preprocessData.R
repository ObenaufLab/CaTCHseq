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


#### Additional methods ####
activity_by_logmeans <- function(signature, sce, fn.scale = scale, ...) {
    # check if logcounts is present
    stopifnot("logcounts" %in% names(assays(sce)))

    # determine shared features between sce and signature
    shared_features <- intersect(signature, rownames(sce))
    stopifnot(length(shared_features) > 0)

    # calculate scaled log-means of signature
    return(sce[shared_features, ] %>%
        assay("logcounts") %>%
        Matrix::colMeans() %>%
        # use scale by default, but allow other functions and additional arguments
        fn.scale(...) %>%
        as.numeric())
}


assignCategoryByMarker <- function(sce, markers = NULL, col.name = NULL, fn.scale = scale, ..., center = TRUE) {
    if (is.null(sce)) {
        stop("'sce' must be specified and cannot be NULL")
    }
    if (is.null(markers)) {
        stop("'markers' must be specified and cannot be NULL. You can use the provided dataset 'stagemarkers_xue2020'")
    }
    if (is.null(col.name)) {
        stop("'col.name' must be specified and cannot be NULL")
    }

    categories <- map_df(markers, activity_by_logmeans, sce = sce, fn.scale = fn.scale, ...)

    # categories should be a matrix with N rows and M columns, whereas N is the
    # number of cells in the SingleCellExperiment object `sce` and M is the
    # number of names of the markers list.
    # If this is not the case after scaling, transpose it.
    dims <- dim(categories)
    if (dims[1] != dim(sce)[2]) {
        categories <- categories %>%
            t()
    }

    # If the results should be centered around 0, use `scale` to do that on a
    # per-cell basis and transpose the result.
    if (center) {
        categories <- categories %>%
            apply(1, scale) %>%
            t()
    }
    colnames(categories) <- names(markers)
    colData(sce)[col.name] <- colnames(categories)[apply(categories, 1, which.max)]
    return(sce)
}


### Create sce objects with corresponding CaTCH barcodes ###
create_SCEs <- function(smpl, data10X, bc) {
    samplelist <- str_split(smpl, ",")[[1]]
    data10Xs <- str_split(data10X, ",")[[1]]
    bcs <- str_split(bc, ",")[[1]]

    if (length(samplelist) == 1) {
        print(paste("Processing the sample ", smpl, sep = ""))
        print("   Loading 10X data ...")
        sceData <- Seurat::Read10X(
            data.dir = data10X,
            gene.column = 2,
            cell.column = 1,
            strip.suffix = TRUE
        )
        sce <- SingleCellExperiment(assays = list(counts = sceData))
        sce$Sample <- smpl
        sce$CellID <- colnames(sce)
        print(paste("   Loaded expression data for ", ncol(sce), " cells", sep = ""))

        #### Now load the CaTCH barcodes ####
        data.catch <- load_BCs(bc)
        print(paste0(" Loading Barcodes from ", bc))
        # print(head(data.catch))
        # print(head(colData(sce) %>%
        #    as_tibble() %>%
        #    left_join(y = data.catch, by = "CellID")))

        colData(sce)[, c("CaTCH.BCs", "CaTCH.BC.counts")] <- colData(sce) %>%
            as_tibble() %>%
            left_join(y = data.catch, by = "CellID") %>%
            select(CaTCH.BCs, CaTCH.BC.counts)


        #### Count the CaTCH barcode reads and find the abundance of the two most abundant CaTCH barcodes ####
        bc.data <- lapply(X = sce$CaTCH.BC.counts, FUN = function(x) {
            v <- str_split(string = x, pattern = ";", simplify = TRUE) %>%
                as.numeric()
            return(c(sum(v), v[1], v[2]))
        }) %>%
            unlist() %>%
            matrix(data = ., ncol = 3, byrow = TRUE)
        colnames(bc.data) <- c("CaTCH.Sum", "CaTCH.BC1", "CaTCH.BC2")
        colData(sce) <- cbind(colData(sce), bc.data)
        colData(sce)["CaTCH.Status"] <- colData(sce) %>%
            as_tibble() %>%
            mutate(
                CaTCH.Status = if_else(is.na(CaTCH.Sum),
                    "No barcode",
                    if_else(is.na(CaTCH.BC2),
                        "Singlet",
                        if_else(CaTCH.BC1 / CaTCH.Sum >= 0.8,
                            "Putative singlet",
                            "Multiplet"
                        )
                    )
                ),
                CaTCH.Status = factor(CaTCH.Status, levels = c(
                    "Singlet",
                    "Putative singlet",
                    "Multiplet",
                    "No barcode"
                ))
            ) %>%
            select(CaTCH.Status)

        sce <- as.Seurat(sce, data = NULL, assay = NULL)
        sce <- RenameAssays(sce, assay.name = "originalexp", new.assay.name = "RNA")
        sce[["RNA"]] <- as(object = sce[["RNA"]], Class = "Assay5")

        return(sce)
    } else {
        scetomerge <- list()
        for (i in 1:length(samplelist)) {
            if (i == 1) {
                print("Loading first sce dataset")
                firstsce <- create_SCEs(samplelist[i], data10Xs[i], bcs[i])
            } else {
                print(paste0("Loading sce dataset ", i))
                scetomerge[[i - 1]] <- create_SCEs(samplelist[i], data10Xs[i], bcs[i])
            }
        }

        sce <- merge(firstsce, scetomerge, add.cell.ids = samplelist, project = "scCaTCH")

        return(sce)
    }
}

load_BCs <- function(bc) {
    print("   Loading CaTCH barcodes ...")
    data.catch <- read_tsv(
        file = bc,
        col_names = c("CellID", "CaTCH.BCs", "CaTCH.BC.counts"),
        col_types = "ccc",
        skip = 1
    )
    return(data.catch)
}

normalize_SCE <- function(sce) {
    print("   Normalizing and scaling SCE ...")

    sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
    sce <- FindVariableFeatures(sce)
    sce <- ScaleData(sce)

    return(sce)
}

sctransform_SCE <- function(sce) {
    sce <- SCTransform(sce, vst.flavor = "v2", vars.to.regress = "percent.mt", verbose = FALSE)
    # sce[["SCT"]] <- as(object = sce[["SCT"]], Class = "Assay5")  # this actually breaks downstream code as SCTAssay is still v3 (https://github.com/satijalab/seurat/issues/7542)
    return(sce)
}

reduceDims_SCE <- function(sce, assay = "RNA", reduction.name = "pca") {
    print("   Reducing dimensions ...")

    sce <- RunPCA(sce, assay = assay, reduction.name = reduction.name)
    return(sce)
}

cluster_SCE <- function(sce, assay = "RNA_integrated.cca", reduction = "integrated.cca", cluster.name = "integrated.cca_cluster", resolution = .2) {
    print("   Clustering ...")

    sce <- FindNeighbors(sce, assay = assay, reduction = reduction, compute.SNN = TRUE, graph.name = c(paste0(assay, "_nn"), paste0(assay, "_snn")))
    sce <- FindClusters(sce, resolution = .2, cluster.name = cluster.name, graph.name = paste0(assay, "_snn"))

    return(sce)
}

integrate_SCE <- function(sce, assay = "RNA", method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", normalization.method = "LogNormalize", ...) {
    print("   Integrating samples ...")

    sce <- IntegrateLayers(object = sce, assay = assay, method = method, orig.reduction = orig.reduction, new.reduction = new.reduction, normalization.method = normalization.method, verbose = FALSE, ...)
    # re-join layers after integration
    sce[[paste0(assay, "_", new.reduction)]] <- JoinLayers(sce[[assay]])

    return(sce)
}

umap_SCE <- function(sce, assay = "RNA", reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca", n.neighbors = 25L, min.dist = 0.1, spread = 5) {
    print("   Preparing UMAP ...")

    sce <- RunUMAP(sce, assay = assay, reduction = reduction, dims = dims, reduction.name = reduction.name, n.neighbors = n.neighbors, min.dist = min.dist, spread = spread)
    return(sce)
}

# Split for integrated analysis
split_SCE <- function(sce, by = "Sample") {
    print("   Splitting SCE by Sample ...")

    sce[["RNA"]] <- split(sce[["RNA_split"]], f = sce[by])
    return(sce)
}

# Join for DE
join_SCE <- function(sce, assay = "RNA", layers = "data", new = "joined_data") {
    print("   Joining SCE layers  ...")

    sce <- JoinLayers(sce, assay = assay, layers = layers, new = new)
    return(sce)
}

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
sce <- PercentageFeatureSet(sce, pattern = "^MT-", col.name = "percent.mt")
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

pdf(file = paste0(opt$out, "_QC_Scatter_MT_content.pdf"), width = 30, height = 10)
FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Sample")
dev.off()

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
# sce <- sctransform_SCE(sce)

### Run initial dim reduction for norm
sce <- reduceDims_SCE(sce)

### Run initial dim reduction for STC
# sce <- reduceDims_SCE(sce, assay = "SCT", reduction.name = "pca_sct")

### Integrate the SCE layers for integrative analysis
sce <- integrate_SCE(sce, assay = "RNA", orig.reduction = "pca", new.reduction = "integrated.cca", normalization.method = "LogNormalize")

### Integrate the SCE layers after SCT for integrative analysis
# sce <- integrate_SCE(sce, assay = "SCT", orig.reduction = "pca_sct", new.reduction = "sct_integrated.cca", normalization.method = "SCT")

### Cluster integrated SCE for plotting
# sce <- cluster_SCE(sce, assay = "SCT_integrated.cca", reduction = "pca_sct", cluster.name = "sct_integrated.cca_cluster")

sce <- cluster_SCE(sce, assay = "RNA_integrated.cca", reduction = "integrated.cca", cluster.name = "integrated.cca_cluster")
sce <- cluster_SCE(sce, assay = "RNA", reduction = "pca", cluster.name = "pca_cluster")
## Run UMAP not integrated
sce <- umap_SCE(sce, assay = "RNA_integrated.cca", reduction = "integrated.cca", dims = 1:30, reduction.name = "umap_integrated.cca", n.neighbors = 30L, min.dist = 0.1, spread = 5)
sce <- umap_SCE(sce, assay = "RNA", reduction = "pca", dims = 1:30, reduction.name = "umap_pca", n.neighbors = 30L, min.dist = 0.1, spread = 5)
## Run UMAP integrated
# sce <- umap_SCE(sce, assay = "SCT_integrated.cca", reduction = "sct_integrated.cca", dims = 1:10, reduction.name = "umap_integrated.cca", n.neighbors = 30L, min.dist = 0.01, spread = 5)

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