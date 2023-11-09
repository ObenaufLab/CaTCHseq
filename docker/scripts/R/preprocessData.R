.libPaths("")

library(optparse)

option_list = list(
  make_option(c("--sample"), type="character", default = NULL, 
              help = "sample name"),
  make_option(c("--data10X"), type="character", default = NULL, 
              help = "path to the 10X data"),
  make_option(c("--catchBC"), type="character", default = NULL, 
              help = "path to the file with CaTCH barcodes"),
  make_option(c("--max_mt"), type = "numeric", default = 10, 
              help = "maximum percent of mitochondrial reads in a valid cell"),
  make_option(c("--min_features"), type = "numeric", default = 500, 
              help = "minimum number of detected features in a valid cell"),
  make_option(c("--hvg_cutoff"), type = "numeric", default = 0.1, 
              help = "cutoff value to call percentage of high variable genes (must be between 0 and 1)"),
  make_option(c("--out"), type="character", default = NULL, 
              help = "path to the output file"),
  make_option(c("--marker"), type="character", default = '/tools/data/R/stagemarkers_xue2020.rda', 
              help = "path to the stagemarker file")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if (is.null(opt$sample) || is.null(opt$out)){
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
  return (sce[shared_features, ] %>%
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
  return (sce)
}


### Create sce objects with corresponding CaTCH barcodes ###
create_SCEs <- function(smpl, data10X, bc){
    samplelist <- str_split(smpl, ',')[[1]]
    data10Xs <- str_split(data10X, ',')[[1]]
    bcs <- str_split(bc, ',')[[1]]
    
    if (length(samplelist) == 1){
        print(paste("Processing the sample ", smpl, sep = ""))
        print("   Loading 10X data ...")
        sceData <- Seurat::Read10X(data.dir = data10X,
                                   gene.column = 2,
                                   cell.column = 1,
                                   strip.suffix = TRUE)
        sce <- SingleCellExperiment(assays = list(counts = sceData))
        sce$Sample <- smpl
        sce$CellID <- colnames(sce)
        print(paste("   Loaded expression data for ", ncol(sce), " cells", sep = ""))
        
        #### Now load the CaTCH barcodes ####
        data.catch <- load_BCs(bc)
        
        colData(sce)[,c("CaTCH.BCs", "CaTCH.BC.counts")] <- colData(sce) %>%
            as_tibble() %>%
            left_join(y = data.catch, by = "CellID") %>%
            select(CaTCH.BCs, CaTCH.BC.counts)
        
        
        #### Count the CaTCH barcode reads and find the abundance of the two most abundant CaTCH barcodes ####
        bc.data <- lapply(X = sce$CaTCH.BC.counts, FUN = function(x){
            v <- str_split(string = x, pattern = ";", simplify = TRUE) %>% 
                as.numeric()
            return (c(sum(v), v[1], v[2]))
        }) %>% 
            unlist() %>% 
            matrix(data = ., ncol = 3, byrow = TRUE)
        colnames(bc.data) <- c("CaTCH.Sum", "CaTCH.BC1", "CaTCH.BC2")
        colData(sce) <- cbind(colData(sce), bc.data)
        colData(sce)["CaTCH.Status"] <- colData(sce) %>%
            as_tibble() %>%
            mutate(CaTCH.Status = if_else(is.na(CaTCH.Sum), 
                                          "No barcode", 
                                          if_else(is.na(CaTCH.BC2), 
                                                  "Singlet",
                                                  if_else(CaTCH.BC1 / CaTCH.Sum >= 0.8, 
                                                          "Putative singlet", 
                                                          "Multiplet"))),
                   CaTCH.Status = factor(CaTCH.Status, levels = c("Singlet", 
                                                                  "Putative singlet", 
                                                                  "Multiplet",
                                                                  "No barcode"))) %>%
            select(CaTCH.Status)
        
        sce <- as.Seurat(sce, data=NULL, assay=NULL)
        sce <- RenameAssays(sce, assay.name = "originalexp", new.assay.name = "RNA")
        sce[["RNA"]] <- as(object = sce[["RNA"]], Class = "Assay5")
        
        return (sce)
        
    }else{
        scetomerge <- list()
        for (i in 1:length(samplelist)){
            if (i == 1){
                print("Loading first sce dataset")
                firstsce <- create_SCEs(samplelist[i], data10Xs[i], bcs[i])
            }else{
                print(paste0("Loading sce dataset ",i))
                scetomerge[[i-1]] <- create_SCEs(samplelist[i], data10Xs[i], bcs[i])
            }
        }
        
        sce <- merge(firstsce, scetomerge, add.cell.ids=samplelist, project = "scCaTCH")
        
        return(sce)
    }
}

load_BCs <- function(bc){
    print("   Loading CaTCH barcodes ...")
    data.catch <- read_tsv(file = bc, 
                           col_names = c("CellID", "CaTCH.BCs", "CaTCH.BC.counts"),
                           col_types = "ccc",
                           skip = 1)
    return (data.catch)
    
}

reduceDims_SCE <- function(sce){
    print("   Reducing dimensions ...")
        
    sce <- RunPCA(sce)
    return(sce)
}

cluster_SCE <- function(sce){
    print("   Clustering ...")
    
    sce <- FindNeighbors(sce, reduction = "integrated.cca", dims = 1:30)
    sce <- FindClusters(sce, resolution = 1)
    
    return(sce)
    
}

normalize_SCE <- function(sce){
    print("   Normalizing and scaling SCE ...")
        
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce)
    sce <- ScaleData(sce)
    
    return(sce)
    
}

integrate_SCE <- function(sce){
    print("   Integrating samples ...")
    
    sce <- Seurat::IntegrateLayers(object = sce, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                    verbose = FALSE)
    # re-join layers after integration
    sce[["RNA"]] <- JoinLayers(sce[["RNA"]])
    
    return(sce)
}

umap_SCE <- function(sce){
    print("   Preparing UMAP ...")
    
    sce <- RunUMAP(sce, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
    return(sce)
}

# Split for integrated analysis
split_SCE <- function(sce){
    print("   Splitting SCE by Sample ...")
    
    sce[["RNA"]] <- split(sce[["RNA"]], f=sce$Sample)
    return(sce)
}

# Join for DE
join_SCE <- function(sce){
    print("   Joining SCE layers  ...")
    
    sce <- JoinLayers(sce)
    return(sce)
}

############################

library(tidyverse)
library(scater)
library(scran)
library(SingleCellExperiment)
library(Seurat)
library(SeuratWrappers)
#options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")

load(opt$marker)

### Load all experiments, add CaTCH barcodes as layer and converting to Seuratv5 Object
sce <- create_SCEs(opt$sample, opt$data10X, opt$catchBC)

### Split SCE for integrative analysis
### ALREADY DONE BY MERGE
# sce <- split_SCE(sce)

### Normalize and scale counts
sce <- normalize_SCE(sce)

### Run initial dim reduction
sce <- reduceDims_SCE(sce)

### Integrate the SCE layers for integrative analysis
sce <- integrate_SCE(sce)

### Cluster SCE for plotting
sce <- cluster_SCE(sce)



#### Normalize the data and assign cell categories ####
is.mitochondrial <- grepl(pattern = "^MT-",
                          x = rownames(sce),
                          ignore.case = FALSE)
cell.categories <- c("Good", "Damaged", "Few features", "Damaged AND few features")

#### Split into subobjects ####
sce.list <- SplitObject(sce, split.by = "Sample")

print("Normalize the counts ...")
sce.list <- lapply(X = sce.list, FUN = function(x) {
    x <- as.SingleCellExperiment(x) %>% 
        scater::logNormCounts() %>%
        scater::addPerCellQC(subsets = list(MT = is.mitochondrial)) %>%
        scater::addPerFeatureQC() %>%
        assignCategoryByMarker(markers = stagemarkers_xue2020, col.name = "CellStage")
})

print("Categorize the cells ...")
for (i in 1:length(sce.list)){
    rowData(sce.list[[i]])$is.mitochondrial <- is.mitochondrial
    colData(sce.list[[i]])$is.damaged <- colData(sce.list[[i]])[, "subsets_MT_percent"] > opt$max_mt
    colData(sce.list[[i]])$is.low_yield <- colData(sce.list[[i]])[, "detected"] < opt$min_features
    colData(sce.list[[i]])$Category <- colData(sce.list[[i]]) %>%
        tibble::as_tibble() %>%
        dplyr::mutate(tmp = as.integer(is.damaged) * 1 + as.integer(is.low_yield) * 2) %>%
        dplyr::select("tmp") %>%
        purrr::map(.f = ~ cell.categories[.x + 1]) %>%
        tibble::as_tibble() %>%
        dplyr::mutate(Class = factor(tmp, levels = cell.categories)) %>%
        dplyr::select(Class)
}
        
#### Run PCA, tSNE and UMAP ####
print("Identify the top variable genes...")
for (i in 1:length(sce.list)){
    gene.var <- scran::modelGeneVar(sce.list[[i]])
    hvg <- scran::getTopHVGs(stats = gene.var, prop = opt$hvg_cutoff)
    sce.list[[i]] <- scater::runPCA(sce.list[[i]], subset_row = hvg)
}

set.seed(101)

print("Run tSNE and UMAP analyses ...")
for (i in 1:length(sce.list)){
    sce.list[[i]] <- scater::runTSNE(sce.list[[i]], dimred = "PCA")
    sce.list[[i]] <- scater::runUMAP(sce.list[[i]], dimred = "PCA")
}

print("Clustering ...")
for (i in 1:length(sce.list)){
    g <- scran::buildSNNGraph(sce.list[[i]], use.dimred = "PCA")
    cluster <- igraph::cluster_walktrap(g)$membership
    
    colData(sce.list[[i]])$Cluster <- factor(cluster)
}

# Prepare the data for the report plots
for (i in 1:length(sce.list)){
    tmp <- reducedDim(sce.list[[i]], "TSNE", withDimnames = FALSE)
    tibble(tSNE1 = tmp[, 1], 
                tSNE2 = tmp[, 2], 
                Sample = sce.list[[i]]$Sample) %>%
                write.table(x = ., file = gzfile(paste0(opt$out, sce.list[[i]]$Sample[1],".tsne.gz")), quote = FALSE, row.names = FALSE)

    colData(sce.list[[i]]) %>%
        write.table(x = ., file = gzfile(paste0(opt$out, sce.list[[i]]$Sample[1],".metadata.gz")), quote = FALSE, row.names = FALSE)
}

for (i in 1:length(sce.list)){
    sce.list[[i]] <- as.Seurat(sce.list[[i]])    
}

sce <- merge(sce.list[[1]], y=sce.list[-1], add.cell.ids=names(sce.list), merge.data=TRUE, merge.dr=TRUE)

### Reduce Dims for overall sce ###
gene.var <- scran::modelGeneVar(as.SingleCellExperiment(sce))
hvg <- scran::getTopHVGs(stats = gene.var, prop = opt$hvg_cutoff)
sce <- scater::runPCA(as.SingleCellExperiment(sce), subset_row = hvg)

set.seed(101)

print("Run tSNE and UMAP analyses ...")
sce <- scater::runTSNE(sce, dimred = "PCA")
sce <- scater::runUMAP(sce, dimred = "PCA")

print("Clustering ...")
g <- scran::buildSNNGraph(sce, use.dimred = "PCA")
cluster <- igraph::cluster_walktrap(g)$membership
colData(sce)$Cluster <- factor(cluster)

### Convert to Seurat Object for later
sce <- as.Seurat(sce)

save(sce, sce.list, file = paste0(opt$out,'SCE.rda.gz'), compress = "gzip")
