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
              help = "cutoff value to call high variable genes (must be between 0 and 1)"),
  make_option(c("--out"), type="character", default = NULL, 
              help = "path to the output file")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)


#opt <- list()
#opt$sample <- "/data/gcbds/users/nowoshil/Projects/CaTCH2.0/KRASi_ON_vs_OFF/Exp1220/CaTCH/12d/manual/samples"
#opt$csv <- "/data/gcbds/users/nowoshil/Projects/CaTCH2.0/KRASi_ON_vs_OFF/Exp1220/CaTCH/12d/libraries.tsv"
#opt$reference <- "Day0"
#opt$max_mt <- 10
#opt$min_features <- 500
#opt$hvg_cutoff <- 0.1
#opt$out <- "~/tmp/sce_CaTCH_Exp1220.rda"



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
############################



library(tidyverse)
library(scater)
library(scran)
library(SingleCellExperiment)
library(Seurat)

load("/tools/data/R/stagemarkers_xue2020.rda")
#load("/home/nowoshil/Repositories/nf-pipelines/pipelines-singlecell-catch-nf/docker/data/R/stagemarkers_xue2020.rda")

print(paste("Processing the sample ", opt$sample, sep = ""))
print("   Loading 10X data ...")
sceData <- Seurat::Read10X(data.dir = opt$data10X,
                           gene.column = 2,
                           cell.column = 1,
                           strip.suffix = TRUE)
sce <- SingleCellExperiment(assays = list(counts = sceData))
sce$Sample <- opt$sample
sce$CellID <- colnames(sce)
print(paste("   Loaded expression data for ", ncol(sce), " cells", sep = ""))

print("   Loading CaTCH barcodes ...")
data.catch <- read_tsv(file = opt$catchBC, 
                       col_names = c("CellID", "CaTCH.BCs", "CaTCH.BC.counts"),
                       col_types = "ccc",
                       skip = 1)

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

#### Normalize the data and assign cell categories ####
is.mitochondrial <- grepl(pattern = "^MT-",
                          x = rownames(sce),
                          ignore.case = FALSE)
cell.categories <- c("Good", "Damaged", "Few features", "Damaged AND few features")

print("Normalize the counts ...")
sce <- sce %>%
       scater::logNormCounts() %>%
       scater::addPerCellQC(subsets = list(MT = is.mitochondrial)) %>%
       scater::addPerFeatureQC() %>%
       assignCategoryByMarker(markers = stagemarkers_xue2020, col.name = "CellStage")

print("Categorize the cells ...")
rowData(sce)["is.mitochondrial"] <- is.mitochondrial
colData(sce)["is.damaged"] <- colData(sce)[, "subsets_MT_percent"] > opt$max_mt
colData(sce)["is.low_yield"] <- colData(sce)[, "detected"] < opt$min_features
colData(sce)["Category"] <- colData(sce) %>%
                            tibble::as_tibble() %>%
                            dplyr::mutate(tmp = as.integer(is.damaged) * 1 + as.integer(is.low_yield) * 2) %>%
                            dplyr::select("tmp") %>%
                            purrr::map(.f = ~ cell.categories[.x + 1]) %>%
                            tibble::as_tibble() %>%
                            dplyr::mutate(Class = factor(tmp, levels = cell.categories)) %>%
                            dplyr::select(Class)

save(sce, sce, file = opt$out)