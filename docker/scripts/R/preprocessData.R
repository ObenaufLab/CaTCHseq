.libPaths("")

library(optparse)

option_list = list(
  make_option(c("--sample"), type="character", default = NULL, 
              help = "sample name"),
  make_option(c("--data10X"), type="character", default = NULL, 
              help = "path to the 10X data"),
  make_option(c("--catchBC"), type="character", default = NULL, 
              help = "path to the file with CaTCH barcodes"),
  make_option(c("--reference"), type="character", default = NULL, 
              help = "name of the treatment to use as the reference"),
  make_option(c("--gtf"), type="character", default = NULL, 
              help = "path to the annotation GTF file"),        
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



if (is.null(opt$samples) || is.null(opt$out)){
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

load("/tools/data/R/stagemarkers_xue2020.rda")

print(paste("Processing the sample ", opt$sample, sep = ""))
print("   Loading 10X data ...")
sceData <- Seurat::Read10X(data.dir = opt$data10X,
                           gene.column = 2,
                           cell.column = 1,
                           strip.suffix = TRUE,
                           unique.features = TRUE)
sce <- SingleCellExperiment(assays = list(counts = sceData))

# Add both gene symbol and gene ID to rowData
rowData(sce)[c("GeneID", "Symbol")] <- read_tsv(file = paste(data10X, "features.tsv.gz", sep = "/"),
                                                col_names = c("GeneID", "Symbol"),
                                                col_types = "cc-")

# Add the biotype and the locus data
gtf <- read_tsv(opt$gtf,
                comment = "#",
                col_names = c("chr", "type", "start", "end", "strand", "attr"),
                col_types = "c-cdd-c-c") %>%
       filter(type == "gene") %>%
       mutate(gene_id = str_extract(attr, 'gene_id "[^"]+"'),
              biotype = str_extract(attr, 'gene_biotype "[^"]+"')) %>%
       select(-attr, -type) %>%
       mutate(gene_id = str_extract(gene_id, '"[^"]+"'),
              gene_id = str_remove_all(gene_id, "\""),
              biotype = str_extract(biotype, '"[^"]+"'),
              biotype = str_remove_all(biotype, "\""),
              locus = sprintf("%s%s:%d-%d", strand, chr, start, end))

rowData(sce)[c("Biotype", "Locus")] <- rowData(sce) %>%
                                       as_tibble() %>%
                                       left_join(y = gtf, by = c("GeneID" = "gene_id")) %>%
                                       select(biotype, locus)


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




#### Annotate the treatments and replicates ####
annot = read_tsv(file = opt$csv,
              col_names = c("Sample", "Treatment", "Replicate", "Type"),
              col_types = "cccc---",
              skip = 1) %>%
        filter(Type == "GEX") %>%
        select(-Type) %>%
        distinct()
# Sanity check: the number of samples in the variable `annot` must be the same as
# the number of distinct samples in the SingleCellExperiment object
if (nrow(annot) != length(unique(sce$Sample))) {
  warn("The number of samples in the CSV file does not match the number of distinct samples in SingleCellExperiment object")
}
colData(sce)[, c("Treatment", "Replicate")] <- colData(sce) %>%
                                               as_tibble() %>%
                                               left_join(y = annot, by = "Sample") %>%
                                               select(Treatment, Replicate) %>%
                                               mutate(Treatment = factor(Treatment, levels = c(opt$reference, setdiff(unique(Treatment), opt$reference))))

tmp <- colData(sce) %>% 
       as_tibble() %>% 
       select(Sample, Treatment) %>% 
       distinct() %>% 
       arrange(Treatment) %>% 
       select(Sample) %>% 
       pull()
sce$Sample <- factor(sce$Sample, levels = tmp)



sce.unfiltered <- sce
mask <- sce$Category == "Good"
sce <- sce[, mask]



#### Run PCA, tSNE and UMAP ####
print("Identify the top variable genes...")
gene.var <- modelGeneVar(sce)
hvg <- getTopHVGs(stats = gene.var, prop = opt$hvg_cutoff)

print("Run the PCA ...")
sce <- runPCA(sce, subset_row = hvg)
set.seed(101)

print("Run tSNE and UMAP analyses ...")
sce <- runTSNE(sce, dimred = "PCA")
sce <- runUMAP(sce, dimred = "PCA")

print("Clustering ...")
g <- buildSNNGraph(sce, use.dimred = "PCA")
cluster <- igraph::cluster_walktrap(g)$membership
colData(sce)["Cluster"] <- factor(cluster)


#### Assign CaTCH barcode indices based on their abundance in the reference samples ####
tmp <- colData(sce) %>% 
       as_tibble() %>%
       filter(CaTCH.Status == "Singlet", Treatment == opt$reference) %>%
       select(CaTCH.BCs, Sample) %>%
       group_by(CaTCH.BCs, Sample) %>%
       mutate(n = n()) %>%
       ungroup() %>%
       distinct() %>%
       pivot_wider(names_from = Sample, values_from = n, values_fill = 0) %>%
       filter(rowSums(across(starts_with(opt$reference))) > 0) %>% 
       mutate(.Means = rowMeans(across(starts_with(opt$reference)))) %>%
       arrange(by = desc(.Means)) %>%
       rowid_to_column(".ID") %>%
       mutate(BC_ID = paste0("BC_", .ID)) %>%
       select(-.Means, -.ID) %>%
       relocate(BC_ID, .after = CaTCH.BCs) %>%
       select(CaTCH.BCs, BC_ID)

colData(sce)["BC_ID"] <-  colData(sce) %>%
                          as_tibble() %>%
                          left_join(y = tmp, by = "CaTCH.BCs") %>%
                          mutate(BC_ID = factor(BC_ID, levels = str_sort(unique(BC_ID), numeric = TRUE))) %>%
                          select(BC_ID)

save(sce.unfiltered, sce, file = opt$out)

# Prepare the data for the report plots
tmp <- reducedDim(sce, "TSNE", withDimnames = FALSE)
data.tsne <- tibble(tSNE1 = tmp[, 1], 
                    tSNE2 = tmp[, 2], 
                    Sample = sce$Sample) %>%
             write.table(x = ., file = paste0(opt$out, ".tsne"), quote = FALSE, row.names = FALSE)

colData(sce) %>%
  write.table(x = ., file = paste0(opt$out, ".metadata"), quote = FALSE, row.names = FALSE)

