parse_object_flag <- function(value, name) {
    value <- tolower(as.character(value))
    if (length(value) != 1L || is.na(value) || !value %in% c("true", "false")) {
        stop(sprintf("--%s must be true or false", name))
    }
    identical(value, "true")
}

anndata_transfer_object <- function(object) {
    rna <- SeuratObject::JoinLayers(object[["RNA"]])
    counts <- SeuratObject::LayerData(rna, layer = "counts")
    data <- SeuratObject::LayerData(rna, layer = "data")
    cells <- colnames(counts)
    genes <- rownames(counts)
    if (!identical(dimnames(counts), dimnames(data))) {
        stop("RNA counts and normalized data must have matching genes and cells")
    }
    result <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts, logcounts = data),
        colData = S4Vectors::DataFrame(object[[]][cells, , drop = FALSE]),
        rowData = S4Vectors::DataFrame(rna[[]][genes, , drop = FALSE])
    )
    for (name in SeuratObject::Reductions(object)) {
        embedding <- SeuratObject::Embeddings(object[[name]])
        SingleCellExperiment::reducedDim(result, paste0("X_", name)) <- embedding[cells, , drop = FALSE]
    }
    result
}

write_anndata <- function(object, path) {
    if (!requireNamespace("anndataR", quietly = TRUE)) {
        stop("AnnData output requires the Bioconductor package 'anndataR'; rebuild the CaTCHseq container")
    }
    anndataR::write_h5ad(anndata_transfer_object(object), path,
        x_mapping = "logcounts", compression = "gzip"
    )
}
