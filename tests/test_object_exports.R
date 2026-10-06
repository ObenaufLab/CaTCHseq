source("docker/scripts/R/object_exports.R")

stopifnot(identical(parse_object_flag("true", "x"), TRUE),
    identical(parse_object_flag("TRUE", "x"), TRUE),
    identical(parse_object_flag("True", "x"), TRUE),
    identical(parse_object_flag(TRUE, "x"), TRUE),
    identical(parse_object_flag("false", "x"), FALSE),
    identical(parse_object_flag("FALSE", "x"), FALSE),
    identical(parse_object_flag("False", "x"), FALSE),
    identical(parse_object_flag(FALSE, "x"), FALSE))
for (bad in list(NULL, NA, c("true", "false"), "junk", 1)) {
    err <- tryCatch(parse_object_flag(bad, "x"), error = function(e) conditionMessage(e))
    stopifnot(is.character(err), grepl("must be true or false", err))
}
cat("parse_object_flag OK\n")

suppressPackageStartupMessages({
    library(Seurat)
    library(SingleCellExperiment)
})
options(Seurat.object.assay.version = "v5")

set.seed(1)
ngenes <- 10L
ncells <- 8L
counts <- as(matrix(rpois(ngenes * ncells, 2), nrow = ngenes, ncol = ncells,
    dimnames = list(paste0("G", seq_len(ngenes)), paste0("C", seq_len(ncells)))), "CsparseMatrix")
obj <- CreateSeuratObject(counts = counts, assay = "RNA")
obj <- NormalizeData(obj)
data_orig <- GetAssayData(obj, layer = "data")
obj[["RNA"]][["feat_meta"]] <- setNames(seq_len(ngenes), rownames(obj[["RNA"]]))
obj$cell_meta <- paste0("m", seq_len(ncells))
obj$sample <- rep(c("s1", "s2"), each = 4L)
pca <- matrix(rnorm(ncells * 5L), nrow = ncells, ncol = 5L,
    dimnames = list(colnames(obj), paste0("PC_", 1:5)))
obj[["pca"]] <- CreateDimReducObject(embeddings = pca, key = "PC_", assay = "RNA")
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sample)
stopifnot(identical(sort(Layers(obj[["RNA"]])), sort(c("counts.s1", "counts.s2", "data.s1", "data.s2"))))

sce <- anndata_transfer_object(obj)

stopifnot(is(sce, "SingleCellExperiment"),
    is(assay(sce, "counts"), "sparseMatrix"), is(assay(sce, "logcounts"), "sparseMatrix"))
stopifnot(identical(dimnames(assay(sce, "counts")), dimnames(counts)))
stopifnot(identical(as.matrix(assay(sce, "counts")), as.matrix(counts)))
stopifnot(identical(as.matrix(assay(sce, "logcounts")), as.matrix(data_orig)))
stopifnot(identical(unname(sce$cell_meta), unname(obj$cell_meta)), identical(unname(sce$sample), unname(obj$sample)),
    all(rowData(sce)$feat_meta == seq_len(ngenes)))
stopifnot(identical(reducedDimNames(sce), "X_pca"),
    isTRUE(all.equal(reducedDim(sce, "X_pca"), Embeddings(obj[["pca"]]))))
stopifnot(identical(sort(Layers(obj[["RNA"]])), sort(c("counts.s1", "counts.s2", "data.s1", "data.s2"))),
    identical(as.matrix(LayerData(obj[["RNA"]], layer = "counts.s1")), as.matrix(counts[, 1:4])))
cat("anndata_transfer_object OK\n")

cat("preprocessData.R BRANCH TESTING (not full preprocessing)\n")
src <- parse("docker/scripts/R/preprocessData.R")
targets <- c("saveRDS", "write_anndata", "sctransform_Seurat", "reduceDims_Seurat", "cluster_Seurat", "umap_Seurat")
top <- as.list(src)
sel <- Filter(function(e) is.call(e) && identical(e[[1]], as.name("if")) && any(all.names(e) %in% targets), top)
final_save <- Filter(function(e) is.call(e) && identical(e[[1]], as.name("saveRDS")), top)
stopifnot(length(sel) == 8L, length(final_save) == 1L)

preprocess_branch_test <- function(output_sce, output_seurat, output_anndata, run_sct) {
    env <- new.env(parent = baseenv())
    env$opt <- list(output_sce = output_sce, output_seurat = output_seurat,
        output_anndata = output_anndata, run_sct = run_sct, out = "/tmp/out")
    env$sce <- "FAKE_SCE"
    env$seurat_sce <- "FAKE_SEURAT"
    calls <- list()
    mk <- function(nm) { force(nm); function(...) calls[[length(calls) + 1L]] <<- list(name = nm, args = list(...)) }
    for (nm in targets) env[[nm]] <- mk(nm)
    for (e in sel) if (isTRUE(eval(e[[2]], env))) eval(e[[3]], env)
    for (e in final_save) eval(e, env)
    list(calls = calls)
}

flags <- expand.grid(output_sce = c(FALSE, TRUE), output_seurat = c(FALSE, TRUE),
    output_anndata = c(FALSE, TRUE), run_sct = c(FALSE, TRUE))
stopifnot(nrow(flags) == 16L)
for (i in seq_len(nrow(flags))) {
    f <- flags[i, ]
    res <- preprocess_branch_test(f$output_sce, f$output_seurat, f$output_anndata, f$run_sct)
    save_calls <- Filter(function(c) c$name %in% c("saveRDS", "write_anndata"), res$calls)
    files <- vapply(save_calls, function(c) if (c$name == "saveRDS") c$args$file else c$args[[2]], character(1))
    expected <- c(
        if (f$output_sce) c("/tmp/out_unfiltered_sce.rds.gz", "/tmp/out_filtered_sce.rds.gz"),
        if (f$output_seurat) "/tmp/out_unfiltered_seurat_sce.rds.gz",
        if (f$output_anndata) "/tmp/out_unfiltered.h5ad",
        "/tmp/out_filtered_seurat_sce.rds.gz"
    )
    stopifnot(identical(sort(files), sort(expected)))
    nms <- vapply(res$calls, function(c) c$name, character(1))
    stopifnot(sum(nms == "sctransform_Seurat") == as.integer(f$run_sct))
    sct_assay <- vapply(res$calls, function(c) c$name %in% c("reduceDims_Seurat", "cluster_Seurat", "umap_Seurat") && identical(c$args$assay, "SCT"), logical(1))
    stopifnot(sum(sct_assay) == 3L * as.integer(f$run_sct))
    if (f$run_sct) {
        stopifnot(sum(nms == "reduceDims_Seurat") == 1L, sum(nms == "cluster_Seurat") == 1L, sum(nms == "umap_Seurat") == 1L)
    }
}
cat("preprocessData.R branch testing OK (16/16 combinations)\n")
cat("ALL TESTS PASSED\n")
