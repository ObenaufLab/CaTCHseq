####### FUNCTIONS #########
activity_by_logmeans <- function(signature, sce, fn.scale = scale, ...) {
    # check if logcounts is present
    if (class(sce)[1] == "SingleCellExperiment") {
        if (!("logcounts" %in% names(assays(sce)))) {
            stop()
        }
    } else {
        if (!("data" %in% Layers(sce))) {
            stop()
        }
    }

    # determine shared features between sce and signature
    signature <- unlist(signature)
    shared_features <- intersect(signature, rownames(sce))
    if (length(shared_features) == 0) {
        print("Trying to match case title, no intersect found")
        shared_features <- intersect(str_to_title(signature), rownames(sce))
    }
    if (length(shared_features) == 0) {
        print("Trying to match case lower, no intersect found")
        shared_features <- intersect(str_to_lower(signature), rownames(sce))
    }
    stopifnot(length(shared_features) > 0)

    if (class(sce)[1] == "SingleCellExperiment") {
        # calculate scaled log-means of signature
        res <- sce[shared_features, ] %>%
            assay("logcounts") %>% # use scale by default, but allow other functions and additional arguments
            fn.scale(...)
    } else {
        res <- LayerData(sce, assay = "RNA", layer = "data") %>%
            .[shared_features, ] %>% # use scale by default, but allow other functions and additional arguments
            fn.scale(...)
    }
    if (is.matrix(res)) {
        res <- Matrix::colMeans(res)
    } else {
        res <- mean(res)
    }
    return(res)
}


assignCategoryByMarker <- function(sce, markers = NULL, col.name = NULL, fn.scale = scale, center = TRUE, ...) {
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
    categories[is.na(categories)] <- 0.0 # Replace NANs by 0 otherwise length of intersect error

    if (class(sce)[1] == "SingleCellExperiment") {
        colData(sce)[col.name] <- colnames(categories)[unlist(apply(categories, 1, which.max))]
    } else {
        sce@meta.data[[col.name]] <- colnames(categories)[unlist(apply(categories, 1, which.max))]
    }
    return(sce)
}


# Custom version of progeny that supports sparse matrix as input
progeny.CaTCH <- function(expr, scale = TRUE, organism = "Human", top = 100, verbose = FALSE, ...) {
    if (!is.logical(scale)) {
        stop("scale should be a logical value")
    }

    if (!is.logical(verbose)) {
        stop("verbose should be a logical value")
    }

    model <- getModel(organism, top = top)
    common_genes <- intersect(rownames(expr), rownames(model))

    if (verbose) {
        number_genes <- apply(model, 2, function(x) {
            sum(rownames(model)[which(x != 0)] %in% unique(rownames(expr)))
        })
        message("Number of genes used per pathway to compute progeny scores:")
        message(paste(names(number_genes), ": ", number_genes, " (",
            (number_genes / top) * 100, "%)",
            sep = "", "\n"
        ))
    }

    result <- t(expr[common_genes, , drop = FALSE]) %*%
        as.matrix(model[common_genes, , drop = FALSE])

    if (scale && nrow(result) > 1) {
        rn <- rownames(result)
        result <- apply(result, 2, scale)
        rownames(result) <- rn
    }

    return(result)
}

# create Plots
createValueDistrPlot <- function(sce, grp.col = "Cluster", val.col = "CMO", colors = NULL, grp.order = NULL, xlab = NULL, ylab = NULL, title = NULL) {
    if (is.null(sce)) {
        stop("'sce' must be specified and cannot be NULL")
    }
    if (class(sce)[1] == "SingleCellExperiment") {
        plot.data <- colData(sce) %>%
            tibble::as_tibble() %>%
            dplyr::count(.data[[grp.col]], .data[[val.col]]) %>%
            dplyr::group_by(.data[[grp.col]]) %>%
            dplyr::mutate(Proportion = n / sum(n) * 100)
    } else {
        plot.data <- sce@meta.data %>%
            tibble::as_tibble() %>%
            dplyr::count(.data[[grp.col]], .data[[val.col]]) %>%
            dplyr::group_by(.data[[grp.col]]) %>%
            dplyr::mutate(Proportion = n / sum(n) * 100)
    }

    if (!is.null(grp.order)) {
        tmp <- tibble()
        for (x in grp.order) {
            tmp <- rbind(tmp, plot.data[which(plot.data[[grp.col]] == x), ])
        }
        tmp[[grp.col]] <- factor(tmp[[grp.col]], levels = grp.order)
        plot.data <- tmp
    }
    plot <- ggplot2::ggplot(data = plot.data, ggplot2::aes(fill = .data[[val.col]], x = .data[[grp.col]], y = Proportion)) +
        ggplot2::geom_bar(position = "stack", stat = "identity")
    if (!is.null(colors)) {
        plot <- plot + ggplot2::scale_fill_manual(values = colors)
    }
    if (!is.null(xlab)) {
        plot <- plot + ggplot2::xlab(xlab)
    }
    if (!is.null(ylab)) {
        plot <- plot + ggplot2::ylab(ylab)
    }
    if (!is.null(title)) {
        plot <- plot + ggplot2::ggtitle(title)
    }
    return(plot)
}

# Generates the heatmap based on the provided output of findMarkers
generateHeatmaps <- function(markers,
                             name = "log2FC",
                             max.FDR = NULL,
                             min.FC = NULL,
                             cluster.rows = TRUE,
                             cluster.cols = TRUE,
                             row.title = NULL,
                             col.title = NULL,
                             col.annotation = NULL,
                             row.annotation = NULL) {
    # Reformat the list of markers so that each entry (DataFrame) becomes a tibble with the rownames (gene symbols) stored in the column "symbol"
    # and the item name. i.e. cluster name, in the column "cluster" of the resulting tibble.
    # Finally, all tibbles are merged into a single big tibble.
    stats_tfa <- map_df(as.list(markers), ~ as_tibble(.x, rownames = "symbol"), .id = "cluster") %>%
        arrange(symbol)

    tmp <- stats_tfa
    if (!is.null(max.FDR)) {
        tmp <- tmp %>%
            filter(FDR < max.FDR)
    }
    if (!is.null(min.FC)) {
        tmp <- tmp %>%
            filter(abs(summary.logFC) > min.FC)
    }
    tf_candidates <- tmp %>%
        distinct(symbol) %>%
        pull(symbol)

    stats_tfa <- stats_tfa %>%
        filter(symbol %in% tf_candidates)

    mat_effect <- stats_tfa %>%
        dplyr::select(cluster, symbol, summary.logFC) %>%
        pivot_wider(names_from = cluster, values_from = summary.logFC) %>%
        data.frame(row.names = "symbol", check.names = FALSE) %>%
        as.matrix()
    mat_stats <- stats_tfa %>%
        dplyr::select(cluster, symbol, FDR) %>%
        pivot_wider(names_from = cluster, values_from = FDR) %>%
        data.frame(row.names = "symbol", check.names = FALSE) %>%
        as.matrix()

    stats_value <- 0.05
    mat <- mat_effect
    mat[mat_stats >= stats_value] <- 0
    max_mat <- max(abs(mat_effect), na.rm = TRUE)

    col_pwa <- colorRamp2(breaks = c(-max_mat, 0, max_mat), colors = c("blue", "white", "red"))

    p_hm_raw <- Heatmap(
        matrix = mat,
        name = name,
        border = TRUE,
        cluster_rows = cluster.rows,
        cluster_columns = cluster.cols,
        col = col_pwa,
        column_names_rot = 90,
        column_names_gp = gpar(fontsize = 12),
        column_title = col.title,
        row_names_gp = gpar(fontsize = 10),
        row_title = row.title,
        cell_fun = function(j, i, x, y, width, height, fill) {
            if (mat_stats[i, j] <= stats_value) {
                grid.text(sprintf("%.1f", mat_effect[i, j]), x, y, gp = gpar(fontsize = 8, col = "black"))
            } else {
                grid.text(sprintf("%.1f", mat_effect[i, j]), x, y, gp = gpar(fontsize = 8, col = "grey"))
            }
        },
        top_annotation = col.annotation,
        left_annotation = row.annotation
    )
    return(p_hm_raw)
}


### Create sce objects with corresponding CaTCH barcodes ###
create_SCEs <- function(smpl, data10X, bc, annotation, minBC = 10, singletCut = 0.9, bc1Cut = 0.4, bc2Cut = 0.3) {
    samplelist <- str_split(smpl, ",")[[1]]
    data10Xs <- str_split(data10X, ",")[[1]]
    bcs <- str_split(bc, ",")[[1]]
    anno <- annotation
    minBC <- minBC
    singletCut <- singletCut
    bc1Cut <- bc1Cut
    bc2Cut <- bc2Cut

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

        # Add both gene symbol and gene ID to rowData
        rowData(sce)[c("GeneID", "Symbol")] <- read_tsv(
            file = paste(data10X, "features.tsv.gz", sep = "/"),
            col_names = c("GeneID", "Symbol"),
            col_types = "cc-"
        )

        # Add the biotype and the locus data
        gtf <- read_gtf(annotation)

        rowData(sce)[c("Biotype", "Locus")] <- rowData(sce) %>%
            tibble::as_tibble() %>%
            dplyr::left_join(y = gtf, by = c("GeneID" = "gene_id")) %>%
            dplyr::select(biotype, locus)


        sce$Sample <- smpl
        sce$CellID <- colnames(sce)
        print(paste("   Loaded expression data for ", ncol(sce), " cells", sep = ""))

        #### Now load the CaTCH barcodes ####
        data.catch <- load_BCs(bc)
        print(paste0(" Loading Barcodes from ", bc))

        colData(sce)[, c("CaTCH.BCs", "CaTCH.BC.counts")] <- colData(sce) %>%
            tibble::as_tibble() %>%
            dplyr::left_join(y = data.catch, by = "CellID") %>%
            dplyr::select(CaTCH.BCs, CaTCH.BC.counts)

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
            tibble::as_tibble() %>%
            dplyr::mutate(
                tmp = if_else(is.na(CaTCH.Sum) | (is.na(CaTCH.BC1) & is.na(CaTCH.BC2)) | CaTCH.Sum < minBC,
                    "No_barcode",
                    if_else(CaTCH.BC1 / CaTCH.Sum >= singletCut,
                        "Singlet",
                        if_else(CaTCH.BC1 / CaTCH.Sum >= bc1Cut & CaTCH.BC2 / CaTCH.Sum >= bc2Cut,
                            "Double_Integration",
                            "Multiplet"
                        )
                    )
                )
            ) %>%
            dplyr::select(tmp) %>%
            dplyr::mutate(
                CaTCH.Status = factor(as.character(tmp), levels = c(
                    "Singlet",
                    "Double_Integration",
                    "Multiplet",
                    "No_barcode",
                    "NA"
                ))
            ) %>%
            dplyr::select(CaTCH.Status)

        seurat_sce <- as.Seurat(sce, data = NULL, assay = NULL)
        seurat_sce <- RenameAssays(seurat_sce, assay.name = "originalexp", new.assay.name = "RNA")
        seurat_sce[["RNA"]] <- as(object = seurat_sce[["RNA"]], Class = "Assay5")

        return(list(sce = sce, seurat_sce = seurat_sce))
    } else {
        scetomerge <- list()
        seurattomerge <- list()
        for (i in 1:length(samplelist)) {
            if (i == 1) {
                print("Loading first sce dataset")
                sl <- create_SCEs(samplelist[i], data10Xs[i], bcs[i], anno)
                firstsce <- sl$sce
                firstseurat <- sl$seurat_sce
            } else {
                print(paste0("Loading sce dataset ", i))
                sl <- create_SCEs(samplelist[i], data10Xs[i], bcs[i], anno)
                scetomerge[[i - 1]] <- sl$sce
                seurattomerge[[i - 1]] <- sl$seurat_sce
            }
        }

        print("Merging SCE")
        sce <- firstsce
        for (x in scetomerge) {
            sce <- cbind(sce, x)
        }

        print("Merging Seurat")
        seurat_sce <- merge(firstseurat, seurattomerge, add.cell.ids = samplelist, project = "CaTCHseq", merge.data = TRUE, merge.dr = FALSE) # merge all the seurat datasets
        # Cleanup missing factor levels for seurat object
        seurat_sce[[]]$CaTCH.Status <- factor(seurat_sce[[]]$CaTCH.Status, levels = c("Singlet", "Double_Integration", "Multiplet", "No_barcode", "NA"))

        return(list(sce = sce, seurat_sce = seurat_sce))
    }
}


### Create sce objects with corresponding CaTCH barcodes ###
create_SCE_only <- function(smpl, data10X, bc, annotation) {
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

        # Add both gene symbol and gene ID to rowData
        rowData(sce)[c("GeneID", "Symbol")] <- read_tsv(
            file = paste(data10X, "features.tsv.gz", sep = "/"),
            col_names = c("GeneID", "Symbol"),
            col_types = "cc-"
        )

        # Add the biotype and the locus data
        gtf <- read_gtf(annotation)

        rowData(sce)[c("Biotype", "Locus")] <- rowData(sce) %>%
            tibble::as_tibble() %>%
            dplyr::left_join(y = gtf, by = c("GeneID" = "gene_id")) %>%
            dplyr::select(biotype, locus)


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
            tibble::as_tibble() %>%
            dplyr::left_join(y = data.catch, by = "CellID") %>%
            dplyr::select(CaTCH.BCs, CaTCH.BC.counts)

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
            tibble::as_tibble() %>%
            dplyr::mutate(
                tmp = if_else(is.na(CaTCH.Sum) | (is.na(CaTCH.BC1) & is.na(CaTCH.BC2)) | CaTCH.Sum < minBC,
                    "No_barcode",
                    if_else(CaTCH.BC1 / CaTCH.Sum >= singletCut,
                        "Singlet",
                        if_else(CaTCH.BC1 / CaTCH.Sum >= bc1Cut & CaTCH.BC2 / CaTCH.Sum >= bc2Cut,
                            "Double_Integration",
                            "Multiplet"
                        )
                    )
                )
            ) %>%
            dplyr::select(tmp) %>%
            tibble::as_tibble() %>%
            dplyr::mutate(
                CaTCH.Status = factor(tmp, levels = c(
                    "Singlet",
                    "Double_Integration",
                    "Multiplet",
                    "No_barcode",
                    "NA"
                ))
            ) %>%
            dplyr::select(CaTCH.Status)

        return(sce)
    } else {
        scetomerge <- list()
        for (i in 1:length(samplelist)) {
            if (i == 1) {
                print("Loading first sce dataset")
                firstsce <- create_SCE_only(samplelist[i], data10Xs[i], bcs[i])
            } else {
                print(paste0("Loading sce dataset ", i))
                sce <- create_SCE_only(samplelist[i], data10Xs[i], bcs[i])
                scetomerge[[i - 1]] <- sce
            }
        }

        sce <- cbind(firstsce, scetomerge) # merge all the sce datasets
        return(sce)
    }
}


ggplotColours <- function(n = 6, h = c(0, 360) + 15) {
    if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360 / n
    hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
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


read_gtf <- function(anno) {
    gtf <- read_tsv(anno,
        comment = "#",
        col_names = c("chr", "type", "start", "end", "strand", "attr"),
        col_types = "c-cdd-c-c"
    ) %>%
        filter(type == "gene") %>%
        mutate(
            gene_id = str_extract(attr, 'gene_id "[^"]+"'),
            biotype = str_extract(attr, 'gene_biotype "[^"]+"')
        ) %>%
        dplyr::select(-attr, -type) %>%
        mutate(
            gene_id = str_extract(gene_id, '"[^"]+"'),
            gene_id = str_remove_all(gene_id, "\""),
            biotype = str_extract(biotype, '"[^"]+"'),
            biotype = str_remove_all(biotype, "\""),
            locus = sprintf("%s%s:%d-%d", strand, chr, start, end)
        )
    return(gtf)
}


normalize_Seurat <- function(sce, normalization.method = "LogNormalize", scale.factor = 10000, selection.method = "vst", nfeatures = 2000, ...) {
    print("   Normalizing and scaling SCE ...")

    sce <- NormalizeData(sce, normalization.method = normalization.method, scale.factor = scale.factor)
    sce <- FindVariableFeatures(sce, selection.method = selection.method, nfeatures = nfeatures)
    sce <- ScaleData(sce, ...)

    return(sce)
}

sctransform_Seurat <- function(sce, vst.flavor = "v2", vars.to.regress = "percent.mt", ...) {
    sce <- SCTransform(sce, vst.flavor = vst.flavor, vars.to.regress = vars.to.regress, verbose = FALSE, ...)
    # sce[["SCT"]] <- as(object = sce[["SCT"]], Class = "Assay5")  # this actually breaks downstream code as SCTAssay is still v3 (https://github.com/satijalab/seurat/issues/7542)
    return(sce)
}

reduceDims_Seurat <- function(sce, assay = "RNA", reduction.name = "pca", ...) {
    print("   Reducing dimensions ...")

    sce <- RunPCA(sce, assay = assay, reduction.name = reduction.name, ...)
    return(sce)
}

cluster_Seurat <- function(sce, assay = "RNA_integrated.cca", reduction = "integrated.cca", cluster.name = "integrated.cca_cluster", resolution = .2, nn.method = "annoy", annoy.metric = "euclidean", n.trees = 50, random.seed = 42, algorithm = 1, method = "matrix", n.start = 10, n.iter = 10, group.singletons = TRUE, initial.membership = NULL, ...) {
    print("   Clustering ...")

    sce <- FindNeighbors(sce, assay = assay, reduction = reduction, compute.SNN = TRUE, graph.name = c(paste0(assay, "_nn"), paste0(assay, "_snn")), nn.method = nn.method, annoy.metric = annoy.metric, n.trees = n.trees)
    sce <- FindClusters(sce, resolution = resolution, cluster.name = cluster.name, graph.name = paste0(assay, "_snn"), random.seed = random.seed, algorithm = algorithm, method = method, n.start = n.start, n.iter = n.iter, group.singletons = group.singletons, initial.membership = initial.membership, ...)

    return(sce)
}

integrate_Seurat <- function(sce, assay = "RNA", method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", normalization.method = "LogNormalize", ...) {
    print("   Integrating samples ...")

    sce <- IntegrateLayers(object = sce, assay = assay, method = method, orig.reduction = orig.reduction, new.reduction = new.reduction, normalization.method = normalization.method, verbose = FALSE, ...)
    # re-join layers after integration
    sce[[paste0(assay, "_", new.reduction)]] <- JoinLayers(sce[[assay]])

    return(sce)
}

umap_Seurat <- function(sce, assay = "RNA", reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca", n.neighbors = 25L, min.dist = 0.1, spread = 5, ...) {
    print("   Preparing UMAP ...")

    sce <- RunUMAP(sce, assay = assay, reduction = reduction, dims = dims, reduction.name = reduction.name, n.neighbors = n.neighbors, min.dist = min.dist, spread = spread, ...)
    return(sce)
}

# Split for integrated analysis
split_Seurat <- function(sce, assay = "RNA", by = "Sample") {
    print("   Splitting SCE by Sample ...")

    sce[[assay]] <- split(sce[[assay]], f = by)
    return(sce)
}

# Join for DE
join_Seurat <- function(sce, assay = "RNA", layers = "data", new = "joined_data") {
    print("   Joining SCE layers  ...")

    sce <- JoinLayers(sce, assay = assay, layers = layers, new = new)
    return(sce)
}

##### Loading of the SCE files generated by 'preprocessData.R' #####
# The pipeline runs 'preprocessData.R' once per library, so the downstream DE
# analyses receive one SCE file per sample/condition/replicate and have to
# combine all of them to see every condition of an experiment. The helpers
# below take a single comma separated string of file paths, which is how the
# Nextflow processes hand over their staged input files.
#
# All of them process the files ONE AT A TIME and only keep the small summary
# they need (cell metadata or an aggregated pseudobulk matrix), so peak memory
# stays at roughly the size of the largest single SCE instead of the sum of all
# of them.

split_SCE_paths <- function(sces) {
    paths <- trimws(unlist(strsplit(as.character(sces), ",")))
    paths <- paths[nzchar(paths)]
    if (length(paths) == 0) {
        stop("No SCE files provided")
    }
    missing <- paths[!file.exists(paths)]
    if (length(missing) > 0) {
        stop(paste0("SCE file(s) not found: ", paste(missing, collapse = ", ")))
    }
    paths
}

# Cell metadata of all SCE files in one data frame, restricted to the columns
# shared by all of them. Cell names are made unique across files so cells
# occurring in several libraries stay separate rows. Only the metadata is kept,
# the SCE objects themselves are discarded right after loading.
collect_SCE_metadata <- function(sces) {
    paths <- split_SCE_paths(sces)
    print(sprintf("   Collecting cell metadata of %d SCE file(s) ...", length(paths)))

    mds <- vector("list", length(paths))
    for (i in seq_along(paths)) {
        sce <- loadRDS(paths[i])
        md <- sce@meta.data
        rownames(md) <- paste("sce", i, rownames(md), sep = "_")
        mds[[i]] <- md
        rm(sce)
        gc(verbose = FALSE)
    }

    shared <- Reduce(intersect, lapply(mds, colnames))
    if (length(shared) == 0) {
        stop("The provided SCE files share no metadata columns")
    }
    md <- do.call(rbind, lapply(mds, function(x) x[, shared, drop = FALSE]))
    print(sprintf("   Collected metadata of %d cells from %d SCE file(s).", nrow(md), length(paths)))
    md
}

# 'preprocessData.R' assigns CaTCH.BC_ID per library from the cells of the
# reference condition only. A library holds a single condition, so every library
# of a non-reference condition can never see a reference cell and ends up with
# 'BC_0' everywhere, and the IDs of different libraries would not be comparable
# even if they were assigned. The IDs are therefore rebuilt here, where the
# metadata of ALL libraries is available: barcodes detected in the reference
# condition are ranked by their mean count per reference sample and numbered
# ('BC_<rank>', 'BC*_<rank>' for double integrations), every other barcode keeps
# 'BC_0'. This is the same ranking 'preprocessData.R' intends, just computed
# over the whole run instead of a single library.
build_BC_ID_map <- function(cellmeta, baseCond) {
    needed <- c("CaTCH.Status", "CaTCH.BC_unique", "Condition", "Replicate")
    missing <- setdiff(needed, colnames(cellmeta))
    if (length(missing) > 0) {
        stop(paste0("Cell metadata is missing the column(s): ", paste(missing, collapse = ", ")))
    }

    ranked <- cellmeta %>%
        as_tibble() %>%
        dplyr::filter(
            CaTCH.Status %in% c("Singlet", "Double_Integration"),
            as.character(Condition) == baseCond
        ) %>%
        dplyr::mutate(Replicate = as.character(Replicate)) %>%
        dplyr::count(CaTCH.BC_unique, Replicate, name = "n") %>%
        tidyr::pivot_wider(names_from = Replicate, values_from = n, values_fill = 0)

    refcols <- setdiff(colnames(ranked), "CaTCH.BC_unique")
    if (length(refcols) == 0) {
        stop(paste0("No cells of the reference condition '", baseCond, "' found in the provided SCE files"))
    }

    ranked <- ranked %>%
        dplyr::mutate(.Means = rowMeans(dplyr::across(dplyr::all_of(refcols)))) %>%
        dplyr::filter(.Means > 0) %>%
        dplyr::arrange(dplyr::desc(.Means)) %>%
        tibble::rowid_to_column(".ID") %>%
        dplyr::mutate(CaTCH.BC_ID = ifelse(grepl("+", CaTCH.BC_unique, fixed = TRUE),
            paste0("BC*_", .ID), paste0("BC_", .ID)
        )) %>%
        dplyr::select(CaTCH.BC_unique, CaTCH.BC_ID)

    print(sprintf(
        "   Ranked %d CaTCH barcode(s) of the reference condition '%s'.",
        nrow(ranked), baseCond
    ))
    ranked
}

# Apply a barcode ID map built by 'build_BC_ID_map' to a vector of barcodes.
# Barcodes that are not part of the map, i.e. that were not detected in the
# reference condition, become 'BC_0'.
apply_BC_ID_map <- function(barcodes, bc_map) {
    ids <- bc_map$CaTCH.BC_ID[match(as.character(barcodes), bc_map$CaTCH.BC_unique)]
    ids[is.na(ids)] <- "BC_0"
    factor(ids, levels = str_sort(unique(ids), numeric = TRUE))
}

# Convenience wrapper for the DE scripts: build the map from the pooled metadata
# and write the IDs back into it.
assign_global_BC_IDs <- function(cellmeta, baseCond) {
    bc_map <- build_BC_ID_map(cellmeta, baseCond)
    cellmeta$CaTCH.BC_ID <- apply_BC_ID_map(cellmeta$CaTCH.BC_unique, bc_map)
    print(sprintf(
        "   Assigned %d CaTCH barcode ID(s), %d cell(s) remain 'BC_0'.",
        nrow(bc_map), sum(as.character(cellmeta$CaTCH.BC_ID) == "BC_0")
    ))
    cellmeta
}

# The IDs are assigned by 'assign_barcode_ids.R' for the whole run, so they are
# used as they are. Objects that never went through that step carry 'BC_0' for
# every cell, in which case the IDs are ranked here instead, which keeps the DE
# scripts usable on their own.
ensure_global_BC_IDs <- function(cellmeta, baseCond) {
    ids <- as.character(cellmeta$CaTCH.BC_ID)
    if (!is.null(ids) && any(ids != "BC_0" & !is.na(ids))) {
        print(sprintf(
            "   Using the CaTCH barcode IDs of the provided objects, %d of %d cell(s) carry a ranked ID.",
            sum(ids != "BC_0" & !is.na(ids)), length(ids)
        ))
        return(cellmeta)
    }
    print("   The provided objects carry no ranked CaTCH barcode IDs, ranking them from the reference condition ...")
    assign_global_BC_IDs(cellmeta, baseCond)
}

# Pseudobulk counts of all SCE files. Every file is aggregated on its own and
# only the resulting (genes x groups) matrix is kept, which is orders of
# magnitude smaller than the SCE it came from. Genes are unioned across files,
# missing genes count as 0.
#
# Aggregation groups by Replicate (which is '<Condition>_<replicate>') instead
# of Condition, so replicates of the same condition stay SEPARATE pseudobulk
# samples instead of being summed into one column. Keeping the condition as the
# column name prefix is required by 'run_deseq'/'run_edger', which select the
# samples of a contrast via startsWith(<Condition>_).
# Returns a list with the count matrix and the matching sample metadata.
aggregate_SCE_pseudobulk <- function(sces, assay = "RNA", group.by = "pca_cluster") {
    paths <- split_SCE_paths(sces)
    print(sprintf("   Aggregating pseudobulk counts of %d SCE file(s) ...", length(paths)))

    mats <- vector("list", length(paths))
    metas <- vector("list", length(paths))

    for (i in seq_along(paths)) {
        sce <- loadRDS(paths[i])
        if (DefaultAssay(sce) != assay) {
            DefaultAssay(sce) <- assay
        }
        if (!is.null(group.by) && group.by %in% colnames(sce[[]])) {
            Idents(sce) <- sce[[]][[group.by]]
        }

        # Aggregate over one explicit grouping column. Grouping by several
        # variables is not used here because AggregateExpression drops the ones
        # that are constant within a file, which would lose the replicate.
        sce$pb_group <- paste(as.character(sce$Replicate), as.character(Idents(sce)), sep = "__")
        groups <- unique(as.character(sce$pb_group))

        ps <- AggregateExpression(sce,
            assays = assay, return.seurat = FALSE,
            group.by = "pb_group"
        )[[assay]]

        # Seurat may replace '_' by '-' in the group names, so map the returned
        # column names back to the original grouping values.
        idx <- match(colnames(ps), groups)
        if (anyNA(idx)) {
            idx <- match(colnames(ps), gsub("_", "-", groups))
        }
        if (anyNA(idx)) {
            stop(paste0(
                "Could not map the pseudobulk column(s) back to a Replicate: ",
                paste(colnames(ps)[is.na(idx)], collapse = ", ")
            ))
        }

        parts <- strsplit(groups[idx], "__", fixed = TRUE)
        colrep <- vapply(parts, function(x) x[1], character(1))
        colident <- vapply(parts, function(x) paste(x[-1], collapse = "__"), character(1))
        # '<Condition>_<replicate>_<cluster>' keeps the condition as the prefix.
        colnames(ps) <- paste(colrep, colident, sep = "_")

        # Replicate -> Condition of this library, used to annotate the columns.
        lookup <- unique(data.frame(
            Replicate = as.character(sce$Replicate),
            Condition = as.character(sce$Condition),
            stringsAsFactors = FALSE
        ))

        mats[[i]] <- ps
        metas[[i]] <- data.frame(
            Sample = colnames(ps),
            Condition = lookup$Condition[match(colrep, lookup$Replicate)],
            Replicate = colrep,
            stringsAsFactors = FALSE
        )

        rm(sce, ps)
        gc(verbose = FALSE)
    }

    genes <- sort(Reduce(union, lapply(mats, rownames)))
    countData <- do.call(cbind, lapply(mats, function(m) {
        full <- matrix(0,
            nrow = length(genes), ncol = ncol(m),
            dimnames = list(genes, colnames(m))
        )
        full[rownames(m), ] <- as.matrix(m)
        full
    }))
    metadata <- do.call(rbind, metas)

    # Libraries of the same replicate would collide, keep the columns unique.
    if (anyDuplicated(colnames(countData))) {
        uniq <- make.unique(colnames(countData), sep = ".")
        colnames(countData) <- uniq
        metadata$Sample <- uniq
    }

    print(sprintf(
        "   Aggregated %d genes x %d pseudobulk samples covering %d condition(s): %s",
        nrow(countData), ncol(countData), length(unique(metadata$Condition)),
        paste(sort(unique(metadata$Condition)), collapse = ", ")
    ))

    list(counts = countData, metadata = metadata)
}

run_deseq <- function(contrast, sampleData_all, countData_all, ...) {
    contrast_name <- contrast
    contrast_groups <- strsplit(contrast, "-vs-")
    print(paste("Comparing ", contrast_name, sep = ""))

    # determine contrast
    A <- unlist(strsplit(contrast_groups[[1]][1], "\\+"), use.names = FALSE)
    B <- unlist(strsplit(contrast_groups[[1]][2], "\\+"), use.names = FALSE)

    print(paste("A: ", A, "B: ", B))

    # subset Datasets for pairwise comparison
    countData <- countData_all %>%
        dplyr::select(starts_with(c(paste0(A, "_"), paste0(B, "_")))) %>%
        filter_at(vars(starts_with(c(paste0(A, "_"), paste0(B, "_")))), all_vars(. > 0))
    sampleData <- sampleData_all %>%
        filter(grepl(paste0("^", A, "$"), Condition) | grepl(paste0("^", B, "$"), Condition))
    sampleData$Condition <- factor(sampleData$Condition, levels = unique(sampleData$Condition))
    # 'select(starts_with(...))' above returns the columns of A BEFORE the columns
    # of B, while 'sampleData' keeps the order of the input metadata. DESeq2 and
    # edgeR match the count columns to the sample annotation by position, so the
    # counts are reordered to the sample order here. Without this the annotation
    # of one condition is attached to the counts of the other and the reported
    # fold changes come out inverted.
    countData <- countData %>% dplyr::select(dplyr::all_of(as.character(sampleData$Sample)))
    samples <- sampleData$Sample
    sampleData <- sampleData %>% add_column(type = "none")
    sampleData <- sampleData %>% add_column(batch = "none")

    ## Create design-table considering different types (paired, unpaired) and batches
    if (length(unique(subset(sampleData, A == Condition)$type)) > 1 | length(unique(subset(sampleData, B == Condition)$type)) > 1) {
        if (length(unique(subset(sampleData, A == Condition)$batch)) > 1 | length(unique(subset(sampleData, B == Condition)$batch)) > 1) {
            design <- ~ type + batch + Condition
        } else {
            design <- ~ type + Condition
        }
    } else {
        if (length(unique(subset(sampleData, A == Condition)$batch)) > 1 | length(unique(subset(sampleData, B == Condition)$batch)) > 1) {
            design <- ~ batch + Condition
        } else {
            design <- ~Condition
        }
    }
    print(design)

    # Create DESeqDataSet
    dds <- DESeqDataSetFromMatrix(countData = countData, colData = sampleData, design = design)

    # filter low counts
    # keep <- rowSums(counts(dds)) >= 10
    # dds <- dds[keep, ]

    # drop unused samples
    dds$Condition <- droplevels(dds$Condition)

    # relevel to base condition B
    dds$Condition <- relevel(dds$Condition, ref = B[[1]])

    # run for each pair of conditions
    vsd <- NULL
    dds <- tryCatch(
        {
            DESeq(dds, parallel = FALSE, betaPrior = FALSE, minReplicatesForReplace = Inf)
        },
        error = function(e) {
            print("Error, will run manual testing")
            dds <- estimateSizeFactors(dds)
            dds <- estimateDispersionsGeneEst(dds)
            dispersions(dds) <- mcols(dds)$dispGeneEst
            return(nbinomWaldTest(dds))
        }
    )

    # Now we want to transform the raw discretely distributed counts so that we can do clustering. (Note: when you expect a large treatment effect you should actually set blind=FALSE (see https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).

    rld <- NULL
    if (is.null(vsd)) {
        tryCatch(
            {
                rld <- rlogTransformation(dds, blind = TRUE)
                vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
            },
            error = function(e) {
                rld <- NULL
                vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
            }
        )
    }
    # We also write the normalized counts to file
    if (!is.null(rld)) {
        write.table(as.data.frame(assay(rld)), gzfile(paste("DE", "DESEQ2", contrast_name, "table", "rld.tsv.gz", sep = "_")), sep = "\t", col.names = NA)
    }

    write.table(as.data.frame(assay(vsd)), gzfile(paste("DE", "DESEQ2", contrast_name, "table", "vsd.tsv.gz", sep = "_")), sep = "\t", col.names = NA)

    # initialize empty objects
    res <- ""
    resOrdered <- ""
    res <- results(dds, contrast = c("Condition", A, B), parallel = TRUE, ...)
    resn <- res
    res_shrink <- lfcShrink(dds = dds, coef = paste("Condition", A, "vs", B, sep = "_"), res = res, type = "apeglm")

    countData <- countData %>%
      as_tibble(rownames = "Gene")
    
    res_shrink <- left_join(res_shrink %>%
        as_tibble(rownames = NA) %>%
        rownames_to_column("Gene"), countData, by = "Gene") %>%
        mutate(p.adj = as.numeric(as.character(padj))) %>%
        dplyr::select(-padj) %>%
        group_by(Gene, p.adj) %>%
        summarise(across(everything(), ~ paste(unique(.x[!is.na(.x)]), collapse = ","))) %>%
        ungroup() %>%
        distinct() %>%
        relocate(p.adj, .after = log2FoldChange) %>%
        mutate(log2FoldChange = as.numeric(as.character(log2FoldChange))) %>%
        arrange(desc(log2FoldChange), p.adj)

    # sort and output
    resOrdered <- res_shrink[order(res_shrink$log2FoldChange), ]

    # write the table to a tsv file
    write.table(as.data.frame(resOrdered), gzfile(paste("DE", "DESEQ2", contrast_name, "table", "results.tsv.gz", sep = "_")), sep = "\t", row.names = FALSE, quote = F)

    # Output no shrink
    res <- resn[order(resn$log2FoldChange), ]
    res <- left_join(res %>%
        as_tibble(rownames = NA) %>%
        rownames_to_column("Gene"), countData, by = "Gene") %>%
        group_by(Gene, log2FoldChange) %>%
        summarise(across(everything(), ~ paste(unique(.x[!is.na(.x)]), collapse = ","))) %>%
        ungroup() %>%
        mutate(p.adj = as.numeric(as.character(padj))) %>%
        dplyr::select(-padj) %>%
        distinct() %>%
        relocate(p.adj, .after = pvalue) %>%
        arrange(desc(log2FoldChange), p.adj)

    write.table(as.data.frame(res), gzfile(paste("DE", "DESEQ2", contrast_name, "table", "results_noshrink.tsv.gz", sep = "_")), sep = "\t", row.names = FALSE, quote = F)

    r <- resn %>%
        as_tibble(rownames = "Gene") %>%
        filter(!is.na(padj), padj <= 0.1) %>%
        mutate(Type = if_else(log2FoldChange > 0, "Enriched", "Depleted")) %>%
        mutate(p.adj = as.numeric(as.character(padj))) %>%
        relocate(p.adj, .after = pvalue) %>%
        dplyr::select(-padj) %>%
        arrange(desc(Type), desc(abs(log2FoldChange)))

    rm(res, resn, resOrdered)
    return(list(dds = r, res = res_shrink))
}


run_edger <- function(contrast, sampleData_all, countData_all, bcv = 0.1) {
    # Typical values for the common BCV (square-root-dispersion) for datasets arising from well-controlled experiments are 0.4 for human data, 0.1 for data on genetically identical model organisms or 0.01 for technical replicates
    # https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
    contrast_name <- contrast
    contrast_groups <- strsplit(contrast, "-vs-")
    print(paste("Comparing ", contrast_name, sep = ""))

    # determine contrast
    A <- unlist(strsplit(contrast_groups[[1]][1], "\\+"), use.names = FALSE)
    B <- unlist(strsplit(contrast_groups[[1]][2], "\\+"), use.names = FALSE)

    # subset Datasets for pairwise comparison
    countData <- countData_all %>%
        dplyr::select(starts_with(c(paste0(A, "_"), paste0(B, "_")))) %>%
        filter_at(vars(starts_with(c(paste0(A, "_"), paste0(B, "_")))), all_vars(. > 0))
    sampleData <- sampleData_all %>%
        filter(grepl(paste0("^", A, "$"), Condition) | grepl(paste0("^", B, "$"), Condition))
    sampleData$Condition <- factor(sampleData$Condition, levels = unique(sampleData$Condition))
    sampleData$Condition <- relevel(sampleData$Condition, ref = B[[1]])
    # 'select(starts_with(...))' above returns the columns of A BEFORE the columns
    # of B, while 'sampleData' keeps the order of the input metadata. DESeq2 and
    # edgeR match the count columns to the sample annotation by position, so the
    # counts are reordered to the sample order here. Without this the annotation
    # of one condition is attached to the counts of the other and the reported
    # fold changes come out inverted.
    countData <- countData %>% dplyr::select(dplyr::all_of(as.character(sampleData$Sample)))
    samples <- colnames(countData) 
    degroups <- colnames(countData) %>% str_remove_all(., "_\\d")
    ## name types and levels for design
    bl <- sapply("batch", paste0, levels(sampleData$batch)[1:length(levels(sampleData$batch)) - 1])
    tl <- sapply("type", paste0, levels(sampleData$type)[1:length(levels(sampleData$type)) - 1])

    ## Create design-table considering different types (paired, unpaired) and batches
    if (length(unique(subset(sampleData, A == Condition)$type)) > 1 | length(unique(subset(sampleData, B == Condition)$type)) > 1) {
        if (length(unique(subset(sampleData, A == Condition)$batch)) > 1 | length(unique(subset(sampleData, B == Condition)$batch)) > 1) {
            des <- ~ type + batch + Condition
            design <- model.matrix(des, data = sampleData)
            # colnames(design) <- c(levels(sampleData$condition), tl, bl)
        } else {
            des <- ~ type + Condition
            design <- model.matrix(des, data = sampleData)
            # colnames(design) <- c(levels(condition), tl)
        }
    } else {
        if (length(unique(subset(sampleData, A == Condition)$batch)) > 1 | length(unique(subset(sampleData, B == Condition)$batch)) > 1) {
            des <- ~ batch + Condition
            design <- model.matrix(des, data = sampleData)
            # colnames(design) <- c(levels(sampleData$condition), bl)
        } else {
            des <- ~Condition
            design <- model.matrix(des, data = sampleData)
            # colnames(design) <- levels(sampleData$condition)
        }
    }
    print(design)

    ## create DGEList
    genes <- rownames(countData)
    dge <- DGEList(counts = countData, group = degroups, samples = samples, genes = genes)

    ## filter low counts
    # keep <- filterByExpr(dge)
    # dge <- dge[keep, , keep.lib.sizes = FALSE]

    ## normalize with TMM
    dgen <- calcNormFactors(dge, method = "TMM")

    ## create file normalized table
    tmm <- as.data.frame(cpm(dge)) %>%
      as_tibble(rownames = "Gene")

    write.table(as.data.frame(tmm), gzfile(paste("DE_EDGER", contrast_name, "Normalized.tsv.gz", sep = "_")), sep = "\t", quote = F, row.names = FALSE)
    rm(dgen)

    ## estimate Dispersion, THIS IS SKIPPED AS WE HAVE TO SET BCV MANUALLY WITHOUT REPLICATES
    # dge <- estimateDisp(dge, design, robust = TRUE)
    bcv <- bcv
    qlf <- exactTest(dge, pair = c(B, A),dispersion = bcv^2)

    # create sorted results Tables
    tops <- topTags(qlf, n = nrow(qlf$table), sort.by = "logFC")
    tops <- tops$table

    countData <- countData %>%
      as_tibble(rownames = "Gene")
    
    tops <- left_join(tops %>%
        as_tibble(rownames = NA) %>%
        rownames_to_column("Gene"), countData, by = "Gene") %>%
        group_by(Gene, logFC) %>%
        summarise(across(everything(), ~ paste(unique(.x[!is.na(.x)]), collapse = ","))) %>%
        ungroup() %>%
        mutate(log2FoldChange = as.numeric(as.character(logFC))) %>%
        dplyr::select(-logFC) %>%
        mutate(p.adj = as.numeric(as.character(FDR))) %>%
        dplyr::select(-FDR, -genes) %>%
        distinct() %>%
        arrange(desc(log2FoldChange), p.adj)
    
    write.table(tops, gzfile(paste("DE_EDGER", contrast_name, "resultsLogFCsorted.tsv.gz", sep = "_")), sep = "\t", quote = F, row.names = FALSE)

    tops <- tops %>%
        filter(!is.na(p.adj), p.adj <= 0.1) %>%
        mutate(Type = if_else(log2FoldChange > 0, "Enriched", "Depleted")) %>%
        arrange(desc(Type), desc(abs(log2FoldChange)))

    return(list(dds = qlf, res = tops))
}


run_deseq_bcs <- function(contrast, sampleData_all, countData_all) {
    contrast_name <- contrast
    contrast_groups <- strsplit(contrast, "-vs-")
    print(paste("Comparing ", contrast_name, sep = ""))

    # determine contrast
    A <- unlist(strsplit(contrast_groups[[1]][1], "\\+"), use.names = FALSE)
    B <- unlist(strsplit(contrast_groups[[1]][2], "\\+"), use.names = FALSE)

    # subset Datasets for pairwise comparison
    countData <- countData_all %>%
        dplyr::select(starts_with(c(paste0(A, "_"), paste0(B, "_")))) %>%
        filter_at(vars(starts_with(c(paste0(A, "_"), paste0(B, "_")))), all_vars(. > 0))
    sampleData <- sampleData_all %>%
        filter(grepl(paste0("^", A, "$"), Condition) | grepl(paste0("^", B, "$"), Condition))
    # 'select(starts_with(...))' above returns the columns of A BEFORE the columns
    # of B, while 'sampleData' keeps the order of the input metadata. DESeq2 and
    # edgeR match the count columns to the sample annotation by position, so the
    # counts are reordered to the sample order here. Without this the annotation
    # of one condition is attached to the counts of the other and the reported
    # fold changes come out inverted.
    countData <- countData %>% dplyr::select(dplyr::all_of(as.character(sampleData$Replicate)))
    samples <- rownames(sampleData)
    sampleData <- sampleData %>% add_column(type = "none")
    sampleData <- sampleData %>% add_column(batch = "none")
    # 'factor(unique(x), levels = unique(x))' assigned only the DISTINCT conditions
    # back into the column and R recycled them over the rows, so with more than one
    # replicate per condition the samples were relabelled in an alternating pattern.
    sampleData$Condition <- factor(as.character(sampleData$Condition), levels = unique(as.character(sampleData$Condition)))
    sampleData$Condition <- relevel(sampleData$Condition, ref = B)

    ## Create design-table considering different types (paired, unpaired) and batches
    if (length(unique(subset(sampleData, A == Condition)$type)) > 1 | length(unique(subset(sampleData, B == Condition)$type)) > 1) {
        if (length(unique(subset(sampleData, A == Condition)$batch)) > 1 | length(unique(subset(sampleData, B == Condition)$batch)) > 1) {
            design <- ~ type + batch + Condition
        } else {
            design <- ~ type + Condition
        }
    } else {
        if (length(unique(subset(sampleData, A == Condition)$batch)) > 1 | length(unique(subset(sampleData, B == Condition)$batch)) > 1) {
            design <- ~ batch + Condition
        } else {
            design <- ~Condition
        }
    }
    print(design)

    # Create DESeqDataSet
    dds <- DESeqDataSetFromMatrix(countData = countData, colData = sampleData, design = design)

    # filter low counts
    # keep <- rowSums(counts(dds)) >= 10
    # dds <- dds[keep, ]

    # drop unused samples
    dds$Condition <- droplevels(dds$Condition)

    # relevel to base condition B
    dds$Condition <- relevel(dds$Condition, ref = B[[1]])

    # run for each pair of conditions
    vsd <- NULL
    dds <- tryCatch(
        {
            DESeq(dds, parallel = FALSE, betaPrior = FALSE)
        },
        error = function(e) {
            dds <- estimateSizeFactors(dds)
            dds <- estimateDispersionsGeneEst(dds)
            dispersions(dds) <- mcols(dds)$dispGeneEst
            return(nbinomWaldTest(dds))
        }
    )

    # Now we want to transform the raw discretely distributed counts so that we can do clustering. (Note: when you expect a large treatment effect you should actually set blind=FALSE (see https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).

    rld <- NULL
    if (is.null(vsd)) {
        tryCatch(
            {
                rld <- rlogTransformation(dds, blind = TRUE)
                vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
            },
            error = function(e) {
                rld <- NULL
                vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
            }
        )
    }
    # We also write the normalized counts to file
    if (!is.null(rld)) {
        write.table(as.data.frame(assay(rld)), gzfile(paste("DE", "DESEQ2", contrast_name, "table", "rld.tsv.gz", sep = "_")), sep = "\t", col.names = NA)
    }

    write.table(as.data.frame(assay(vsd)), gzfile(paste("DE", "DESEQ2", contrast_name, "table", "vsd.tsv.gz", sep = "_")), sep = "\t", col.names = NA)

    # initialize empty objects
    res <- ""
    resOrdered <- ""
    res <- results(dds, contrast = c("Condition", A, B), parallel = TRUE)
    resn <- res
    res_shrink <- lfcShrink(dds = dds, coef = paste("Condition", A, "vs", B, sep = "_"), res = res, type = "apeglm")
    
    countData <- countData %>%
      as_tibble(rownames = "CaTCH.BC_ID")
    
    res_shrink <- left_join(res_shrink %>%
        as_tibble(rownames = NA) %>%
        rownames_to_column("CaTCH.BC_ID"), countData, by = "CaTCH.BC_ID") %>%
        mutate(p.adj = as.numeric(as.character(padj))) %>%
        dplyr::select(-padj) %>%
        group_by(CaTCH.BC_ID, p.adj) %>%
        summarise(across(everything(), ~ paste(unique(.x[!is.na(.x)]), collapse = ","))) %>%
        ungroup() %>%
        distinct() %>%
        # 'summarise' pasted every column into a string, so the fold change is
        # turned back into a number here. Otherwise it is sorted alphabetically
        # and EnhancedVolcano rejects the table.
        mutate(log2FoldChange = as.numeric(as.character(log2FoldChange))) %>%
        relocate(p.adj, .after = log2FoldChange) %>%
        arrange(desc(log2FoldChange), p.adj)

        # write the table to a tsv file
    write.table(as.data.frame(res_shrink), gzfile(paste("DE", "DESEQ2", contrast_name, "table", "results.tsv.gz", sep = "_")), sep = "\t", row.names = FALSE, quote = F)

    # Output no shrink
    res <- resn[order(resn$log2FoldChange), ]
    res <- left_join(res %>%
        as_tibble(rownames = NA) %>%
        rownames_to_column("CaTCH.BC_ID"), countData, by = "CaTCH.BC_ID") %>%
        group_by(CaTCH.BC_ID, log2FoldChange) %>%
        summarise(across(everything(), ~ paste(unique(.x[!is.na(.x)]), collapse = ","))) %>%
        ungroup() %>%
        mutate(p.adj = as.numeric(as.character(padj))) %>%
        dplyr::select(-padj) %>%
        distinct() %>%
        arrange(desc(log2FoldChange), p.adj)


    write.table(as.data.frame(res), gzfile(paste("DE", "DESEQ2", contrast_name, "table", "results_noshrink.tsv.gz", sep = "_")), sep = "\t", row.names = FALSE, quote = F)

    r <- results(dds,
        alpha = 0.1,
        contrast = c("Condition", A, B)
    ) %>%
        as_tibble(rownames = NA) %>%
        rownames_to_column("CaTCH.BC_ID") %>%
        filter(!is.na(padj), padj <= 0.1) %>%
        mutate(Type = if_else(log2FoldChange > 0, "Enriched", "Depleted")) %>%
        arrange(desc(Type), desc(abs(log2FoldChange)))

    r <- r %>%
        group_by(CaTCH.BC_ID, log2FoldChange) %>%
        summarise(across(everything(), ~ paste(unique(.x[!is.na(.x)]), collapse = ","))) %>%
        ungroup() %>%
        mutate(p.adj = as.numeric(as.character(padj))) %>%
        dplyr::select(-padj) %>%
        distinct() 

    rm(res, resn, resOrdered)
    return(list(dds = r, res = res_shrink))
}


run_edger_bcs <- function(contrast, sampleData_all, countData_all, bcv = 0.1) {
    # Typical values for the common BCV (square-root-dispersion) for datasets arising from well-controlled experiments are 0.4 for human data, 0.1 for data on genetically identical model organisms or 0.01 for technical replicates
    # https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
    contrast_name <- contrast
    contrast_groups <- strsplit(contrast, "-vs-")
    print(paste("Comparing ", contrast_name, sep = ""))

    # determine contrast
    A <- unlist(strsplit(contrast_groups[[1]][1], "\\+"), use.names = FALSE)
    B <- unlist(strsplit(contrast_groups[[1]][2], "\\+"), use.names = FALSE)

    # subset Datasets for pairwise comparison
    countData <- countData_all %>%
        dplyr::select(starts_with(c(paste0(A, "_"), paste0(B, "_")))) %>%
        filter_at(vars(starts_with(c(paste0(A, "_"), paste0(B, "_")))), all_vars(. > 0))
    sampleData <- sampleData_all %>%
        filter(grepl(paste0("^", A, "$"), Condition) | grepl(paste0("^", B, "$"), Condition)) 
    sampleData$Condition <- factor(sampleData$Condition, levels = unique(sampleData$Condition))
    sampleData$Condition <- relevel(sampleData$Condition, ref = B)
    # 'select(starts_with(...))' above returns the columns of A BEFORE the columns
    # of B, while 'sampleData' keeps the order of the input metadata. DESeq2 and
    # edgeR match the count columns to the sample annotation by position, so the
    # counts are reordered to the sample order here. Without this the annotation
    # of one condition is attached to the counts of the other and the reported
    # fold changes come out inverted.
    countData <- countData %>% dplyr::select(dplyr::all_of(as.character(sampleData$Replicate)))
    samples <- sampleData$Replicate
    degroups <- colnames(countData) %>% str_remove_all(., "_\\d")
    ## name types and levels for design
    bl <- sapply("batch", paste0, levels(sampleData$batch)[1:length(levels(sampleData$batch)) - 1])
    tl <- sapply("type", paste0, levels(sampleData$type)[1:length(levels(sampleData$type)) - 1])

    ## Create design-table considering different types (paired, unpaired) and batches
    if (length(unique(subset(sampleData, A == Condition)$type)) > 1 | length(unique(subset(sampleData, B == Condition)$type)) > 1) {
        if (length(unique(subset(sampleData, A == Condition)$batch)) > 1 | length(unique(subset(sampleData, B == Condition)$batch)) > 1) {
            des <- ~ type + batch + Condition
            design <- model.matrix(des, data = sampleData)
            # colnames(design) <- c(levels(sampleData$condition), tl, bl)
        } else {
            des <- ~ type + Condition
            design <- model.matrix(des, data = sampleData)
            # colnames(design) <- c(levels(condition), tl)
        }
    } else {
        if (length(unique(subset(sampleData, A == Condition)$batch)) > 1 | length(unique(subset(sampleData, B == Condition)$batch)) > 1) {
            des <- ~ batch + Condition
            design <- model.matrix(des, data = sampleData)
            # colnames(design) <- c(levels(sampleData$condition), bl)
        } else {
            des <- ~Condition
            design <- model.matrix(des, data = sampleData)
            # colnames(design) <- levels(sampleData$condition)
        }
    }
    print(design)
    
    ## create DGEList
    genes <- rownames(countData)
    dge <- DGEList(counts = countData, group = degroups, samples = samples, genes = genes)

    ## filter low counts
    # keep <- filterByExpr(dge, min.count = 1)
    # dge <- dge[keep, keep.lib.sizes = FALSE]

    ## normalize with TMM
    dgen <- calcNormFactors(dge, method = "TMM")

    ## create file normalized table
    tmm <- as.data.frame(cpm(dge))
    tmm$CaTCH.BC_ID <- dgen$genes$genes
    tmm <- tmm[c(ncol(tmm), 1:ncol(tmm) - 1)]

    write.table(as.data.frame(tmm), gzfile(paste("DE_EDGER", contrast_name, "Normalized.tsv.gz", sep = "_")), sep = "\t", quote = F, row.names = FALSE)
    rm(dgen)

    ## estimate Dispersion, THIS IS SKIPPED AS WE HAVE TO SET BCV MANUALLY WITHOUT REPLICATES
    # dge <- estimateDisp(dge, design, robust = TRUE)
    bcv <- bcv
    qlf <- exactTest(dge, pair=c(B, A), dispersion = bcv^2)

    # create sorted results Tables
    tops <- topTags(qlf, n = nrow(qlf$table), sort.by = "logFC")
    tops <- tops$table

    countData <- countData %>%
      as_tibble(rownames = "CaTCH.BC_ID")
    
    tops <- left_join(tops %>%
        as_tibble(rownames = NA) %>%
        rownames_to_column("CaTCH.BC_ID"), countData, by = "CaTCH.BC_ID") %>%
        group_by(CaTCH.BC_ID, logFC) %>%
        summarise(across(everything(), ~ paste(unique(.x[!is.na(.x)]), collapse = ","))) %>%
        ungroup() %>%
        mutate(log2FoldChange = as.numeric(as.character(logFC))) %>%
        dplyr::select(-logFC) %>%
        mutate(p.adj = as.numeric(as.character(FDR))) %>%
        dplyr::select(-FDR) %>%
        distinct() %>%
        arrange(desc(log2FoldChange), p.adj)

    write.table(tops, gzfile(paste("DE_EDGER", contrast_name, "resultsLogFCsorted.tsv.gz", sep = "_")), sep = "\t", quote = F, row.names = FALSE)

    tops <- tops %>%
        filter(!is.na(p.adj), p.adj <= 0.1) %>%
        mutate(Type = if_else(log2FoldChange > 0, "Enriched", "Depleted")) %>%
        arrange(desc(Type), desc(abs(log2FoldChange)))

    return(list(dds = qlf, res = tops))
}


celltype_anno_celldex <- function(counts, ref, refname, assay = 1, labels = "main") {
    if (labels == "main") {
        labs <- ref$label.main
    } else if (labels == "fine") {
        labs <- ref$label.fine
    } else {
        labs <- ref$label.ont
    }
    print(paste("Running SingleR with label", labels, sep = " "))

    predcells <- SingleR(
        test = counts, ref = ref, assay.type.test = assay,
        labels = labs
    )

    # Plot Heatmap
    pdf(
        file = paste("SingleR_SCE_Predcells_Heatmap_", refname, "_", labels, ".pdf", sep = ""),
        width = 15, height = 10
    )
    print(plotScoreHeatmap(predcells))
    dev.off()

    # Plot Deltas
    pdf(
        file = paste("SingleR_SCE_Predcells_Deltas_", refname, "_", labels, ".pdf", sep = ""),
        width = 15, height = 10
    )
    print(plotDeltaDistribution(predcells, ncol = 3))
    dev.off()

    return(predcells)
}

celltype_anno_generic <- function(counts, ref, refname, assay = 1, labels = "main") {
    
    print(paste("Running SingleR with ", refname, sep = " "))

    predcells <- SingleR(
        test = counts, ref = str_to_title(ref), assay.type.test = assay,
        labels = str_to_title(labs)
    )

    # Plot Heatmap
    pdf(
        file = paste("SingleR_SCE_Predcells_Heatmap_", refname, "_", labels, ".pdf", sep = ""),
        width = 15, height = 10
    )
    print(plotScoreHeatmap(predcells))
    dev.off()

    # Plot Deltas
    pdf(
        file = paste("SingleR_SCE_Predcells_Deltas_", refname, "_", labels, ".pdf", sep = ""),
        width = 15, height = 10
    )
    print(plotDeltaDistribution(predcells, ncol = 3))
    dev.off()

    return(predcells)
}
