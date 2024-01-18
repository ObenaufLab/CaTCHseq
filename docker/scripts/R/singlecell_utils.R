####### FUNCTIONS #########

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


createValueDistrPlot <- function(sce, grp.col = "Cluster", val.col = "CMO", colors = NULL, grp.order = NULL, xlab = NULL, ylab = NULL, title = NULL) {
    if (is.null(sce)) {
        stop("'sce' must be specified and cannot be NULL")
    }

    plot.data <- sce@meta.data %>%
        dplyr::count(.data[[grp.col]], .data[[val.col]]) %>%
        dplyr::group_by(.data[[grp.col]]) %>%
        dplyr::mutate(Proportion = n / sum(n) * 100)
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


## Run DESeq2/Edger ####

run_deseq <- function(contrast, sampleData_all, countData_all, ... ) {
    contrast_name <- contrast
    contrast_groups <- strsplit(contrast, "-vs-")
    print(paste("Comparing ", contrast_name, sep = ""))
    
    # determine contrast
    A <- unlist(strsplit(contrast_groups[[1]][1], "\\+"), use.names = FALSE)
    B <- unlist(strsplit(contrast_groups[[1]][2], "\\+"), use.names = FALSE)
    
    print(paste("A: ",A, "B: ", B))
    
    # subset Datasets for pairwise comparison
    countData <- cbind(countData_all[, grepl(paste("^", B, "_", sep = ""), colnames(countData_all))], countData_all[, grepl(paste("^", A, "_", sep = ""), colnames(countData_all))])
    rownames(countData) <- rownames(countData_all)
    sampleData <- droplevels(rbind(subset(sampleData_all, B == Condition), subset(sampleData_all, A == Condition)))
    
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
    
    res_shrink$Gene <- res_shrink %>%
        as_tibble(rownames = NA) %>%
        rownames_to_column("Gene") %>%
        pull(Gene)
    
    res_shrink <- res_shrink %>%
        as_tibble(rownames = NA) %>%
        mutate(p.adj=as.numeric(as.character(padj))) %>%
        select(-padj)
    
    # sort and output
    resOrdered <- res_shrink[order(res_shrink$log2FoldChange), ]
    
    # write the table to a tsv file
    write.table(as.data.frame(resOrdered), gzfile(paste("DE", "DESEQ2", contrast_name, "table", "results.tsv.gz", sep = "_")), sep = "\t", row.names = FALSE, quote = F)
    
    # Output no shrink
    res <- resn[order(resn$log2FoldChange), ]
    res$Gene <- res %>%
        as_tibble(rownames = NA) %>%
        rownames_to_column("Gene") %>%
        pull(Gene)
    
    write.table(as.data.frame(res), gzfile(paste("DE", "DESEQ2", contrast_name, "table", "results_noshrink.tsv.gz", sep = "_")), sep = "\t", row.names = FALSE, quote = F)
    
    r <- resn %>%
        as_tibble(rownames = NA) %>%
        filter(!is.na(padj), padj <= 0.1) %>%
        mutate(Type = if_else(log2FoldChange > 0, "Enriched", "Depleted")) %>%
        arrange(desc(Type), desc(abs(log2FoldChange))) 
    
    r <- r %>%
        mutate(p.adj=as.numeric(as.character(padj))) %>%
        select(-padj)
    
    r$Gene <- r %>%
        as_tibble(rownames = NA) %>%
        rownames_to_column("Gene") %>%
        pull(Gene)
    
    rm(res, resn, resOrdered)
    return(list(dds = r, res = res_shrink))
}


run_edger <- function(contrast, sampleData_all, countData_all, bcv = 0.1) {
    #Typical values for the common BCV (square-root-dispersion) for datasets arising from well-controlled experiments are 0.4 for human data, 0.1 for data on genetically identical model organisms or 0.01 for technical replicates
    #https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
    contrast_name <- contrast
    contrast_groups <- strsplit(contrast, "-vs-")
    print(paste("Comparing ", contrast_name, sep = ""))
    
    # determine contrast
    A <- unlist(strsplit(contrast_groups[[1]][1], "\\+"), use.names = FALSE)
    B <- unlist(strsplit(contrast_groups[[1]][2], "\\+"), use.names = FALSE)
    
    # subset Datasets for pairwise comparison
    countData <- cbind(countData_all[, grepl(paste("^", B, "_", sep = ""), colnames(countData_all))], countData_all[, grepl(paste("^", A, "_", sep = ""), colnames(countData_all))])
    rownames(countData) <- rownames(countData_all)
    sampleData <- droplevels(rbind(subset(sampleData_all, B == Condition), subset(sampleData_all, A == Condition)))
    
    samples <- rownames(sampleData)
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
    dge <- DGEList(counts = countData, group = sampleData$Condition, samples = samples, genes = genes)
    
    ## filter low counts
    #keep <- filterByExpr(dge)
    #dge <- dge[keep, , keep.lib.sizes = FALSE]
    
    # relevel to base condition B
    dge$samples$group <- relevel(dge$samples$group, ref = B[[1]])
    
    ## estimate Dispersion, THIS IS SKIPPED AS WE HAVE TO SET BCV MANUALLY WITHOUT REPLICATES
    #dge <- estimateDisp(dge, design, robust = TRUE)
    bcv <- bcv
    qlf <- exactTest(dge, dispersion=bcv^2)
    
    # create sorted results Tables
    tops <- topTags(qlf, n = nrow(qlf$table), sort.by = "logFC")
    tops <- tops$table
    tops$Gene <- tops %>%
        as_tibble(rownames = NA) %>%
        rownames_to_column("Gene") %>%
        pull(Gene)
    
    tops <- tops %>% 
        mutate(log2FoldChange=as.numeric(as.character(logFC))) %>%
        select(-logFC) %>%
        mutate(p.adj=as.numeric(as.character(FDR))) %>%
        select(-FDR)
    
    tops <- tops  %>%
        filter(!is.na(p.adj), p.adj <= 0.1) %>%
        mutate(Type = if_else(log2FoldChange > 0, "Enriched", "Depleted")) %>%
        arrange(desc(Type), desc(abs(log2FoldChange)))
    
    write.table(tops, gzfile(paste("DE_EDGER", contrast_name, "resultsLogFCsorted.tsv.gz", sep = "_")), sep = "\t", quote = F, row.names = FALSE)
    
    
    return(list(dds = qlf, res = tops))
}


run_deseq_bcs <- function(contrast, sampleData_all, countData_all, ids) {
    contrast_name <- contrast
    contrast_groups <- strsplit(contrast, "-vs-")
    print(paste("Comparing ", contrast_name, sep = ""))
    
    # determine contrast
    A <- unlist(strsplit(contrast_groups[[1]][1], "\\+"), use.names = FALSE)
    B <- unlist(strsplit(contrast_groups[[1]][2], "\\+"), use.names = FALSE)
    
    # subset Datasets for pairwise comparison
    countData <- cbind(countData_all[, grepl(paste("^", B, "_", sep = ""), colnames(countData_all))], countData_all[, grepl(paste("^", A, "_", sep = ""), colnames(countData_all))])
    sampleData <- droplevels(rbind(subset(sampleData_all, B == Condition), subset(sampleData_all, A == Condition)))
    rownames(countData) <- rownames(countData_all)
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
    
    res_shrink <- left_join(res_shrink %>%
                                as_tibble(rownames = NA) %>%
                                rownames_to_column("CaTCH.BC_ID"), ids, by = "CaTCH.BC_ID") %>%
        mutate(p.adj = as.numeric(as.character(padj))) %>%
        select(-padj)
    
    # sort and output
    resOrdered <- res_shrink[order(res_shrink$log2FoldChange), ]
    
    # write the table to a tsv file
    write.table(as.data.frame(resOrdered), gzfile(paste("DE", "DESEQ2", contrast_name, "table", "results.tsv.gz", sep = "_")), sep = "\t", row.names = FALSE, quote = F)
    
    # Output no shrink
    res <- resn[order(resn$log2FoldChange), ]
    res <- left_join(res %>%
                         as_tibble(rownames = NA) %>%
                         rownames_to_column("CaTCH.BC_ID"), ids, by = "CaTCH.BC_ID")
    
    write.table(as.data.frame(res), gzfile(paste("DE", "DESEQ2", contrast_name, "table", "results_noshrink.tsv.gz", sep = "_")), sep = "\t", row.names = FALSE, quote = F)
    
    r <- results(dds,
                 alpha = 0.1,
                 contrast = c("Condition", t, ref.Condition)
    ) %>%
        as_tibble(rownames = NA) %>%
        rownames_to_column("CaTCH.BC_ID") %>%
        filter(!is.na(padj), padj <= 0.1) %>%
        mutate(Type = if_else(log2FoldChange > 0, "Enriched", "Depleted")) %>%
        arrange(desc(Type), desc(abs(log2FoldChange))) 
    
    r <- r %>%
        mutate(p.adj=as.numeric(as.character(padj))) %>%
        select(-padj)
    
    rm(res, resn, resOrdered)
    return(list(dds = r, res = res_shrink))
}


run_edger_bcs <- function(contrast, sampleData_all, countData_all, ids, bcv = 0.1) {
    # Typical values for the common BCV (square-root-dispersion) for datasets arising from well-controlled experiments are 0.4 for human data, 0.1 for data on genetically identical model organisms or 0.01 for technical replicates
    # https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
    contrast_name <- contrast
    contrast_groups <- strsplit(contrast, "-vs-")
    print(paste("Comparing ", contrast_name, sep = ""))
    
    # determine contrast
    A <- unlist(strsplit(contrast_groups[[1]][1], "\\+"), use.names = FALSE)
    B <- unlist(strsplit(contrast_groups[[1]][2], "\\+"), use.names = FALSE)
    
    # subset Datasets for pairwise comparison
    countData <- cbind(countData_all[, grepl(paste("^", B, "_", sep = ""), colnames(countData_all))], countData_all[, grepl(paste("^", A, "_", sep = ""), colnames(countData_all))])
    rownames(countData) <- rownames(countData_all)
    sampleData <- droplevels(rbind(subset(sampleData_all, B == Condition), subset(sampleData_all, A == Condition)))
    
    samples <- rownames(sampleData)
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
    dge <- DGEList(counts = countData, group = sampleData$Condition, samples = samples, genes = genes)
    
    ## filter low counts
    # keep <- filterByExpr(dge)
    # dge <- dge[keep, , keep.lib.sizes = FALSE]
    
    # relevel to base condition B
    dge$samples$group <- relevel(dge$samples$group, ref = B[[1]])
    
    ## estimate Dispersion, THIS IS SKIPPED AS WE HAVE TO SET BCV MANUALLY WITHOUT REPLICATES
    # dge <- estimateDisp(dge, design, robust = TRUE)
    bcv <- bcv
    qlf <- exactTest(dge, dispersion = bcv^2)
    
    # create sorted results Tables
    tops <- topTags(qlf, n = nrow(qlf$table), sort.by = "logFC")
    tops <- tops$table
    tops <- left_join(tops %>%
                          as_tibble(rownames = NA) %>%
                          rownames_to_column("CaTCH.BC_ID"), ids, by = "CaTCHBC_ID") %>%
        mutate(log2FoldChange = as.numeric(as.character(logFC))) %>%
        select(-logFC) %>%
        mutate(p.adj = as.numeric(as.character(FDR))) %>%
        select(-FDR)
    
    write.table(tops, gzfile(paste("DE_EDGER", contrast_name, "resultsLogFCsorted.tsv.gz", sep = "_")), sep = "\t", quote = F, row.names = FALSE)
    
    tops <- tops %>%
        filter(!is.na(p.adj), p.adj <= 0.1) %>%
        mutate(Type = if_else(log2FoldChange > 0, "Enriched", "Depleted")) %>%
        arrange(desc(Type), desc(abs(log2FoldChange)))
    
    return(list(dds = qlf, res = tops))
}


celltype_anno <- function(counts, ref, refname, assay=1, labels='main'){

	if (labels == 'main'){
	   labs = ref$label.main
    }else if(labels == 'fine'){
	   labs = ref$label.fine	
	}else{
	   labs = ref$label.ont	
	}
    print(paste("Running SingleR with label", labels, sep=" "))
        
    predcells <- SingleR(test = counts, ref = ref, assay.type.test = assay,
    labels = labs)

    # Plot Heatmap
    pdf(
        file = paste("SingleR_SCE_Predcells_Heatmap_", refname, "_", labels, ".pdf", sep=""),
        width = 15, height = 10
    )
    print(plotScoreHeatmap(predcells))
	dev.off()

    # Plot Deltas
    pdf(
        file = paste("SingleR_SCE_Predcells_Deltas_", refname, "_", labels, ".pdf", sep=""),
        width = 15, height = 10
    )
    print(plotDeltaDistribution(predcells, ncol = 3))
    dev.off()
    
    return(predcells)
}
