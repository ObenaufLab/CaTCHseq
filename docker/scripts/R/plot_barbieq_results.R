.libPaths("")

library(optparse)

option_list <- list(
    make_option(c("--diversity"),
        type = "character", default = NULL,
        help = "barbieQ_diversity.tsv written by 'barbieQanalysis.R'"
    ),
    make_option(c("--diffprop"),
        type = "character", default = NULL,
        help = "barbieQ_diffProp_results.tsv written by 'barbieQanalysis.R'"
    ),
    make_option(c("--pcut"),
        type = "numeric", default = 0.1,
        help = "Cutoff for adjusted pvalues, default(0.1)"
    ),
    make_option(c("--format"),
        type = "character", default = "pdf",
        help = "format of the plots: pdf (default), jpeg, png, tiff"
    ),
    make_option(c("--width"),
        type = "numeric", default = 20,
        help = "width of the plot. The units depend on the format: inches for PDF, pixels otherwise"
    ),
    make_option(c("--height"),
        type = "numeric", default = 25,
        help = "height of the plot. The units depend on the format: inches for PDF, pixels otherwise"
    ),
    make_option(c("--plots_per_row"),
        type = "numeric", default = 3,
        help = "number of volcano panels per row of the output plot"
    ),
    make_option(c("--out"),
        type = "character", default = NULL,
        help = "prefix for the output file"
    ),
    make_option(c("--libpath"),
        type = "character", default = "/tools/scripts/R/",
        help = "path to R libs, trailing slash is required"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$out)) {
    print_help(opt_parser)
    stop("Not enough parameters")
}

opt$format <- tolower(opt$format)
if (!(opt$format %in% c("pdf", "png", "jpeg", "tiff"))) {
    opt$format <- "pdf"
}

#### Source Functions ####
source(paste0(opt$libpath, "barbieq_report_helpers.R"))

########################################################
library(readr)

open_device <- function(file, width, height, format) {
    if (format == "pdf") {
        pdf(file = file, width = width, height = height)
    } else if (format == "png") {
        png(filename = file, width = width, height = height)
    } else if (format == "jpeg") {
        jpeg(filename = file, width = width, height = height)
    } else {
        tiff(filename = file, width = width, height = height)
    }
}

read_if_present <- function(path) {
    if (is.null(path) || !nzchar(path) || !file.exists(path)) {
        return(NULL)
    }
    tryCatch(bq_read_tsv(path), error = function(e) NULL)
}

divTab <- read_if_present(opt$diversity)
dpTab <- read_if_present(opt$diffprop)

if (is.null(divTab) && is.null(dpTab)) {
    print("No barbieQ result tables found, skipping plots")
    quit(status = 0)
}

#### Diversity and dominance plots ####
if (!is.null(divTab) && nrow(divTab) > 0) {
    plots <- list(
        bq_diversity_per_sample(divTab, cols = c("shannon", "simpson"), ylab = "Diversity"),
        bq_diversity_per_group(divTab, cols = c("shannon", "simpson"), ylab = "Diversity"),
        bq_diversity_per_sample(divTab,
            cols = c("simpson_dominance", "berger_parker"), ylab = "Dominance"
        ),
        bq_diversity_per_group(divTab,
            cols = c("simpson_dominance", "berger_parker"), ylab = "Dominance"
        )
    )
    plots <- plots[!vapply(plots, is.null, logical(1))]

    if (length(plots) > 0) {
        open_device(
            paste0(opt$out, "_barbieQ_diversity.", opt$format),
            opt$width, opt$height, opt$format
        )
        for (p in plots) print(p)
        dev.off()
    }
}

#### Differential enrichment plots ####
if (!is.null(dpTab) && nrow(dpTab) > 0) {
    volcano <- bq_volcano_faceted(dpTab, padj_thresh = opt$pcut, ncol = opt$plots_per_row)
    if (!is.null(volcano) && identical(bq_effect_col(dpTab), "meanDiff")) {
        # ASCII axis label, the unicode label of the helper is not renderable
        # in the C locale of the container
        volcano <- volcano + ggplot2::xlab("meanDiff (asin-sqrt proportion difference)")
    }

    plots <- list(
        volcano,
        bq_sig_barplot(dpTab, padj_thresh = opt$pcut)
    )
    plots <- plots[!vapply(plots, is.null, logical(1))]

    heatmap <- tryCatch(bq_top_heatmap(dpTab), error = function(e) NULL)

    if (length(plots) > 0 || !is.null(heatmap)) {
        open_device(
            paste0(opt$out, "_barbieQ_diffProp.", opt$format),
            opt$width, opt$height, opt$format
        )
        for (p in plots) print(p)
        if (!is.null(heatmap)) print(heatmap)
        dev.off()
    }
}

print("barbieQ plots done")
