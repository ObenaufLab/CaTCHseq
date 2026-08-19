# barbieQ report helpers
#
# Shared plotting/summarising logic for the barbieQ "Barcode diversity &
# differential enrichment" section of the CaTCHv2 reports. Sourced by both
# R/CaTCHv2_report.Rmd and R/CaTCHv2_report_and_tables.Rmd so the two reports
# stay in sync.
#
# All functions return ggplot objects (or, for the heatmap, a ComplexHeatmap
# object) and do NOT decide between interactive/static rendering; the caller
# wraps them with to_plotly() for HTML or print() for PDF.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# ---------------------------------------------------------------------------
# IO
# ---------------------------------------------------------------------------

bq_read_tsv <- function(path) {
  readr::read_tsv(path, show_col_types = FALSE)
}

# Interactive (searchable / sortable / paged) table for HTML output via DT.
# Returns an htmlwidget which - like plotly - MUST be the auto-printed value of
# the chunk (or gathered into an htmltools::tagList) for its dependencies to be
# embedded. Falls back to NULL if DT is unavailable so callers can use kable.
bq_datatable <- function(df, caption = NULL, page_len = 15) {
  if (!requireNamespace("DT", quietly = TRUE)) return(NULL)
  DT::datatable(
    df,
    caption = caption,
    rownames = FALSE,
    filter = "top",
    extensions = "Buttons",
    options = list(
      pageLength = page_len,
      lengthMenu = list(c(10, 15, 25, 50, 100, -1),
                        c("10", "15", "25", "50", "100", "All")),
      scrollX = TRUE,
      dom = "Blfrtip",
      buttons = c("copy", "csv")
    )
  )
}

# Human-readable diagnostic for a missing/unreadable input file. Reports exactly
# what path was checked and, for symlinks, where they point and whether the
# target is reachable - this is the usual failure mode when a report is rendered
# inside a container and a Nextflow publishDir symlink points at an unmounted
# work directory.
bq_missing_file_msg <- function(kind, path, extra = "") {
  if (!nzchar(path)) {
    return(sprintf("No %s results available%s: no file was provided.\n", kind, extra))
  }
  info <- sprintf("No %s results available%s.\n\nLooked for: `%s`\n",
                  kind, extra, path)
  target  <- tryCatch(Sys.readlink(path), error = function(e) NA_character_)
  is_link <- !is.na(target) && nzchar(target)
  if (isTRUE(is_link)) {
    info <- paste0(info, sprintf("This path is a symlink to: `%s`\n", target),
                   sprintf("Symlink target reachable: %s\n", file.exists(target)))
    info <- paste0(info,
                   "\nIf this is a container run, the symlink target is likely ",
                   "outside the paths mounted into the container. Bind-mount the ",
                   "Nextflow `work/` directory, or publish with `mode: 'copy'`, ",
                   "or pass a real (dereferenced) file path.\n")
  } else {
    info <- paste0(info, sprintf("File exists: %s\n", file.exists(path)))
  }
  info
}

# Detect the effect-size column produced by barbieQ (meanDiff) with back-compat
# aliases for other differential methods, and the adjusted p-value column.
bq_effect_col <- function(tab) {
  cand <- c("meanDiff", "logFC", "log2FC", "estimate", "coefficients")
  hit <- intersect(cand, colnames(tab))
  if (length(hit)) hit[1] else NA_character_
}

bq_padj_col <- function(tab) {
  cand <- c("adj.P.Val", "padj", "FDR", "adj.P", "adjP")
  hit <- intersect(cand, colnames(tab))
  if (length(hit)) hit[1] else NA_character_
}

bq_effect_label <- function(effect_col) {
  if (identical(effect_col, "meanDiff")) {
    "meanDiff (asin\u221a proportion difference)"
  } else {
    effect_col
  }
}

# Column used as the human-readable barcode label in plots. Prefer the short
# stable ID (bc_id, e.g. "bc_123") over the long sequence so labels stay legible;
# the full sequence is kept for hover text.
bq_label_col <- function(tab) {
  if ("bc_id" %in% colnames(tab)) "bc_id" else if ("Barcode" %in% colnames(tab)) "Barcode" else NA_character_
}

# ---------------------------------------------------------------------------
# Diversity
# ---------------------------------------------------------------------------

# Per-sample view: how each index moves across samples, kept in samplesheet
# order and coloured by group so group-level shifts are visible. Faceted by
# index. The same helper serves diversity and dominance -- pass the relevant
# columns and axis label.
bq_diversity_per_sample <- function(divTab,
                                    cols = c("shannon", "simpson"),
                                    ylab = "Diversity") {
  idxCols <- intersect(cols, colnames(divTab))
  if (length(idxCols) == 0) return(NULL)

  hasGroup <- "group" %in% colnames(divTab)
  # Keep the row order of the diversity table (samplesheet order) on the x axis
  # rather than sorting alphabetically; grouping is still conveyed by fill.
  divTab <- divTab %>%
    dplyr::mutate(sample = factor(.data[["sample"]], levels = unique(.data[["sample"]])))

  divLong <- divTab %>%
    tidyr::pivot_longer(dplyr::all_of(idxCols), names_to = "index", values_to = "value")

  aesFill <- if (hasGroup) aes(x = sample, y = value, fill = .data[["group"]]) else
    aes(x = sample, y = value, fill = index)

  p <- ggplot(divLong, aesFill) +
    geom_col() +
    facet_wrap(~index, ncol = 1, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    xlab("") + ylab(ylab) +
    labs(fill = if (hasGroup) "group" else "index")

  if (!hasGroup) p <- p + theme(legend.position = "none")
  p
}

# Per-group summary: distribution of each index within each group (boxplot +
# jitter), showing how an index changes between experimental groups. Shared by
# diversity and dominance -- pass the relevant columns and axis label.
bq_diversity_per_group <- function(divTab,
                                   cols = c("shannon", "simpson"),
                                   ylab = "Diversity") {
  if (!("group" %in% colnames(divTab))) return(NULL)

  idxCols <- intersect(cols, colnames(divTab))
  if (length(idxCols) == 0) return(NULL)
  divLong <- divTab %>%
    tidyr::pivot_longer(dplyr::all_of(idxCols), names_to = "index", values_to = "value")

  ggplot(divLong, aes(x = .data[["group"]], y = value, fill = .data[["group"]])) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.15, height = 0, alpha = 0.6, size = 1) +
    facet_wrap(~index, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "none") +
    xlab("") + ylab(ylab)
}

# ---------------------------------------------------------------------------
# Differential enrichment
# ---------------------------------------------------------------------------

# Restrict the differential-proportion table to a user-selected set of
# comparisons, so a manual report run can focus on a few contrasts out of the
# all-vs-all output. `spec` holds comparison labels exactly as they appear in
# the table's `comparison` column (e.g. "treated_vs_control"), separated by
# commas and/or whitespace; an empty spec keeps every comparison (the default,
# i.e. whatever the samplesheet produced). Unknown labels are reported instead
# of being dropped silently, and the requested order is carried into the plots.
bq_filter_comparisons <- function(dpTab, spec = "") {
  spec <- paste(spec, collapse = ",")
  if (!nzchar(trimws(spec))) return(dpTab)

  if (!("comparison" %in% colnames(dpTab))) {
    warning("Results table has no 'comparison' column - ",
            "ignoring the comparison selection.")
    return(dpTab)
  }

  want <- trimws(unlist(strsplit(spec, "[,[:space:]]+")))
  want <- want[nzchar(want)]
  have <- unique(as.character(dpTab$comparison))

  notFound <- setdiff(want, have)
  if (length(notFound)) {
    warning("Comparison(s) not found: ", paste(notFound, collapse = ", "),
            ". Available: ", paste(have, collapse = ", "))
  }

  keep <- intersect(want, have)
  if (length(keep) == 0) {
    warning("None of the requested comparisons are present - ",
            "keeping all comparisons.")
    return(dpTab)
  }

  out <- dpTab[as.character(dpTab$comparison) %in% keep, , drop = FALSE]
  out$comparison <- factor(as.character(out$comparison), levels = keep)
  out
}

# Enhanced volcano, faceted by comparison (one panel per comparison). x = effect
# size (meanDiff), y = -log10(adjusted p). Points coloured up/down/n.s.; the
# Barcode and stats are exposed as hover text so they remain discoverable under
# to_plotly() (rendered as a single, reliable figure rather than per-tab widgets).
bq_volcano_faceted <- function(dpTab, effect_col = NULL, padj_col = NULL,
                               padj_thresh = 0.05, ncol = 3, n_label = 6) {
  if (is.null(effect_col)) effect_col <- bq_effect_col(dpTab)
  if (is.null(padj_col))   padj_col   <- bq_padj_col(dpTab)
  if (is.na(effect_col) || is.na(padj_col)) return(NULL)
  if (!("comparison" %in% colnames(dpTab))) return(NULL)

  x <- dpTab[[effect_col]]
  padj <- dpTab[[padj_col]]
  sig <- !is.na(padj) & padj < padj_thresh

  df <- dpTab
  df$.status <- dplyr::case_when(
    sig & x > 0 ~ "up",
    sig & x < 0 ~ "down",
    TRUE        ~ "n.s."
  )
  df$.neglog <- -log10(padj)
  df$.effect <- x
  df$.padj   <- padj
  labCol   <- bq_label_col(df)
  df$.id    <- if (!is.na(labCol)) as.character(df[[labCol]]) else NA
  df$.seq   <- if ("Barcode" %in% colnames(df)) as.character(df[["Barcode"]]) else NA

  cols <- c(up = "#c0392b", down = "#2c7fb8", "n.s." = "grey70")

  # Label only the top hits per comparison, and only with the short bc_id so the
  # long sequences do not clutter the panels (sequences stay in the hover text).
  labelDat <- NULL
  if (n_label > 0 && !is.na(labCol) && any(sig)) {
    labelDat <- df[sig, , drop = FALSE]
    labelDat <- labelDat %>%
      dplyr::group_by(.data[["comparison"]]) %>%
      dplyr::slice_max(order_by = .neglog, n = n_label, with_ties = FALSE) %>%
      dplyr::ungroup()
  }

  p <- ggplot(df, aes(x = .effect, y = .neglog, color = .status,
                 text = paste0("ID: ", .id,
                               "\nsequence: ", .seq,
                               "\n", effect_col, ": ", round(.effect, 4),
                               "\n", padj_col, ": ", signif(.padj, 3)))) +
    geom_point(alpha = 0.7, size = 1.3) +
    scale_color_manual(values = cols, name = "") +
    geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey60") +
    facet_wrap(~comparison, ncol = ncol, scales = "free") +
    xlab(bq_effect_label(effect_col)) +
    ylab(paste0("-log10(", padj_col, ")"))

  if (!is.null(labelDat) && nrow(labelDat) &&
      requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      data = labelDat,
      aes(x = .effect, y = .neglog, label = .id),
      inherit.aes = FALSE, size = 2.6, max.overlaps = 20,
      color = "grey20", segment.color = "grey60"
    )
  }
  p
}

# Count of significant up/down barcodes per comparison - a quick cross-
# comparison overview.
bq_sig_barplot <- function(dpTab, effect_col = NULL, padj_col = NULL,
                           padj_thresh = 0.05) {
  if (is.null(effect_col)) effect_col <- bq_effect_col(dpTab)
  if (is.null(padj_col))   padj_col   <- bq_padj_col(dpTab)
  if (is.na(effect_col) || is.na(padj_col) || !("comparison" %in% colnames(dpTab)))
    return(NULL)

  x <- dpTab[[effect_col]]
  padj <- dpTab[[padj_col]]
  sig <- !is.na(padj) & padj < padj_thresh

  df <- dpTab
  df$.status <- dplyr::case_when(
    sig & x > 0 ~ "up",
    sig & x < 0 ~ "down",
    TRUE        ~ "n.s."
  )
  df <- df[df$.status != "n.s.", , drop = FALSE]
  if (!nrow(df)) return(NULL)

  counts <- df %>%
    dplyr::count(.data[["comparison"]], .status, name = "n")

  ggplot(counts, aes(x = .data[["comparison"]], y = n, fill = .status)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c(up = "#c0392b", down = "#2c7fb8"), name = "") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    xlab("") + ylab("Significant barcodes")
}

# Heatmap of the top barcodes (by absolute effect size across comparisons) as a
# barcode x comparison matrix of effect sizes. Returns a ComplexHeatmap object.
bq_top_heatmap <- function(dpTab, effect_col = NULL, n_top = 40) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) return(NULL)
  if (is.null(effect_col)) effect_col <- bq_effect_col(dpTab)
  if (is.na(effect_col) ||
      !all(c("Barcode", "comparison") %in% colnames(dpTab))) return(NULL)

  # Row key stays the sequence (unique); label rows with the short bc_id when
  # available so the heatmap stays legible.
  labCol <- bq_label_col(dpTab)
  idMap <- if (!is.na(labCol) && labCol != "Barcode") {
    stats::setNames(as.character(dpTab[[labCol]]), as.character(dpTab[["Barcode"]]))
  } else NULL

  wide <- dpTab %>%
    dplyr::select(dplyr::all_of(c("Barcode", "comparison", effect_col))) %>%
    tidyr::pivot_wider(names_from = "comparison", values_from = dplyr::all_of(effect_col))

  mat <- as.matrix(wide[, -1, drop = FALSE])
  rownames(mat) <- if (!is.null(idMap)) unname(idMap[wide[["Barcode"]]]) else wide[["Barcode"]]

  ord <- order(apply(abs(mat), 1, max, na.rm = TRUE), decreasing = TRUE)
  mat <- mat[utils::head(ord, n_top), , drop = FALSE]
  if (!nrow(mat)) return(NULL)

  rng <- max(abs(mat), na.rm = TRUE)
  col_fun <- circlize::colorRamp2(c(-rng, 0, rng), c("#2c7fb8", "white", "#c0392b"))

  ComplexHeatmap::Heatmap(
    mat,
    name = effect_col,
    col = col_fun,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 7),
    column_names_gp = grid::gpar(fontsize = 8),
    column_names_rot = 45,
    na_col = "grey90"
  )
}
