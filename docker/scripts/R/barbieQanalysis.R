#!/usr/bin/env Rscript

# MIT License
#
# Copyright (c) 2026 Joerg Fallmann
#
# Barcode diversity (Shannon / normalized Simpson) and differential enrichment
# (barbieQ) analysis for the CaTCHv2 pipeline.
#
# Runs downstream of barcode collation. It builds a barbieQ (SummarizedExperiment)
# object from the collated count table, computes per-sample diversity indices and,
# if the samplesheet defines replicated groups, performs pairwise differential
# proportion testing of every group against the control/baseline group.

suppressWarnings(suppressMessages({
  library(getopt)
  library(tidyverse)
  library(duckplyr)
}))

spec <- matrix(c(
  "help"        , "h", 0, "logical",   "Help",
  "countsFile"  , "c", 1, "character", "Collated barcode count table (CaTCHv2_collation.txt).",
  "designFile"  , "d", 1, "character", "Samplesheet / design file (csv) with sample,group,control columns.",
  "outputDir"   , "o", 1, "character", "Directory in which to write results.",
  "prefix"      , "p", 1, "character", "Prefix for output files.",
  "minReplicates", "r", 1, "integer",  "Minimum replicates in a group to run differential testing (default 2).",
  "threads"      , "t", 1, "integer",  "Parallel threads for differential testing (default 1).",
  "minCount"     , "m", 1, "integer",  "Per comparison, drop barcodes whose summed count within the two contrasted groups is below this (default 10; 0 disables).",
  "rdsIn"        , "R", 1, "character", "Optional precomputed barbieQ .rds to reuse instead of rebuilding the object.",
  "allVsAll"     , "A", 0, "logical",   "Also test every group against every other group (all pairwise). Written to barbieQ_diffProp_results_allVsAll.tsv; the report still uses the baseline-only barbieQ_diffProp_results.tsv."
), byrow = TRUE, ncol = 5)

opt <- getopt(spec)

if (!is.null(opt$help)) {
  cat(getopt(spec, usage = TRUE))
  q(status = 1)
}

if (is.null(opt$countsFile) || is.null(opt$designFile)) {
  stop("Missing input. Need --countsFile and --designFile.")
}
if (is.null(opt$outputDir)) opt$outputDir <- "."
if (is.null(opt$prefix)) opt$prefix <- ""
if (is.null(opt$minReplicates)) opt$minReplicates <- 2
if (is.null(opt$threads) || is.na(opt$threads) || opt$threads < 1) opt$threads <- 1
if (is.null(opt$minCount) || is.na(opt$minCount) || opt$minCount < 0) opt$minCount <- 10
if (is.null(opt$allVsAll)) opt$allVsAll <- FALSE

dir.create(opt$outputDir, recursive = TRUE, showWarnings = FALSE)
out_path <- function(name) file.path(opt$outputDir, paste0(opt$prefix, name))

status_file <- out_path("barbieQ_status.txt")
write_status <- function(msg) {
  cat(msg, "\n", sep = "")
  writeLines(msg, status_file)
}

# Timestamped progress logging to stderr so Nextflow captures it in .command.err
# and the user can follow long-running steps.
verbose_log <- function(msg) {
  message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg)
}

verbose_log("===== CaTCHv2 barbieQ analysis =====")
verbose_log(paste0("Threads requested: ", opt$threads))

## ---------------------------------------------------------------------------
## Optional runtime bootstrap of barbieQ (safety net; normally supplied by the
## conda env conf/barbieq_env.yml or the catchvis container)
## ---------------------------------------------------------------------------
# When CATCHV2_BARBIEQ_NO_BOOTSTRAP is set, never attempt a network install
# (the package is expected to be provided by the conda env / container). This
# also keeps automated tests offline and fast.
no_bootstrap <- nzchar(Sys.getenv("CATCHV2_BARBIEQ_NO_BOOTSTRAP"))
ensure_pkg <- function(pkg, bioc = FALSE) {
  if (requireNamespace(pkg, quietly = TRUE)) return(TRUE)
  if (bioc && !no_bootstrap) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      try(install.packages("BiocManager", repos = "https://cloud.r-project.org"), silent = TRUE)
    }
    if (requireNamespace("BiocManager", quietly = TRUE)) {
      try(BiocManager::install(pkg, update = FALSE, ask = FALSE), silent = TRUE)
    }
  }
  requireNamespace(pkg, quietly = TRUE)
}

have_barbieQ <- ensure_pkg("barbieQ", bioc = TRUE)
have_vegan   <- ensure_pkg("vegan", bioc = FALSE)

## ---------------------------------------------------------------------------
## Read inputs
## ---------------------------------------------------------------------------
counts <- read_delim(opt$countsFile, skip_empty_rows = T, trim_ws = T)
barcodeCol <- colnames(counts)[1]
dropCols <- c(barcodeCol, "bc_id")
sampleCols <- setdiff(colnames(counts), dropCols)

countMatrix <- as.matrix(counts[, sampleCols, drop = FALSE])
mode(countMatrix) <- "numeric"
countMatrix[is.na(countMatrix)] <- 0
rownames(countMatrix) <- counts[[barcodeCol]]

# Map barcode sequence -> short stable ID (bc_id from the collation) so the
# result tables can carry the readable ID alongside the full sequence. The
# sequences remain the row keys; bc_id is added as an extra column downstream.
bcIdMap <- if ("bc_id" %in% colnames(counts)) {
  stats::setNames(as.character(counts[["bc_id"]]), counts[[barcodeCol]])
} else {
  NULL
}

# Detect the samplesheet delimiter: nf-core uses ',' by default but some
# sheets use ';'. Pick whichever appears more often on the header line.
headerLine <- readLines(opt$designFile, n = 1)
nSemi  <- lengths(regmatches(headerLine, gregexpr(";", headerLine)))
nComma <- lengths(regmatches(headerLine, gregexpr(",", headerLine)))
designSep <- if (nSemi > nComma) ";" else ","
design <- read.csv(opt$designFile, header = TRUE, stringsAsFactors = FALSE, sep = designSep)
if (!"sample" %in% colnames(design)) {
  stop("Design file must contain a 'sample' column.")
}
# Collapse technical replicate rows (same sample appears multiple times).
design <- design[!duplicated(design$sample), , drop = FALSE]
rownames(design) <- design$sample

verbose_log(paste0("Inputs loaded: ", nrow(countMatrix), " barcodes, ",
                   length(sampleCols), " count columns, ",
                   nrow(design), " samples in design."))

## ---------------------------------------------------------------------------
## Diversity indices (always computed when vegan is available)
## ---------------------------------------------------------------------------
if (have_vegan) {
  # Compute Shannon and normalized (Gini) Simpson for one count matrix.
  # Samples with zero total barcode counts otherwise get a misleading
  # shannon = 0 / simpson = 1: their proportions are NaN and vegan's internal
  # na.rm collapses the sums to 0. Report such empty samples as NA instead.
  computeDiversity <- function(mat) {
    x <- t(mat)  # vegan expects samples in rows, barcodes in columns
    shannon <- vegan::diversity(x, index = "shannon")
    simpson <- vegan::diversity(x, index = "simpson")
    # Dominance measures (complementary to diversity):
    #  - simpson_dominance: Simpson's D = sum(p_i^2) = 1 - Gini-Simpson. Higher
    #    values mean a few barcodes dominate. Derived from vegan's simpson.
    #  - berger_parker: proportion of the single most abundant barcode
    #    (max(count) / total), the most intuitive dominance index.
    simpson_dominance <- 1 - simpson
    totals <- rowSums(x)
    berger_parker <- apply(x, 1, max) / totals
    emptySample <- totals == 0
    shannon[emptySample] <- NA_real_
    simpson[emptySample] <- NA_real_
    simpson_dominance[emptySample] <- NA_real_
    berger_parker[emptySample] <- NA_real_
    df <- data.frame(
      sample  = rownames(x),
      shannon = shannon,
      simpson = simpson,
      simpson_dominance = simpson_dominance,
      berger_parker = berger_parker,
      stringsAsFactors = FALSE
    )
    if ("group" %in% colnames(design)) {
      df$group <- design[df$sample, "group"]
    }
    # Preserve samplesheet order (design row order) rather than the counts-file
    # column order; samples absent from the design are appended at the end.
    df <- df[order(match(df$sample, design$sample)), , drop = FALSE]
    rownames(df) <- NULL
    attr(df, "nEmpty") <- sum(emptySample)
    df
  }

  # Raw indices: computed on the full, unfiltered barcode matrix. This includes
  # the long tail of low-count/noise barcodes and therefore tends to overstate
  # diversity.
  diversityRaw <- computeDiversity(countMatrix)
  if (attr(diversityRaw, "nEmpty") > 0) {
    verbose_log(paste0("Note: ", attr(diversityRaw, "nEmpty"),
                       " sample(s) had zero barcode counts; diversity set to NA."))
  }
  write.table(diversityRaw, out_path("barbieQ_diversity_raw.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # Filtered indices: drop barcodes whose TOTAL count across all samples is
  # below --minCount before computing diversity, removing the low-count noise
  # tail that artificially inflates the raw indices. This is the default output.
  if (opt$minCount > 0) {
    keepBc <- rowSums(countMatrix) >= opt$minCount
  } else {
    keepBc <- rowSums(countMatrix) > 0
  }
  verbose_log(paste0("Diversity filter: kept ", sum(keepBc), " of ",
                     length(keepBc), " barcodes (global count ",
                     if (opt$minCount > 0) paste0(">= ", opt$minCount) else "> 0",
                     ") for the filtered indices."))
  diversityFiltered <- computeDiversity(countMatrix[keepBc, , drop = FALSE])
  write.table(diversityFiltered, out_path("barbieQ_diversity.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  verbose_log("Diversity indices (Shannon, Simpson) and dominance (Simpson's D, Berger-Parker) computed and saved (raw + filtered).")
} else {
  verbose_log("Package 'vegan' not available - skipping diversity indices.")
}

## ---------------------------------------------------------------------------
## Low-abundance filter mode.
## The filter is applied PER COMPARISON (see run_single_comparison below) over
## only the samples of the two groups actually being contrasted -- never over
## the full matrix. With many conditions a barcode that is near-zero in the two
## tested groups could otherwise survive on counts contributed by unrelated
## conditions. The filter statistic (summed count within the contrasted groups)
## is independent of which of the two groups a sample belongs to, so it does not
## bias the differential proportion model.
## ---------------------------------------------------------------------------
if (opt$minCount > 0) {
  verbose_log(paste0("minCount filter active (>= ", opt$minCount,
                     " summed counts within each contrasted group pair; ",
                     "applied per comparison)."))
} else {
  verbose_log("minCount filter disabled (--minCount 0).")
}

## ---------------------------------------------------------------------------
## Guard: only run differential enrichment with grouped, replicated samples
## ---------------------------------------------------------------------------
run_diff <- TRUE
reason <- ""

if (!have_barbieQ) {
  run_diff <- FALSE; reason <- "barbieQ package not available"
} else if (!"group" %in% colnames(design)) {
  run_diff <- FALSE; reason <- "no 'group' column in design file"
} else {
  # restrict to samples present in the count table
  commonSamples <- intersect(colnames(countMatrix), design$sample)
  design <- design[design$sample %in% commonSamples, , drop = FALSE]
  cm <- countMatrix[, commonSamples, drop = FALSE]
  groups <- design[commonSamples, "group"]
  groupSizes <- table(groups)
  if (length(groupSizes) < 2) {
    run_diff <- FALSE; reason <- "fewer than 2 groups"
  } else if (max(groupSizes) < opt$minReplicates) {
    run_diff <- FALSE
    reason <- sprintf("no group has >= %d replicates", opt$minReplicates)
  }
}

if (!run_diff) {
  write_status(paste0("SKIPPED differential enrichment: ", reason))
  quit(status = 0)
}

## ---------------------------------------------------------------------------
## Determine baseline group (the group flagged as control == 1)
## ---------------------------------------------------------------------------
if (!"control" %in% colnames(design)) {
  write_status("SKIPPED differential enrichment: no 'control' column to define baseline")
  quit(status = 0)
}
controlFlag <- suppressWarnings(as.integer(design[commonSamples, "control"]))
baselineGroups <- unique(groups[which(controlFlag == 1)])
baselineGroups <- baselineGroups[!is.na(baselineGroups)]
if (length(baselineGroups) == 0) {
  write_status("SKIPPED differential enrichment: no sample flagged as control (control == 1)")
  quit(status = 0)
}
baseline <- baselineGroups[1]
if (length(baselineGroups) > 1) {
  message("Multiple control groups found; using '", baseline, "' as baseline.")
}

## ---------------------------------------------------------------------------
## Build (or load) the barbieQ object and run pairwise differential tests
## ---------------------------------------------------------------------------
sampleMetadata <- S4Vectors::DataFrame(group = factor(groups))
rownames(sampleMetadata) <- commonSamples

# The object is always (re)written here so a subsequent run can reuse it.
rds_out <- out_path("barbieQ.rds")
# Preferred reuse source: an explicit --rdsIn (e.g. staged by Nextflow from
# params.de_rds); otherwise fall back to a barbieQ.rds already in the output
# dir (handy for local re-runs). Building/tagging the object is the dominant
# cost for large barcode sets, so we skip it whenever a usable cache exists.
rds_in <- NULL
if (!is.null(opt$rdsIn) && nzchar(opt$rdsIn) && file.exists(opt$rdsIn)) {
  rds_in <- opt$rdsIn
} else if (file.exists(rds_out)) {
  rds_in <- rds_out
}

## Reuse a previously built object when it is consistent with the current
## input: identical samples (columns) and overlapping barcodes; otherwise
## rebuild. cm is aligned to the cached object so the per-comparison minCount
## filter (which subsets bbq by a logical over cm rows) stays row-consistent.
bbq <- NULL
if (!is.null(rds_in)) {
  verbose_log(paste0("Existing barbieQ object found at ", rds_in,
                     "; attempting to load it and skip the (slow) build."))
  cached <- tryCatch(readRDS(rds_in), error = function(e) {
    verbose_log(paste0("Could not read cached object (", conditionMessage(e),
                       "); rebuilding.")); NULL
  })
  if (!is.null(cached)) {
    sameSamples <- identical(as.character(colnames(cached)),
                             as.character(commonSamples))
    sharedBarcodes <- intersect(rownames(cached), rownames(cm))
    if (!sameSamples) {
      verbose_log("Cached object samples differ from current input; rebuilding.")
    } else if (!length(sharedBarcodes)) {
      verbose_log("Cached object shares no barcodes with current input; rebuilding.")
    } else {
      bbq <- cached[sharedBarcodes, ]
      cm <- cm[sharedBarcodes, , drop = FALSE]
      verbose_log(paste0("Loaded cached barbieQ object (", length(sharedBarcodes),
                         " barcodes); skipping build and tagging."))
    }
  }
}

if (is.null(bbq)) {
  ## Global low-abundance pre-reduction (speed only; never overrides the
  ## per-comparison filter). A barcode whose TOTAL across all samples is below
  ## minCount can never reach minCount within a single contrasted pair (a pair
  ## sum is always <= the global sum), so dropping such barcodes here removes
  ## only barcodes that every per-comparison filter would drop anyway. When the
  ## filter is disabled (minCount 0) we still drop all-zero barcodes, which
  ## contribute nothing and can never be tested. The per-comparison minCount
  ## filter still runs below, unchanged.
  globalTotals <- rowSums(cm)
  prefilterKeep <- if (opt$minCount > 0) globalTotals >= opt$minCount else globalTotals > 0
  nBefore <- nrow(cm)
  cm <- cm[prefilterKeep, , drop = FALSE]
  verbose_log(paste0("Object-construction pre-reduction: kept ", nrow(cm), " of ",
                     nBefore, " barcodes (global count ",
                     if (opt$minCount > 0) paste0(">= ", opt$minCount) else "> 0",
                     "); the per-comparison minCount filter still applies below."))

  verbose_log("Building barbieQ object and tagging top barcodes ...")
  bbq <- barbieQ::createBarbieQ(object = cm, sampleMetadata = sampleMetadata)
  bbq <- tryCatch(barbieQ::tagTopBarcodes(bbq), error = function(e) bbq)
}

# Persist the object so a subsequent re-run can skip the build.
saveRDS(bbq, rds_out)

# Save the tagged top barcodes for the report / downstream use.
topTab <- tryCatch({
  rd <- as.data.frame(SummarizedExperiment::rowData(bbq))
  rd$Barcode <- rownames(rd)
  rd
}, error = function(e) NULL)
if (!is.null(topTab)) {
  if (!is.null(bcIdMap) && "Barcode" %in% colnames(topTab)) {
    topTab$bc_id <- unname(bcIdMap[topTab$Barcode])
    topTab <- topTab[, unique(c("bc_id", "Barcode",
                                setdiff(colnames(topTab), c("bc_id", "Barcode"))))]
  }
  write.table(topTab, out_path("barbieQ_topBarcodes.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# Coefficient terms follow paste0(columnName, level); sanitize like R does.
term <- function(level) paste0("group", make.names(level))

extract_results <- function(res) {
  # testBarcodeSignif() returns the barbieQ object and stores the per-barcode
  # statistics as a NESTED DataFrame column, rowData(res)$testingBarcode, whose
  # columns are meanDiff, Amean, t, P.Value, adj.P.Val, direction, tendencyTo.
  rd <- tryCatch(SummarizedExperiment::rowData(res), error = function(e) NULL)
  if (is.null(rd)) return(NULL)

  # Primary path: pull the nested stats DataFrame directly.
  if ("testingBarcode" %in% colnames(rd)) {
    tab <- tryCatch(as.data.frame(rd$testingBarcode), error = function(e) NULL)
    if (!is.null(tab) && nrow(tab)) {
      tab$Barcode <- rownames(rd)
      return(tab)
    }
  }

  # Fallback: a flattened rowData layout. Match both plain and nested
  # (e.g. 'testingBarcode.P.Value') p-value / stat column names.
  df <- tryCatch(as.data.frame(rd), error = function(e) NULL)
  if (!is.null(df) && nrow(df)) {
    pcols <- c("P.Value", "adj.P.Val", "PValue", "adj.P", "p.value", "padj", "FDR", "logFC")
    hit <- any(pcols %in% colnames(df)) ||
      length(grep("(^|\\.)(P\\.Value|adj\\.P\\.Val|PValue|padj|FDR|logFC)$",
                  colnames(df))) > 0
    if (hit) {
      df$Barcode <- rownames(df)
      return(df)
    }
  }
  NULL
}

# Run a single groupA-vs-groupB comparison and return the extracted table (or
# NULL). meanDiff is oriented as groupA relative to groupB.
run_pair <- function(gA, gB, idx = NA, total = NA) {
  t0 <- Sys.time()
  tag <- if (!is.na(idx) && !is.na(total)) paste0("[", idx, "/", total, "] ") else ""
  cf <- sprintf("(%s) - (%s)", term(gA), term(gB))

  # Low-abundance filter restricted to the two groups being contrasted: keep
  # only barcodes whose summed count within this pair's samples reaches
  # --minCount. Counts from any other condition never enter this decision.
  pairSamples <- commonSamples[groups %in% c(gA, gB)]
  if (opt$minCount > 0) {
    keep <- rowSums(cm[, pairSamples, drop = FALSE]) >= opt$minCount
  } else {
    keep <- rep(TRUE, nrow(cm))
  }
  bbqSub <- bbq[keep, ]

  verbose_log(paste0("  ", tag, "Testing differential proportion: ",
                     gA, " vs ", gB, "  [", cf, "]  (",
                     sum(keep), "/", length(keep), " barcodes pass minCount)"))
  res <- tryCatch(
    barbieQ::testBarcodeSignif(barbieQ = bbqSub, contrastFormula = cf,
                               method = "diffProp", transformation = "asin-sqrt"),
    error = function(e) { verbose_log(paste0("  ", tag, "test failed: ",
                          conditionMessage(e))); NULL }
  )
  elapsed <- format(round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), nsmall = 1)
  if (is.null(res)) {
    verbose_log(paste0("  ", tag, "FAILED after ", elapsed, "s"))
    return(NULL)
  }
  tab <- extract_results(res)
  if (!is.null(tab)) {
    tab$comparison <- paste0(gA, "_vs_", gB)
    verbose_log(paste0("  ", tag, "OK in ", elapsed, "s"))
    return(tab)
  }
  verbose_log(paste0("  ", tag, "OK in ", elapsed, "s but no extractable result table"))
  NULL
}

# Run a list of group pairs (each c(gA, gB)) and return the combined, bc_id-
# annotated results data frame (or NULL when nothing was extractable).
run_pairs <- function(pairs, label = "") {
  total <- length(pairs)
  runner <- function(k) run_pair(pairs[[k]][1], pairs[[k]][2], idx = k, total = total)
  nWorkers <- min(opt$threads, total)
  useParallel <- nWorkers > 1 &&
    requireNamespace("parallel", quietly = TRUE) &&
    .Platform$OS.type == "unix"
  if (useParallel) {
    verbose_log(paste0(label, "Running ", total,
                       " comparison(s) in parallel with ", nWorkers, " worker(s)."))
    lst <- parallel::mclapply(seq_len(total), runner, mc.cores = nWorkers)
  } else {
    verbose_log(paste0(label, "Running ", total, " comparison(s) sequentially."))
    lst <- lapply(seq_len(total), runner)
  }
  lst <- lst[vapply(lst, is.data.frame, logical(1))]
  if (!length(lst)) return(NULL)
  combined <- do.call(rbind, lapply(lst, function(d) {
    d[, unique(c("Barcode", "comparison", setdiff(colnames(d), c("Barcode", "comparison"))))]
  }))
  if (!is.null(bcIdMap)) {
    combined$bc_id <- unname(bcIdMap[combined$Barcode])
    combined <- combined[, unique(c("bc_id", "Barcode", "comparison",
                                    setdiff(colnames(combined),
                                            c("bc_id", "Barcode", "comparison"))))]
  }
  combined
}

allGroups <- levels(factor(groups))

# ---------------------------------------------------------------------------
# Baseline comparisons: every group vs the control baseline. This is what the
# report consumes (barbieQ_diffProp_results.tsv).
# ---------------------------------------------------------------------------
otherGroups <- setdiff(allGroups, baseline)
verbose_log(paste0("Baseline group: '", baseline, "'"))
verbose_log(paste0("Groups to test against baseline: ",
                   paste(otherGroups, collapse = ", ")))
verbose_log(paste0("Baseline comparisons planned: ", length(otherGroups)))

baselinePairs <- lapply(otherGroups, function(g) c(g, baseline))
baselineRes <- run_pairs(baselinePairs, "Baseline: ")

nBaseline <- 0
if (!is.null(baselineRes)) {
  write.table(baselineRes, out_path("barbieQ_diffProp_results.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  nBaseline <- length(unique(baselineRes$comparison))
}

# ---------------------------------------------------------------------------
# Optional all-vs-all: every unique group pair. Written to a SEPARATE file so
# the report keeps showing only the baseline comparisons.
# ---------------------------------------------------------------------------
nAllVsAll <- 0
if (isTRUE(opt$allVsAll)) {
  if (length(allGroups) < 2) {
    verbose_log("--allVsAll requested but fewer than 2 groups; skipping.")
  } else {
    allPairs <- utils::combn(allGroups, 2, simplify = FALSE)
    verbose_log(paste0("All-vs-all comparisons planned: ", length(allPairs)))
    allVsAllRes <- run_pairs(allPairs, "All-vs-all: ")
    if (!is.null(allVsAllRes)) {
      write.table(allVsAllRes, out_path("barbieQ_diffProp_results_allVsAll.tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      nAllVsAll <- length(unique(allVsAllRes$comparison))
    }
  }
}

if (nBaseline > 0 || nAllVsAll > 0) {
  status <- paste0("OK: differential enrichment for ", nBaseline,
                   " comparison(s) vs baseline '", baseline, "'")
  if (isTRUE(opt$allVsAll)) {
    status <- paste0(status, "; all-vs-all wrote ", nAllVsAll, " comparison(s)")
  }
  write_status(status)
  verbose_log("===== barbieQ analysis complete =====")
} else {
  write_status("WARNING: differential testing ran but produced no extractable results")
  verbose_log("===== barbieQ analysis complete (with warnings) =====")
}
