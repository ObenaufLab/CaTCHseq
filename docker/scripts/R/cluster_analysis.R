library(tidyverse)
library(scater)
library(scran)
library(SingleCellExperiment)
library(dorothea)
library(progeny)
library(ComplexHeatmap)
library(circlize)
library(DESeq2)

source("lib.R")

#load("/data/gcbds/users/nowoshil/Projects/CaTCH2.0/KRASi_ON_vs_OFF/Exp1010/_NF_run_/outputs/sce/sce.rda")
#rm(sce.unfiltered)

opt <- list(min_size = 20,
            plots_out = "~/plots/",
            bin_out = "~/",
            subsample = 35000)


#### Data preparation ####
##### Load the DoRothEA regulons #####
dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea")) %>%
                          filter(confidence %in% c("A", "B", "C"))

##### Reference treatment is the first treatment #####
ref.treatment <- levels(sce$Treatment)[1]


#### Use DeSeq2 to identify over- and underrepresented barcodes ####
metadata <- colData(sce) %>% 
            as_tibble() %>% 
            select(Sample, Treatment) %>% 
            distinct()

bc.counts <- colData(sce) %>% 
             as_tibble() %>%
             filter(CaTCH.Status == "Singlet") %>%
             select(CaTCH.BCs, Sample, BC_ID) %>%
             group_by(CaTCH.BCs, Sample) %>%
             mutate(n = n()) %>%
             ungroup() %>%
             distinct() %>%
             pivot_wider(names_from = Sample, values_from = n, values_fill = 0) %>%
             filter(rowSums(across(starts_with(ref.treatment))) > 0)

idx <- match(metadata$Sample, setdiff(colnames(bc.counts), c("CaTCH.BCs", "BC_ID")))
metadata <- metadata[idx, ]

dds <- DESeqDataSetFromMatrix(countData = bc.counts %>%
                                          select(-CaTCH.BCs) %>%
                                          column_to_rownames(var = "BC_ID"),
                              colData = metadata,
                              design= ~ Treatment)
dds <- DESeq(dds)

#### PROGENy and DoRothEA ####
for (t in setdiff(levels(sce$Treatment), ref.treatment)) {
  print(sprintf("   Processing the treatment '%s'...", t))
  
  ##### Extract the DE results for the given treatment
  print("      ... identifying differentially represented CaTCH barcodes...")
  r <- results(dds, 
               alpha = 0.05,
               contrast = c("Treatment", t, ref.treatment)) %>%
       as_tibble(rownames = NA) %>% 
       rownames_to_column("BC_ID") %>%
       filter(!is.na(padj), padj <= 0.05) %>%
       mutate(Type = if_else(log2FoldChange > 0, "Enriched", "Depleted")) %>%
       arrange(desc(Type), desc(abs(log2FoldChange)))
  
  print("      ... filtering the data...")
  mask <- (sce$CaTCH.Status == "Singlet") & (sce$Treatment %in% c(ref.treatment, t))
  sce.tmp <- sce[, mask]
  
  # Subsample the data if there are more than `opt$subsample` cells
  if ((opt$subsample > 0) && (ncol(sce.tmp) > opt$subsample)) {
    print(sprintf("      ... sampling the data (%d => %d)...", ncol(sce.tmp), opt$subsample))
    idx <- 1:ncol(sce.tmp)
    # The cells are naturally grouped by the sample. Therefore, shuffle the
    # cells before subsampling
    idx.shuffle <- sample(idx, length(idx))
    idx <- sample(idx.shuffle, opt$subsample)
    sce.tmp <- sce.tmp[, idx]
  }
  
  # Remove extremely small clusters (less than `opt$min_size`), since
  # no reasonable statistics can be calculated.
  mask <- colData(sce.tmp) %>%
          as_tibble() %>%
          select(Cluster) %>%
          group_by(Cluster) %>%
          mutate(mask = n() >= opt$min_size) %>%
          ungroup() %>%
          select(mask) %>%
          pull()
  sce.tmp <- sce.tmp[, mask]
  sce.tmp$Cluster <- factor(sce.tmp$Cluster, levels = sort(unique(sce.tmp$Cluster)))
  logcounts <- assay(sce.tmp, "logcounts")
  
  print("      ... assigning barcode types...")
  colData(sce.tmp)["BarcodeType"] <- colData(sce.tmp) %>%
                                     as_tibble() %>%
                                     left_join(y = (r %>% 
                                                    select(BC_ID, Type)), 
                                               by = "BC_ID") %>%
                                     mutate(Type = if_else(is.na(BC_ID), "Invalid barcode", Type),
                                            Type = if_else(is.na(Type), "Not significant", Type)) %>%
                                     select(Type)
  
  ##### Prepare the annotation blocks for both PROGENy and DoRothEA #####
  ###### Heatmaps with the clusters ######
  # Annotation line: fraction of cells in each cluster that are either the 
  # reference or the current treatment
  print("      ... preparing the heatmap annotations...")
  ha.data.frac <- colData(sce.tmp) %>%
                  as_tibble() %>%
                  select(Treatment, Cluster) %>%
                  group_by(Cluster) %>%
                  mutate(ClusterSize = n()) %>%
                  group_by(Treatment, Cluster) %>%
                  mutate(n = n()) %>%
                  ungroup() %>%
                  distinct() %>%
                  mutate(Frac = n / ClusterSize * 100) %>%
                  select(-ClusterSize, -n) %>%
                  pivot_wider(names_from = Treatment, values_from = Frac, values_fill = 0) %>%
                  arrange(Cluster) %>%
                  relocate(starts_with(ref.treatment)) %>%
                  column_to_rownames(var = "Cluster") %>%
                  as.matrix()
  
  # Annotation line: cluster sizes irrespective of the treatment
  ha.data.size <- colData(sce.tmp) %>%
                  as_tibble() %>%
                  select(Cluster) %>%
                  group_by(Cluster) %>%
                  mutate(n = n()) %>%
                  ungroup() %>%
                  distinct() %>%
                  arrange(Cluster) %>%
                  column_to_rownames(var = "Cluster") %>%
                  as.matrix()
  max.y <- max(ha.data.size)
  breaks.y <- c(round(max.y / 4), round(max.y / 2), round(max.y / 4 * 3), max.y)
  lgd.list <- list(Legend(labels = c(ref.treatment, t), 
                          title = "Treatment",
                          type = "points",
                          pch = 20,
                          legend_gp = gpar(col = c("forestgreen", "orange"))))
  ha <- HeatmapAnnotation(Size = anno_barplot(ha.data.size, 
                                              border = FALSE, 
                                              axis_param = list(at = breaks.y, labels = sprintf("%d", breaks.y)),
                                              height = unit(3, "cm")),
                          Treatment = anno_barplot(ha.data.frac, border = FALSE, gp = gpar(fill = c("forestgreen", "orange"))),
                          show_legend = c(FALSE, TRUE))
  
  ###### Heatmaps with individual cells ######
  annotations = data.frame(Treatment = sce.tmp$Treatment)
  colors.treatment <- c()
  colors.treatment[[ref.treatment]] <- "forestgreen"
  colors.treatment[[t]] <- "orange"

  annot.colors <- list(Treatment = unlist(colors.treatment))
  
  # Add a separate line for each cluster for better visualization
  for (cl in levels(sce.tmp$Cluster)) {
    annotations[paste0("Cluster ", cl)] <- sce.tmp$Cluster == cl
    annot.colors[[paste0("Cluster ", cl)]] <- c("TRUE" = "black", "FALSE" = "white")
  }

  annotations$Type <- sce.tmp$BarcodeType
  annotations$VPR <- logcounts["dCas9-VPR", ]
  annotations$CellCycle <- sce.tmp$CellStage

  colors.type <- c("Depleted" = "blue", 
                   "Enriched" = "firebrick2", 
                    "Invalid barcode" = "yellow", 
                    "Not significant" = "grey75")
  colors.cellcycle <- c("G0" = "gray75",
                        "G1S" = "#FDB916",
                        "G2M" = "#11B09C",
                        "M" = "#D22248",
                        "MG1" = "#F17724",
                        "S" = "#BECE2D")

  annot.colors$Type <- colors.type
  annot.colors$CellCycle <- colors.cellcycle
  annot.colors$VPR <- colorRamp2(c(min(annotations$VPR), max(annotations$VPR)), c("blue", "red"))

  # Add a separate line for each enriched or depleted CaTCH barcode for better visualization
  for (i in 1:nrow(r)) {
    bc <- r$BC_ID[i]
    mask <- sce.tmp$BC_ID == bc
    mask[mask] <- r$Type[i]
    annotations[bc] <- mask
    annot.colors[[bc]] <- c("Enriched" = "red", "Depleted" = "blue", "FALSE" = "white")
  }
  
  ##### PROGENy #####
  print("      ... running PROGENy...")
  pw_activities <- progeny.CaTCH(logcounts, 
                                 scale = TRUE, 
                                 organism = "Human", 
                                 verbose = TRUE) %>%
                   t()
  
  # source("_analyses.R")
  
  
  ###### Draw the PROGENy heatmap ######
  print("         ... generating the cluster heatmap...")
  tmp <- pw_activities %>%
         as_tibble(rownames = NA) %>%
         rownames_to_column(var = "Pathway") %>%
         pivot_longer(cols = c(-Pathway), names_to = "CellID", values_to = "Activity") %>%
         left_join(y = colData(sce.tmp) %>%
                       as_tibble(rownames = NA) %>%
                       rownames_to_column(var = ".CellID") %>%
                       select(.CellID, Cluster),
                   by = c("CellID" = ".CellID")) %>%
         group_by(Pathway, Cluster) %>%
         mutate(MeanActivity = mean(Activity)) %>%
         ungroup() %>%
         select(-CellID, -Activity) %>%
         distinct() %>%
         arrange(Pathway, Cluster) %>%
         pivot_wider(names_from = Cluster, values_from = MeanActivity) %>%
         column_to_rownames(var = "Pathway") %>%
         as.matrix()
  
  hm.progeny <- Heatmap(tmp,
                        name = "PROGENy score",
                        column_title = sprintf("PROGENy pathway analysis for %s vs %s", t, ref.treatment),
                        cluster_columns = FALSE,
                        cluster_rows = FALSE,
                        top_annotation = ha,
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(sprintf("%.1f", tmp[i, j]), x, y, gp = gpar(fontsize = 8, col = "black"))
                        })
  
  pdf(file = paste0(opt$plots_out, "PROGENy_", t, ".pdf"), width = 14, height = 5)
  draw(hm.progeny, annotation_legend_list = lgd.list)
  dev.off()
  
  
  ###### Draw the heatmap of all cells using the PROGENy scores ######
  print("         ... generating the single cell heatmap...")
  hm.progeny.sc <- Heatmap(matrix = t(pw_activities),
                           name = "PROGENy score",
                           cluster_rows = TRUE,
                           cluster_columns = TRUE,
                           show_row_names = FALSE,
                           show_column_names = TRUE,
                           show_row_dend = FALSE,
                           show_column_dend = TRUE,
                           right_annotation = rowAnnotation(df = annotations,
                                                            col = annot.colors,
                                                            show_legend = str_starts(string = names(annotations), pattern = "^Cluster|BC_", negate = TRUE),
                                                            annotation_height = unit(5, "mm")),
                           use_raster = TRUE)
  w <- 0.2 * (length(annotations) + (2 * nrow(pw_activities)))
  pdf(file = paste0(opt$plots_out, "PROGENy_", t, "_cells.pdf"), width = ceiling(w), height = 15)
  draw(hm.progeny.sc)
  dev.off()

  # Try to free up as much RAM as possible
  gc()

  
  ##### DoRothEA #####
  print("      ... running DoRothEA...")
  print("         ... generating the input data...")
  dorothea.input <- logcounts %>%
                    t() %>%
                    scale() %>%
                    t()
  
  print("         ... running VIPER...")
  dorothea.output <- viper::viper(eset = dorothea.input, 
                                  regulon = df2regulon(dorothea_regulon_human),
                                  method = "none", 
                                  minsize = 10, 
                                  eset.filter = FALSE, 
                                  cores = 1, 
                                  verbose = TRUE)

  #sce.dorothea <- run_viper(sce.tmp,
  #                          dorothea_regulon_human,
  #                          options = list(method = "scale", minsize = 10, eset.filter = FALSE, cores = 1, verbose = FALSE))
  sce.dorothea <- SingleCellExperiment(assays = list(tf_activities = dorothea.output),
                                       colData = colData(sce.tmp))

  #sce.dorothea <- altExp(sce.dorothea, "dorothea")
  #colData(sce.dorothea) <- colData(sce.tmp)

  print("         ... saving the data...")
  save(dorothea.input, sce.dorothea, file = paste0(opt$bin_out, "DoRothEA_", t, "_data.RData"))
  next

  
  ###### Draw DoRothEA heatmap ######
  tf.list <- rownames(dorothea.output) %>%
             str_extract(pattern = "^SMAD[0-9]+.*$|^TEAD[0-9]+.*$|^STAT[0-9]+.*$|^CDH.*$") %>%
             discard(is.na)
  
  tmp <- dorothea.output %>%
         as_tibble(rownames = NA) %>%
         rownames_to_column(var = "TF") %>%
         filter(TF %in% tf.list) %>%
         pivot_longer(cols = c(-TF), names_to = "CellID", values_to = "Activity") %>%
         left_join(y = colData(sce) %>%
                       as_tibble(rownames = NA) %>%
                       rownames_to_column(var = ".CellID") %>%
                       select(.CellID, Cluster),
                   by = c("CellID" = ".CellID")) %>%
         group_by(TF, Cluster) %>%
         mutate(MeanActivity = mean(Activity)) %>%
         ungroup() %>%
         select(-CellID, -Activity) %>%
         distinct() %>%
         arrange(TF, Cluster) %>%
         pivot_wider(names_from = Cluster, values_from = MeanActivity) %>%
         column_to_rownames(var = "TF") %>%
         as.matrix()
  
  # Min and max values
  mm <- ceiling(max(abs(c(max(tmp), min(tmp)))))
  col_fun <- colorRamp2(c(-mm, 0, mm), c("blue", "white", "red"))
  
  hm.dorothea <- Heatmap(tmp[tf.list,],
                         name = "DoRothEA score",
                         column_title = sprintf("DoRothEA TF analysis for %s vs %s", t, ref.treatment),
                         cluster_columns = FALSE,
                         cluster_rows = FALSE,
                         top_annotation = ha,
                         cell_fun = function(j, i, x, y, width, height, fill) {
                           grid.text(sprintf("%.1f", tmp[i, j]), x, y, gp = gpar(fontsize = 8, col = "black"))
                         },
                         col = col_fun)
  
  hm.dorothea.sc <- Heatmap(matrix = t(dorothea.output),
                            name = "DoRothEA score",
                            cluster_rows = TRUE,
                            cluster_columns = FALSE,
                            show_row_names = FALSE,
                            show_column_names = TRUE,
                            show_row_dend = FALSE,
                            show_column_dend = FALSE,
                            column_title_gp = gpar(fontsize = 8),
                            right_annotation = rowAnnotation(df = annotations,
                                                             col = annot.colors,
                                                             show_legend = str_starts(string = names(annotations), pattern = "^Cluster|BC_", negate = TRUE),
                                                             annotation_height = unit(3, "mm"),
                                                             annotation_name_gp = gpar(fontsize = 8)),
                            use_raster = TRUE)
  
  pdf(file = "~/fff.pdf", width = 50, height = 10)
  draw(hm.dorothea.sc)
  dev.off()
  
  p_dorothea.cluster <- generateHeatmaps(markers.cluster,
                                         max.FDR = 0.05,
                                         min.FC = 2,
                                         col.title = "Cluster",
                                         cluster.cols = TRUE,
                                         cluster.rows = TRUE,
                                         col.annotation = ha)
  pdf(file = paste0(opt$plots_out, "Dorothea_", t, ".pdf"), width = 20, height = 35)
  draw(p_dorothea.cluster, annotation_legend_list = lgd.list)
  dev.off()
  
  ###### Draw the heatmap of all cells using the DoRothEA scores ######
  
  pdf(file = paste0(opt$plots_out, "Dorothea_", t, "_cells.pdf"), width = 20, height = 45)
  draw(hm.dorothea.sc)
  dev.off()
  
  ##### Free up resources #####
  rm(sce.tmp, logcounts, sce.dorothea, dorothea.output)
  gc()
}

