# Custom version of progeny that supports sparse matrix as input
progeny.CaTCH = function(expr, scale = TRUE, organism = "Human", top = 100, verbose = FALSE,...) {
  
  if (!is.logical(scale)){
    stop("scale should be a logical value")
  }
  
  if (!is.logical(verbose)){
    stop("verbose should be a logical value")
  }
  
  model <- getModel(organism, top=top)
  common_genes <- intersect(rownames(expr), rownames(model))
  
  if (verbose){
    number_genes <- apply(model, 2, function (x) {
      sum(rownames(model)[which (x != 0)] %in% unique(rownames(expr)))
    })
    message("Number of genes used per pathway to compute progeny scores:")
    message(paste(names(number_genes),": ", number_genes, " (", 
                  (number_genes/top)*100,"%)",sep = "","\n"))
  }
  
  result <- t(expr[common_genes,,drop=FALSE]) %*% 
    as.matrix(model[common_genes,,drop=FALSE])
  
  if (scale && nrow(result) > 1) {
    rn <- rownames(result)
    result <- apply(result, 2, scale)
    rownames(result) <- rn
  }
    
  return(result)    
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
    select(cluster, symbol, summary.logFC) %>%
    pivot_wider(names_from = cluster, values_from = summary.logFC) %>%
    data.frame(row.names = "symbol", check.names = FALSE) %>%
    as.matrix()
  mat_stats <- stats_tfa %>%
    select(cluster, symbol, FDR) %>%
    pivot_wider(names_from = cluster, values_from = FDR) %>%
    data.frame(row.names = "symbol", check.names = FALSE) %>%
    as.matrix()
  
  stats_value <- 0.05
  mat <- mat_effect
  mat[mat_stats >= stats_value] <- 0
  max_mat <- max(abs(mat_effect), na.rm = TRUE)
  
  col_pwa <- colorRamp2(breaks = c(-max_mat, 0, max_mat), colors = c("blue", "white", "red"))
  
  p_hm_raw <- Heatmap(matrix = mat,
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
                      left_annotation = row.annotation)
  return (p_hm_raw)
}
