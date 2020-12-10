#' Run HDBSCAN clustering on the resulting embedding from umap
#'
#' @param umap_df Data frame with the resulting embedding from \code{run_umap}.
#' @param k_snn_major nearest neighbors K value. Passed to build an SNN graph that will be used to determine superclones
#' @param k_snn_minor nearest neighbors K value. Passed to \code{dbscan::hdbscan} minPts argument.
#'
#' @return A data frame containing the superclone and subclone assignment for each cell
#' @export
#'
#' @examples
#'
#'
run_clustering <- function(umap_df,
                           k_snn_major = 35,
                           k_snn_minor = 17) {

  library(scran)
  # building a snn graph for superclones
  message("Building SNN graph.")
  g_major <- scran::buildSNNGraph(umap_df[,c(1:2)], k = k_snn_major, transposed = T)

  g_clusters <- igraph::membership(igraph::components(g_major))
  g_clusters <- paste0("s", g_clusters)

  # Clustering
  message("Running hdbscan.")
  subclones <- dbscan::hdbscan(umap_df[,c(1:2)],
                               minPts = k_snn_minor)
  umap_df$subclones <- paste0("c",subclones$cluster)

  # for hdb
  # adding the ones classified as outliers to the closest cluster possible according to euclidean distance
  dist_umap <- dist(umap_df[,c(1:2)]) %>% as.matrix() %>% as.data.frame() %>%
    rownames_to_column("cell2") %>%
    gather(key = "cell1",
           value = "dist",
           -cell2) %>%
    dplyr::filter(cell1 != cell2)

  dist_min <- dist_umap %>%
    right_join(umap_df %>% dplyr::select(cell, subclones), by = c("cell2" = "cell")) %>%
    filter(subclones != "c0") %>%
    group_by(cell1) %>%
    slice_min(dist) %>%
    ungroup()


  for (i in 1:nrow(umap_df)) {

    if(umap_df$subclones[i] == "c0") {
      cellname <- rownames(umap_df)[i]
      closest_cell <- filter(dist_min, cell1 == rownames(umap_df)[i])$cell2
      closest_cell_cluster <- filter(umap_df, cell == closest_cell)$subclones
      umap_df$subclones[i] <- closest_cell_cluster

    }

  }

  cl_df <- tibble::tibble(
    superclones = g_clusters,
    subclones =  umap_df$subclones,
    cells = rownames(umap_df)
  )

  cl_df <- cl_df %>%
    arrange(superclones, subclones)

  # calculating number of cells in every cluster
  freq_df <- janitor::tabyl(umap_df$subclones)
  names(freq_df)[1] <- "cluster"
  print(freq_df)

  classification <- cl_df

  message("Done.")

  return(classification)

}
