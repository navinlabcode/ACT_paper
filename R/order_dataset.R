#' Orderes the segmentation dataset for downstream usage
#'
#' @param popseg_long The resulting data from segmentation
#' @param clustering Resulting data frame from
#' @param mc.cores Number of threads to be used.
#'
#' @return A list with 2 objects. 1. The ordered segmentation dataset and the ordered clustering dataset
#' @export
#'
#' @examples
order_dataset <- function(popseg_long,
                          clustering,
                          mc.cores = 10) {

  or <- popseg_long[dplyr::pull(clustering, cells),]

  clustering <- clustering %>%
    mutate(subclones = factor(subclones,
                              levels = unique(subclones)))

  # performing hclust within the clusters
  list_m <- split(or, clustering$subclones)

  list_hc <- lapply(list_m, function(m) {
    fastcluster::hclust(amap::Dist(
      m, method = "manhattan", nbproc = mc.cores
    ), method = "ward.D2")
  })

  list_o <- lapply(seq_along(list_hc), function(i) {
    list_m[[i]] <- list_m[[i]][list_hc[[i]]$order, ]
  })

  names(list_o) <- names(list_m)

  popseg_hc <- as.data.frame(do.call(rbind, unname(list_o)))

  popseg_long_o <- or[rownames(popseg_hc),]

  #ordering clusters
  clustering <- clustering[match(rownames(popseg_hc), clustering$cells),]

  ordered_list <- list(dataset_ordered = popseg_long_o,
                       clustering_ordered = clustering)

  return(ordered_list)

}
