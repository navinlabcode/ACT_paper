#' Run Umap dimension reduction
#'
#' @param popseg_long The resulting data from segmentation
#' @param ploidy_VAL Optional ploidy value to scale the segment ratios values by its ploidy
#' @param umap_dist Distance metric used by umap. Defaults to manhattan
#' @param umap_min_dist UMAP min dist parameter
#' @param umap_spread UMAP spread parameter
#' @param umap_n_neighbors UMAP n_neighbors parameter
#' @param mc.cores Number of threads to be used
#' @param seed Seed, defaults to 55
#' @param round Optional in case data is ploidy scaled values are rounded to the nearest integer
#'
#' @return A data frame with the resulting embedding from umap and a column containing the cell names.
#' @export
#'
#' @examples
#'
run_umap <- function(popseg_long,
                     ploidy_VAL = NULL,
                     umap_dist = "manhattan",
                     umap_min_dist = 0,
                     umap_spread = 1,
                     umap_n_neighbors = 40,
                     mc.cores = 10,
                     seed = 55,
                     round = FALSE) {

  if (round == TRUE) {
    popseg_long <- ploidy_scale(ploidy_VAL = ploidy_VAL, popseg_long, round = round)
  }

  message("Constructing UMAP embedding.")
  set.seed(seed)
  dat_umap <- uwot::umap(popseg_long,
                         metric = umap_dist,
                         min_dist = umap_min_dist,
                         n_neighbors = umap_n_neighbors,
                         spread = umap_spread,
                         n_components = 2,
                         n_thread = mc.cores)

  umap_df <- as.data.frame(dat_umap)

  rownames(umap_df) <- rownames(popseg_long)
  umap_df$cell <- rownames(umap_df)

  return(umap_df)

}
