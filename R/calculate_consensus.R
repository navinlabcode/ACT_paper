#' Calculate the median consensus profile of every subclone
#'
#' @param df The ordered population segmentation profile
#' @param clusters The cluster assignment of each single-cell
#' @param mc.cores Number of threads.
#'
#' @return
#' @export
#'
#' @examples
calculate_consensus <- function(df,
                                clusters,
                                mc.cores = 30) {

  ## reading list with clusters
  long_list <- split(df, clusters)

  consensus_list <- parallel::mclapply(long_list, function(x) { apply(x, 2, median) }, mc.cores = mc.cores)

  cs_df <- as.data.frame(do.call(rbind, consensus_list))

  return(cs_df)

}
