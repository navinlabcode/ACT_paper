#' For every bin, returns the inferred genomic classification class (cCNA, sCNA, uCNA)
#'
#' @param consensus_df The data frame with the consensus profile for each superclone/subclone
#' @param ploidy_VAL The inferred ploidy value for the tumor
#'
#' @return A vector where every bin is classified in cCNA, sCNA or uCNA
#' @export
#'
#' @examples
#'
consensus_genomic_classes <- function(consensus_df,
                                        ploidy_VAL) {

  consensus_int <- ploidy_scale(ploidy_VAL, consensus_df)

  percent_clonal <- 1
  percent_extant <- 1 / nrow(consensus_int)

  # for every bin
  ps_percents_list <- future_apply(consensus_int, 2, function(x) {
    perc <- janitor::tabyl(x) %>%
      pull(percent)
  })

  bin_classes <- future_lapply(ps_percents_list, function(x) {
    if (any(x == percent_extant)) {
      return("uCNA")
    } else if (any(x == percent_clonal)) {
      return("cCNA")
    } else
      return("sCNA")
  })

  bin_classes <- unlist(unname(bin_classes))

  return(bin_classes)
}
