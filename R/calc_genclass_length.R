#' From the genomic classes annotation calculates counts and lengths
#'
#' @param annotations Vector with genomic classes annotations
#'
#' @return A data frame with genomic classes annotation summarized by chromosome arm
#' @export
#'
#' @examples
calc_genclass_length <-  function(annotations,
                                  popseg,
                                  popseg_long) {
  popseg_long_t <- as.data.frame(t(popseg_long))
  seg_index = rep.int(1:nrow(popseg), popseg$n.probes)

  popseg_long_anno <- popseg_long_t %>%
    mutate(
      class = annotations,
      chr = rep.int(paste0(popseg$chrom, popseg$arm), popseg$n.probes),
      start = rep.int(popseg$start.pos, popseg$n.probes),
      end = rep.int(popseg$end.pos, popseg$n.probes),
      seg_length = end - start,
      seg_index = seg_index
    ) %>%
    select(chr, start, end, class, seg_length, seg_index)

  return(popseg_long_anno)
}
