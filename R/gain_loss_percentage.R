#' Calculates the percentage of gains and losses relative to the mean ploidy value
#' of a sample
#'
#' @param consensus The consensus profile from a sample
#' @param ploidy_VAL The FACS derived mean ploidy value
#' @param sample Character. Annotation with sample name added to return df
#'
#' @return A data frame containing the percentage of gains/neutral/losses
#' @export
#'
#' @examples
gain_loss_percentage <- function(consensus,
                                 ploidy_VAL,
                                 sample) {
  ps_cs1_int <- as.data.frame(t(ploidy_scale(ploidy_VAL, consensus)))

  ps_int_cnfreq <- purrr::map(ps_cs1_int, janitor::tabyl)

  ps_int_cnfreq_df <- bind_rows(ps_int_cnfreq, .id = "clone")

  names(ps_int_cnfreq_df)[2] <- "CN"

  gain_loss_df <- ps_int_cnfreq_df %>%
    mutate(class = case_when(
      CN > round(ploidy_VAL) ~ "gain",
      CN < round(ploidy_VAL) ~ "loss",
      TRUE ~ "ground_state"
    )) %>%
    group_by(clone, class) %>%
    tally(percent) %>%
    mutate(sample = sample)

  return(gain_loss_df)

}
