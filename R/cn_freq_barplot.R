#' Copy number gene barplot
#'
#' @param act_df Data frame with the segment log ratios for ACT.
#' @param tenx_df Data frame with the segment log ratios for 10XCNA data.
#' @param keep_genes Genes that will be shown on the plot.
#' @param ploidy_VAL Tumor ploidy value derived from FACS.
#'
#' @return A ggplot object with the frequency for each copy number state from the dataset.
#' @export
#'
#' @examples
cn_freq_barplot <- function(act_df,
                            tenx_df,
                            keep_genes = NULL,
                            ploidy_VAL) {


  # ~~~~~~~~~~~~~~~~~~
  # grep bins
  # ~~~~~~~~~~~~~~~~~~

  grep_bins <- function(df) {
    # obtains bins positions from popseg

    bins_bed <- bins_in_cna_pipeline

    bins_gr <- bins_bed %>%
      makeGRangesFromDataFrame(keep.extra.columns = T, ignore.strand = T)

    bed <- genes(Homo.sapiens, columns = "SYMBOL")

    olaps <- findOverlaps(bed, bins_gr)

    bed_df <- as.data.frame(bed)

    mk_df <- tibble(chr = unlist(unname(bed_df$seqnames[queryHits(olaps)])),
                    gene = unlist(unname(bed_df$SYMBOL[queryHits(olaps)])),
                    pos = subjectHits(olaps)) %>%
      distinct(gene, .keep_all = T)


    mk_df <- mk_df %>%
      filter(gene %in% keep_genes)


  }

  # ACT
  int <- ploidy_scale(ploidy_VAL = ploidy_VAL, df = act_df)

  act_mk <- grep_bins(act_df)

  int_markers <- int[, act_mk$pos]

  names(int_markers) <- act_mk$gene

  int_long_act <- int_markers %>%
    gather(key = "gene",
           value = "CN") %>%
    group_by(gene) %>%
    count(CN) %>%
    mutate(freq = n/sum(n),
           gene = as.factor(gene)) %>%
    mutate(tech = "ACT",
           tech_gene = paste(gene, "(ACT)"))

  # 10X

  tenx_int <- ploidy_scale(ploidy_VAL = ploidy_VAL, df = tenx_df)

  tenx_mk <- grep_bins(tenx_df)

  if (!identical(act_mk$gene, tenx_mk$gene)) {
    tenx_mk <- tenx_mk %>%
      filter(gene %in% act_mk$gene)
  }

  tenx_int_markers <- tenx_int[, tenx_mk$pos]

  names(tenx_int_markers) <- tenx_mk$gene

  int_long_tenx <- tenx_int_markers %>%
    gather(key = "gene",
           value = "CN") %>%
    group_by(gene) %>%
    count(CN) %>%
    mutate(freq = n/sum(n),
           gene = as.factor(gene),
           tech = "10X CNA",
           tech_gene = paste(gene, "(10X CNA)"))


  #binding
  # trick to intercalate 2 vectors https://stackoverflow.com/questions/25961897/how-to-merge-2-vectors-alternating-indexes

  df <- bind_rows(int_long_act,
                  int_long_tenx) %>%
    mutate(tech = as.factor(tech)) %>%
    mutate(tech = fct_relevel(tech,
                              c(rbind(unique(int_long_act$tech), unique(int_long_tenx$tech)))))


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # stacked bar plot CN color
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  if (round(ploidy_VAL) == 3) {
    color_heat <- structure(
      c("#212A62",
        "#4B93BA",
        "#166FBB",
        "white",
        "#DFBCB2",
        "#D59E8C",
        "#CB7F65",
        "#C16043",
        "#B53E29",
        "#A01D25",
        "#810E28",
        "#5E0D20",
        "#3C0912",
        "#3C0912",
        "#3C0912",
        "#3C0912",
        "#3C0912",
        "#3C0912",
        "#3C0912"),
      names = 0:18
    )
  }

  if (!is.null(keep_genes)) {

    df <- df %>%
      filter(gene %in% keep_genes)

  }


  p_freq <-  df %>%
    mutate(gene = as.factor(gene)) %>%
    mutate(gene = fct_relevel(gene, keep_genes)) %>%
    ggplot() +
    geom_col(aes(x = tech,
                 y = freq,
                 fill = fct_relevel(as.character(CN), rev(as.character(sort(unique(df$CN)))))),
             position = "fill",
             color = "black",
             width = .8) +
    facet_wrap(vars(gene), nrow = 2) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    theme(strip.background = element_blank(),
          legend.text = element_text(size = 8),
          axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1)) +
    ylab("frequency") +
    xlab("") +
    scale_fill_manual(values = color_heat) +
    labs(fill = "copy number")

  p_freq

}




