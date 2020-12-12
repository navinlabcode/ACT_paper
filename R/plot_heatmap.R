#' Plot the single-cell heatmap
#'
#' @param df The ordered segmented single-cell data frame
#' @param ploidy_VAL The inferred sample ploidy value
#' @param ploidy_trunc Numeric. Value that will truncate levels higer than it for visualization
#' @param clusters Data frame with clustering information
#' @param genomic_classes Vector with genomic regions classification
#' @param keep_gene Genes that will be kept on the annotation. Must be provided otherwise it will annotate all genes in hg19
#' @param tree_order Ordering of subclones resulting from the me tree. Will be used to order the heatmap
#' @param show_legend Boolean. Display heatmap legend
#'
#' @return A ComplexHeatmap with the single-cells annotated by clustering.
#' @export
#'
#' @examples
plot_heatmap <- function(df,
                         ploidy_VAL,
                         ploidy_trunc,
                         clusters,
                         genomic_classes = NULL,
                         keep_gene,
                         tree_order = NULL,
                         show_legend = FALSE) {

  if (!is.null(tree_order)) {
    clusters <- clusters[order(match(clusters$subclones, tree_order)), ]
  }

  # check to avoid crashing with too many genes annotations
  if (length(keep_gene) > 40) {
    stop("Reduce the list of genes to be annotated. Max = 40")
  }

  ###### chromosome bar top and annotations basic vectors
  # loading data containing vector chr
  chr_lengths <-
    c(
      953,
      1038,
      875,
      838,
      782,
      724,
      646,
      635,
      480,
      561,
      568,
      575,
      431,
      391,
      332,
      310,
      311,
      336,
      220,
      266,
      149,
      139,
      607,
      38
    )
  chr_binary <- rep(c(2, 1), 12)
  chr <-
    data.frame(chr = rep.int(x = chr_binary, times = chr_lengths))

  # getting lengths for chr numbers annotation
  chr_rl_c <- c(1, cumsum(chr_lengths))

  # creating a data frame to calculate rowMeans
  chr_df <-
    data.frame(a = chr_rl_c[1:length(chr_rl_c) - 1], b = chr_rl_c[2:length(chr_rl_c)])
  chr_l_means <- round(rowMeans(chr_df))

  chrom.names <- c(1:22, "X", "Y")

  # creating the vector for chr number annotations
  v <- vector(length = sum(chr_lengths), mode = "character")
  v[chr_l_means] <- chrom.names
  v[is.na(v)] <- ""

  # chr bar with the chr names
  chr_bar <- HeatmapAnnotation(
    chr_text = anno_text(v[1:ncol(df)],
                         gp = gpar(fontsize = 16)),
    df = as.character(chr[1:12167, ]),
    show_legend = FALSE,
    show_annotation_name = FALSE,
    which = "column",
    col = list(df = c("1" = "grey88", "2" = "black"))
  )

  # colors
  color_heat <- structure(pals::ocean.balance(length(0:ploidy_trunc)),
                          names = 0:ploidy_trunc)

  # special ploidy colors if ground state rounds to 2
  if (round(ploidy_VAL) == 2) {
    color_heat <- structure(
      c(
        "#3787BA",
        "#95B8C5",
        "#F0ECEB",
        "#D7A290",
        "#BF583B",
        "#8D1128",
        "#3C0912"
      ),
      names = c("0", "1", "2", "3", "4", "5", "6")
    )
  }


  # removing cells from ann_df to be used as an annotation
  ann_df <- clusters %>%
    as.data.frame() %>%
    dplyr::select(superclones,
                  subclones)

  popseg_round <- ploidy_scale(ploidy_VAL, df)
  popseg_heatmap <- popseg_round
  popseg_heatmap[popseg_round > ploidy_trunc] <- ploidy_trunc

  colors_superclones <- colors_vector$superclones
  colors_subclones <- colors_vector$subclones

  ann <- rowAnnotation(
    df = ann_df,
    col = list(subclones = colors_subclones,
               superclones = colors_superclones),
    show_annotation_name = TRUE,
    annotation_name_gp = gpar(fontsize = 17),
    simple_anno_size = unit(0.9, "cm"),
    annotation_legend_param = list(
      subclones = list(
        labels = gtools::mixedsort(unique(as.character(clusters$subclones))),
        at = gtools::mixedsort(unique(as.character(clusters$subclones)))
      ),
      superclones = list(labels = gtools::mixedsort(unique(clusters$superclones)))
    ),
    show_legend = show_legend
  )

  # genes annotation
  browser()
  bins_bed <- bins_in_cna_pipeline %>%
    dplyr::filter(chr != "chrY") %>%
    mutate(regions = genomic_classes)

  bins_gr <- bins_bed %>%
    makeGRangesFromDataFrame(keep.extra.columns = T, ignore.strand = T)

  #extras for mdamb231
  hg19_genes <- GenomicFeatures::genes(txdb, columns = "SYMBOL")

  hg19_genes <- hg19_genes %>%
    as.data.frame() %>%
    filter(SYMBOL %in% keep_gene) %>%
    select(seqnames, start, end, SYMBOL) %>%
    dplyr::rename(chr = "seqnames",
                  symbol = "SYMBOL") %>%
    mutate(symbol = as.character(symbol))

  bed <- hg19_genes %>%
    makeGRangesFromDataFrame(ignore.strand = T, keep.extra.columns = T)

  olaps <- findOverlaps(bed, bins_gr)

  mk_df <- tibble(
    gene = bed$symbol[queryHits(olaps)],
    pos = subjectHits(olaps),
    region = bins_bed$regions[subjectHits(olaps)]
  ) %>%
    dplyr::distinct(gene, .keep_all = TRUE) %>%
    mutate(color = case_when(
      str_detect(region, "cCNA") ~ "#414451",
      str_detect(region, "sCNA") ~ "gray65",
      str_detect(region, "uCNA") ~ "#FF800E"
    ))

  if (!is.null(keep_gene)) {
    mk_df <- mk_df %>%
      dplyr::filter(gene %in% keep_gene)
  }

  mk <-
    ComplexHeatmap::columnAnnotation(
      df = genomic_classes,
      foo = anno_mark(
        at = mk_df$pos,
        labels = mk_df$gene,
        side = "bottom",
        labels_gp = gpar(fontsize = 17, col = mk_df$color)
      ),
      col = list(df = c(
        "cCNA" = "#414451",
        "sCNA" = "gray65",
        "uCNA" = "#FF800E"
      )),
      show_annotation_name = FALSE,
      show_legend = FALSE
    )


  ht <- Heatmap(
    popseg_heatmap,
    use_raster = TRUE,
    column_title = "genomic coordinates",
    column_title_gp = gpar(fontsize = 24),
    column_title_side = "bottom",
    top_annotation = chr_bar,
    cluster_rows = FALSE,
    border = TRUE,
    row_title = paste(nrow(popseg_heatmap), "single cells"),
    row_title_gp = gpar(fontsize = 24),
    row_title_side = "right",
    cluster_columns = FALSE,
    show_column_names = FALSE,
    show_row_names = FALSE,
    col = color_heat,
    heatmap_legend_param = list(title = "copy number"),
    left_annotation = ann,
    bottom_annotation = mk,
    show_heatmap_legend = show_legend
  )

  ht

}


}
