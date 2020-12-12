#' Plot the consensus heatmap with the inferred MRCA profile and gene annotation
#'
#' @param df The data frame with the inferred consensus profiles
#' @param clusters data.frame. Resulting from clustering
#' @param ploidy_VAL Numeric. The inferred ploidy value
#' @param ploidy_trunc Numeric. Value that will truncate levels higer than it for visualization
#' @param keep_gene Genes that will be kept on the annotation. Must be provided otherwise it will annotate all genes in hg19
#' @param tree_order Ordering of subclones resulting from the me tree. Will be used to order the heatmap
#' @param plot_title Title of the plot
#' @param genomic_classes Character. Vector with genomic regions classification
#'
#' @return
#' @export
#'
#' @examples
plot_consensus_heatmap <- function(df,
                                    clusters,
                                    ploidy_VAL,
                                    ploidy_trunc,
                                    keep_gene = NULL,
                                    tree_order = NULL,
                                    plot_title = NULL,
                                    genomic_classes = NULL) {

  # check to avoid crashing with too many genes annotations
  if (length(keep_gene) > 40) {
    stop("Reduce the list of genes to be annotated. Max = 40")
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Thu Nov  5 15:36:32 2020
  # adding an ancestral profile
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Thu Nov  5 15:36:41 2020

  anno_split <- clusters %>%
    distinct(subclones, .keep_all = T)

  anno_split_o <- anno_split[match(rownames(df), anno_split$subclones),]

  cs_list <- split(df, anno_split_o$superclones)

  cs_list_median <- lapply(cs_list, function(x) { apply(x, 2, median) })

  cs_median <- bind_rows(cs_list_median)

  rownames(cs_median) <- names(cs_list_median)

  cs_median_int <- ploidy_scale(ploidy_VAL, cs_median)

  #obtain number closest to the ground state for each
  anc_profile <- apply(cs_median_int, 2, function(x) x[which.min(abs(x-round(ploidy_VAL)))] )

  anc_profile_df <- data.frame(anc_profile = anc_profile)
  anc_profile_df <- as.data.frame(t(anc_profile_df))

  ann_df <- clusters %>%
    dplyr::mutate_all(as.character) %>%
    dplyr::select(-cells) %>%
    distinct(.keep_all = T) %>%
    as.data.frame()

  ann_df <- ann_df[match(rownames(df), ann_df$subclones),]

  # colors
  color_heat <- structure(
    pals::ocean.balance(length(0:ploidy_trunc)),
    names = 0:ploidy_trunc
  )

  # special ploidy for ploidy that rounds to 2
  if (round(ploidy_VAL) == 2) {

    color_heat <- structure(c("#3787BA",
                              "#95B8C5",
                              "#F0ECEB",
                              "#D7A290",
                              "#BF583B",
                              "#8D1128",
                              "#3C0912"),
                            names = c("0","1","2","3","4","5", "6"))
  }

  # chromosome annotation top bar
  # loading data containing vector chr
  chr_lengths <- c(953, 1038,  875,  838,  782,  724,  646,  635,  480,  561,  568,  575,  431,  391,  332,  310,  311,  336,  220, 266,  149 , 139,  607,   38)
  chr_binary <- rep(c(2,1), 12)
  chr <- data.frame(chr = rep.int(x = chr_binary, times = chr_lengths))

  # getting lengths for chr numbers annotation
  chr_rl_c <- c(1, cumsum(chr_lengths))

  # creating a data frame to calculate rowMeans
  chr_df <-  data.frame(a = chr_rl_c[1:length(chr_rl_c)-1],b= chr_rl_c[2:length(chr_rl_c)])
  chr_l_means <- round(rowMeans(chr_df))

  chrom.names <- c(1:22,"X", "Y")

  # creating the vector for chr number annotations
  v <- vector(length = sum(chr_lengths), mode = "character")
  v[chr_l_means] <- chrom.names
  v[is.na(v)] <- ""

  # chr bar with the chr names
  chr_bar <- HeatmapAnnotation(chr_text = anno_text(v[1:ncol(df)],
                                                    gp = gpar(fontsize = 16)),
                               df = as.character(chr[1:12167,]),
                               show_legend = FALSE,
                               show_annotation_name = FALSE,
                               which = "column",
                               col = list(df = c("1" = "grey88", "2" = "black")))

  # genes annotation

  bins_bed <- bins_in_cna_pipeline %>%
    dplyr::filter(chr != "chrY") %>%
    mutate(regions = genomic_classes)

  bins_gr <- bins_bed %>%
    makeGRangesFromDataFrame(keep.extra.columns = T, ignore.strand = T)

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

  mk_df <- tibble(gene = bed$symbol[queryHits(olaps)],
                  pos = subjectHits(olaps),
                  region = bins_bed$regions[subjectHits(olaps)]) %>%
    dplyr::filter(gene %in% keep_gene) %>%
    dplyr::distinct(gene, .keep_all = TRUE) %>%
    mutate(color = case_when(
      str_detect(region, "cCNA") ~ "#414451",
      str_detect(region, "sCNA") ~ "gray65",
      str_detect(region, "uCNA") ~ "#FF800E"
    ))

  mk <-
    ComplexHeatmap::columnAnnotation(
      df = genomic_classes,
      foo = anno_mark(
        at = mk_df$pos,
        labels = mk_df$gene,
        side = "bottom",
        labels_gp = gpar(fontsize = 17, col = mk_df$color)
      ),
      col = list(
        df = c(
          "cCNA" = "#414451",
          "sCNA" = "gray65",
          "uCNA" = "#FF800E"
        )
      ),
      name = "CNA class",
      show_annotation_name = FALSE,
      show_legend = FALSE
    )

  colors_subclones <- colors_vector$subclones

  # ordering by tree
  if (!is.null(tree_order)) {
    ann_df <- ann_df[match(tree_order, ann_df$subclones),]
  }

  ann_df[nrow(ann_df) + 1,] <- c(NA, NA)
  ann_df <- ann_df[match(c("anc_profile", tree_order), ann_df$subclones),]

  ann <- rowAnnotation(df = ann_df,
                       col = list(subclones = c(colors_subclones, "anc_profile" = "black") ,
                                  superclones = c(colors_vector$superclones, "anc_profile" = "black")),
                       show_annotation_name = FALSE,
                       simple_anno_size = unit(0.9, "cm"),
                       show_legend = F)


  popseg_round <- ploidy_scale(ploidy_VAL, df)
  popseg_round[nrow(popseg_round) + 1,] <- anc_profile
  rownames(popseg_round)[nrow(popseg_round)] <- "anc_profile"
  popseg_heatmap <- popseg_round

  popseg_heatmap[popseg_heatmap > ploidy_trunc] <- ploidy_trunc
  rownames(popseg_heatmap) <- rownames(popseg_round)

  if (!is.null(tree_order)) {
    popseg_heatmap <- popseg_heatmap[c("anc_profile", tree_order),]
  }

  ht <- Heatmap(popseg_heatmap,
                use_raster = TRUE,
                row_split = c("anc_profile", rep(1, (nrow(popseg_heatmap)-1))),
                row_title = NULL,
                top_annotation = chr_bar,
                bottom_annotation = mk,
                column_title = "genomic coordinates",
                column_title_gp = gpar(fontsize = 24),
                column_title_side = "bottom",
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                show_column_names = FALSE,
                show_row_names = F,
                border = TRUE,
                col = color_heat,
                show_heatmap_legend = FALSE)

    ComplexHeatmap::draw(ann+ht,
                         row_title = "Consensus",
                         row_title_gp = gpar(fontsize = 24),
                         padding = unit(c(1, .1, 1, .1), "cm"))

}
