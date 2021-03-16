#' Reproduces the moving average plots
#'
#' @param chromosome Chromosome to be plotted
#' @param genes list of genes to be plotted as a boxplot.
#'
#' @return a ggplot object with the three plots using patchwork
#' @export
#'
#' @examples
#'
plot_moving_average <- function(chromosome,
                                genes) {
  # Methods:
  # DNA copy number profiles from the expanded clusters are shown
  # by taking the mode of the ith segment from their profiles
  # according to the co-clustering identities.

  # mode function thanks to https://rpubs.com/Mentors_Ubiqum/using_mode
  mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }

  active_chr <- unique(blk_long_gj_jit$chr)[i]

  blk_long_gj2 <- blk_long_gj_jit %>%
    dplyr::filter(chr == chromosome)

  cnt_long_gene_trip_chr <- cnt_long_gene_trip %>%
    filter(gene %in% genes) %>%
    mutate(gene = fct_relevel(gene, genes))

  message("Calculating moving averages")

  genes_common <-
    blk_long_gj2$gene_id[blk_long_gj2$gene_id %in% rownames(cnt_avg)]
  chr <- as.data.frame(cnt_avg[unique(genes_common),])
  mp <- map_df(chr, function(x) {
    evobiR::SlidingWindow(mean, x, 100, 1)
  })

  mp_l <-
    mp %>% mutate(pos = 1:nrow(mp)) %>% gather(key = "sample", value = "window_exp", -pos)

  mp_l <- inner_join(mp_l, cl_info) %>%
    group_by(subclones, pos) %>%
    summarize(mean_cluster = mean(window_exp))

  p0 <- cnt_long_gene_trip_chr %>%
    filter(gene %in% genes) %>%
    ggplot() +
    geom_boxplot(aes(x = subclones,
                     y = z_score,
                     fill = subclones)) +
    facet_wrap(vars(gene), ncol = 1)  +
    scale_fill_manual(values = colors_vector$subclones) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))  +
    ylab("mean expression \n (subclone)") +
    xlab("expanded cluster") +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 16),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        size = 1
      ),
      axis.text.x = element_text(
        size = 16,
        angle = 90,
        vjust = 0.5
      ),
      axis.text.y = element_text(size = 16),
      axis.title = element_text(size = 16),
      legend.position = "none"
    )

  p1 <- ggplot() +
    geom_line(
      data = blk_long_gj2 %>%
        group_by(pos, subclones) %>%
        summarise(cn = mode(cn)),
      aes(
        x = pos,
        y = cn,
        color = fct_relevel(subclones,
                            gtools::mixedsort(as.character(
                              unique(blk_long_gj2$subclones)
                            )))
      ),
      size = .8
    ) +
    geom_text(
      data = blk_long_gj2   %>%
        filter(gene %in% genes) %>%
        distinct(gene, .keep_all = T),
      aes(x = pos,
          y = 3.7,
          label = gene),
      angle = 90 ,
      size = 8
    ) +
    theme_cowplot() +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      axis.ticks.x = element_blank(),
      axis.title = element_text(size = 20),
      legend.position = "none"
    ) +
    scale_color_manual(name = "cluster", values = colors_vector$subclones) +
    scale_y_continuous(
      breaks = function(x)
        unique(floor(pretty(seq(
          0, (max(x) + 1) * 1.1
        ))))
    ) +
    # ggtitle(active_chr) +
    xlab(paste0("genomic position ", "(", active_chr, ")")) +
    ylab("copy number")


  p2 <-
    mp_l %>%
    ungroup %>%
    mutate(hdb = fct_relevel(subclones,
                             gtools::mixedsort(as.character(
                               unique(cnt_long_gene_trip$subclones)
                             )))) %>%
    ggplot() +
    geom_line(aes(x = pos,
                  y = mean_cluster,
                  color = subclones),
              size = 1.2) +
    scale_color_manual(values = colors_vector$subclones) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
    theme_cowplot() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 20, vjust = 0.5),
      axis.ticks.x = element_blank(),
      axis.title = element_text(size = 20),
      legend.position = "none"
    ) +
    xlab("genomic windows (100 genes)") +
    ylab("expression \n z-score")

  message("plotting.")

  p1 / p2 / p0

}
