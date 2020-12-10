#' Plots the resulting umap and colors with the clustering information
#'
#' @param umap_df
#' @param clustering
#'
#' @return A ggplot object with a umap scatterplot.
#' @export
#'
#' @examples
plot_umap <- function(umap_df,
                      clustering) {

  umap_df <- umap_df %>%
    dplyr::rename(cells = "cell") %>%
    left_join(clustering)

  # make sure to run setup.R to obtain the colors vector
  colors_vector_gg <- c(colors_vector$subclones, colors_vector$superclones)

  # umap plot
  umap_p <- ggplot(umap_df) +
    geom_point(aes(x = V1, y = V2, colour = superclones), alpha = 1, size = 8) +
    geom_point(aes(x = V1, y = V2, colour = subclones), alpha = 1) +
    scale_color_manual(values = colors_vector_gg) +
    theme_classic() +
    theme(axis.title.x=element_text(size = 16),
          axis.text.x= element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y = element_text(size = 16),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          panel.border = element_rect(color = "black",
                                      fill = NA,
                                      size = 2))  +
    xlab("umap1") +
    ylab("umap2")

  print(umap_p)

  return(umap_p)

}

