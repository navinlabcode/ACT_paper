#' Plot the single cell tree
#'
#' @param tree A phylo class object containing the single-cell tree
#' @param anno_y Controls where the annotation with distances will be plot
#' @param title Character. Optional to add a title to the plot
#'
#' @return
#' @export
#'
#' @examples
plot_sctree <- function(tree,
                        anno_y = 1200,
                        title = NULL) {

  dist_nodes <- calc_sctree_dists(tree)

  p <-
    ggtree::ggtree(ape::ladderize(tree),
                   ladderize = F,
                   size = .2) +
    geom_tippoint(color = "black", size = 1) +
    theme_tree2(axis.text.x = element_text(size = 15),
                axis.title.x = element_text(size = 15)) +
    # punctuated
    annotate(
      "segment",
      x = 0,
      xend = dist_nodes$truncal - 30,
      y = anno_y,
      yend = anno_y
    ) +
    annotate(
      "text",
      x = dist_nodes$truncal / 2,
      y = anno_y + 30,
      size = 5,
      label = round(dist_nodes$truncal)
    ) +
    annotate(
      "text",
      x = dist_nodes$truncal / 2,
      y = anno_y - 30,
      size = 5,
      label = "truncal"
    ) +
    # MRCA to last node
    annotate(
      "segment",
      x = dist_nodes$truncal + 30,
      xend = dist_nodes$truncal + dist_nodes$branching,
      y = anno_y,
      yend = anno_y
    ) +
    annotate(
      "text",
      x = (dist_nodes$truncal) + (dist_nodes$branching / 2),
      y = anno_y + 30,
      size = 5,
      label = round(dist_nodes$branching)
    ) +
    annotate(
      "text",
      x = (dist_nodes$truncal) + (dist_nodes$branching / 2),
      y = anno_y - 30,
      label = "branching",
      size = 5
    ) +
    annotate(
      "text",
      x = 1400,
      y = dist_nodes$truncal_node,
      label = "diploid",
      size = 5
    ) +
    xlab("manhattan distance")

  p

  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }

  return(p)

}
