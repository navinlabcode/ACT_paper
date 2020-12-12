#' Runs a minimal evolution tree algorithm for the consensus data frame
#'
#' @param consensus_df The consensus data frame from the segmented dataset
#' @param clusters The results from clustering
#' @param ploidy_VAL The inferred ploidy value. Will be used to scale the data
#' @param rotate_nodes Nodes to be rotated
#'
#' @return
#' @export
#'
#' @examples
run_me_tree <- function(consensus_df,
                        clusters,
                        ploidy_VAL,
                        rotate_nodes = NULL) {

  consensus_int <- ploidy_scale(ploidy_VAL, consensus_df)

  #adding a neutral state, will use as root
  consensus_int[nrow(consensus_int) + 1, ] <- round(ploidy_VAL)
  consensus_int[nrow(consensus_int) + 1, ] <- round(ploidy_VAL)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Nov 20 12:24:27 2020
  # tree ME
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Nov 20 12:24:35 2020

  tree <- ape::fastme.bal(dist(consensus_int, method = "manhattan"))

  tree <-
    root.phylo(tree,
               outgroup = which(tree$tip.label == Ntip(tree)),
               resolve.root = T)

  tree <-
    drop.tip(tree, tip = as.character(c(
      nrow(consensus_int), nrow(consensus_int) - 1
    )))

  tree <- ladderize(tree)

  # getting order
  is_tip <- tree$edge[, 2] <= length(tree$tip.label)
  ordered_tips_index <- tree$edge[is_tip, 2]
  tree_tips_order <- tree$tip.label[ordered_tips_index] %>% rev()

  # adding superclones information to frequencies_df
  clones_df <- clusters %>%
    distinct(superclones, subclones) %>%
    mutate(taxa = subclones)

  clones_df <- clones_df[c("taxa", "superclones","subclones")]

  if (!is.null(rotate_nodes)) {
    tree <-
      phytools::rotateNodes(tree, rotate_nodes)
    is_tip <-
      tree$edge[, 2] <= length(tree$tip.label)
    ordered_tips_index <- tree$edge[is_tip, 2]
    tree_tips_order <-
      tree$tip.label[ordered_tips_index] %>% rev()
    tree_tips_order <-
      tree_tips_order[str_detect(tree_tips_order, "c")]
  }

  p <- ggtree(tree, ladderize = F, size = 2) +
    geom_treescale()

  colors_vector_gg <-
    c(colors_vector$subclones, colors_vector$superclones)

  p <- p %<+% clones_df +
    geom_tippoint(aes(color = superclones),
                  size = 10) +
    geom_tippoint(aes(color = subclones),
                  size = 3) +
    scale_color_manual(values = colors_vector_gg) +
    theme(legend.position = "none")

  print(p)

  results <- list(tree = tree,
                  cs_plot = p,
                  cs_tree_order = tree_tips_order)

  return(results)

}
