#' Return the truncal and branching distances in a data frame
#'
#' @param tree A phylo class object containing the single-cell tree
#'
#' @return A data frame containing the MRCA node, and the truncal and branching branch lengths distances.
#' @export
#'
#' @examples
calc_sctree_dists <- function(tree) {
  # getting info for segment annotation
  tbl_tree <- as_tibble(tree)

  # maximum length will be the truncal jump
  max_length <-
    max(tbl_tree$branch.length[!is.na(tbl_tree$branch.length)])
  dip_label_height <- which(tbl_tree$branch.length == max_length)

  # subsetting a subtree with only nodes after MRCA
  distnodes <- dist.nodes(tree)
  max_mrca_dist <- distnodes %>%
    as_tibble() %>%
    rownames_to_column("node") %>%
    gather(key = "node2", value = "value", -node) %>%
    mutate(node = as.numeric(node),
           node2 = as.numeric(node2)) %>%
    dplyr::filter(node == dip_label_height,
                  node2 >= dip_label_height) %>%
    summarise(max_dist_mrca = max(value)) %>%
    pull(max_dist_mrca)

  # return the truncal and branching dist in a df
  distnodes_df <- tibble(truncal_node = dip_label_height,
                         truncal = max_length,
                         branching = max_mrca_dist)

  return(distnodes_df)

}
