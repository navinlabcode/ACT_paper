#' Calculates the breadth of coverage for samples
#'
#' @param path Filepath with the location of .covhist files from genomeBedCoverage
#'
#' @return
#' @export
#'
#' @examples
calc_coverage <- function(path) {

  inpaths <- Sys.glob(paste0(path, "*.covhist.txt"))

  coverage.stats <- tibble(bed_path=inpaths) %>%
    mutate(cellname = str_extract(basename(bed_path), "^[^.]*")) %>%
    group_by(cellname) %>%
    summarize(.groups="keep",
              read_tsv(bed_path,
                       col_names=c("refname", "depth", "count",
                                   "refsize", "frac"),
                       col_types=cols(col_character(), col_double(),
                                      col_double(), col_double(),
                                      col_double())),
    ) %>%
    filter(refname=="genome") %>%
    summarize(breadth = 1 - frac[depth==0],
              .groups="keep")

}
