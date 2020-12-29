
#' Parsing mutation timer vcf to capture clonal/subclonal status
#'
#' @param vcf_file Vcf file from mutation timer
#'
#' @return A data frame with the parsed vcf
#' @export
#'
#' @examples
parse_mtimer_vcf <- function(vcf_file) {

  skip_lines <- system(paste0("grep '#' ", vcf_file, " | wc -l"), intern = T)

  vcf <- read_delim(vcf_file,
                    delim = "\t",
                    skip = as.numeric(skip_lines) - 1,
                    col_types = cols(`#CHROM` = 'c'))


  vcf <- vcf %>%
    mutate(CLS = str_extract(INFO, "CLS(.*)"),
           mutation_time = str_remove(CLS, "CLS=")) %>%
    separate(mutation_time, into = c("mutation_time", "mutation_time_class"), sep = " ") %>%
    mutate(mutation_time_class = str_remove(mutation_time_class, "\\["),
           mutation_time_class = str_remove(mutation_time_class, "\\]")) %>%
    dplyr::rename(CHROM = `#CHROM`) %>%
    mutate(CHROM = paste0("chr", CHROM)) %>%
    select(CHROM, POS, mutation_time, mutation_time_class)

  return(vcf)

}
