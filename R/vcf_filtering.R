vcf_filtering <- function(mutect,
                          annovar,
                          bulk_normal,
                          bulk_tumor,
                          min_DP = 3,
                          min_AD = 1,
                          min_AF_bulk = 0.1,
                          min_DP_bulk = 10,
                          min_AD_bulk = 5,
                          write_vcf = FALSE,
                          mutect_vcf_path = NULL) {

  gt_table <- extract.gt(mutect) %>%
    as.data.frame() %>%
    rownames_to_column("chr_pos") %>%
    separate(col = chr_pos,
             into = c("CHROM",
                      "POS")) %>%
    gather(key = "sample",
           value = "GT",
           -CHROM,
           -POS)

  dp_table <- extract.gt(mutect,
                         element = "DP",
                         as.numeric = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column("chr_pos") %>%
    separate(col = chr_pos,
             into = c("CHROM",
                      "POS")) %>%
    gather(key = "sample",
           value = "DP",-CHROM,-POS)

  ad_table <- extract.gt(mutect,
                         element = "AD",
                         as.numeric = F) %>%
    as.data.frame() %>%
    rownames_to_column("chr_pos") %>%
    separate(col = chr_pos,
             into = c("CHROM",
                      "POS")) %>%
    gather(key = "sample",
           value = "AD",-CHROM,-POS) %>%
    separate(col = AD,
             into = c("REF_AD",
                      "ALT_AD")) %>%
    mutate(REF_AD = as.numeric(REF_AD),
           ALT_AD = as.numeric(ALT_AD))

  af_table <- extract.gt(mutect,
                         element = "AF",
                         as.numeric = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column("chr_pos") %>%
    separate(col = chr_pos,
             into = c("CHROM",
                      "POS")) %>%
    gather(key = "sample",
           value = "AF",-CHROM,-POS)


  vcf_table <- left_join(gt_table,
                         af_table) %>%
    left_join(dp_table) %>%
    left_join(ad_table) %>%
    mutate(POS = as.numeric(POS))

  # simplifying the annovar table
  annovar_red <- annovar %>%
    select(CHROM, POS, Func.refGene, snp129, Gene.refGene, cosmic70, ExonicFunc.refGene) %>%
    filter(ExonicFunc.refGene != "unknown")

  # rescuing tp53 from pD which has a off by one error
  if (bulk_tumor == "BRCADT") {
    annovar_red <- annovar_red %>%
      mutate(POS = case_when(Gene.refGene == "TP53" ~ POS -1,
                             TRUE ~ POS))
  }


  vcf_table_annovar <- vcf_table %>%
    left_join(annovar_red) %>%
    filter(!str_detect(snp129, "rs"),
           Func.refGene == "exonic") %>%
    mutate(id = paste0(Gene.refGene,":",POS),
           chr_pos = paste(CHROM, POS, sep = "_")) %>%
    group_by(CHROM) %>%
    mutate(pos_dist = abs(lead(POS) - POS)) %>%
    ungroup() %>%
    mutate(AF = ifelse(ALT_AD < min_AD, NA, AF),
           DP = ifelse(DP == 0, NA, DP),
           AD = ifelse(ALT_AD < min_AD, NA, ALT_AD),
           class = case_when(str_detect(sample, !!bulk_tumor) ~ "tumor",
                             str_detect(sample, !!bulk_normal) ~ "normal",
                             TRUE ~ "sc pseudo-bulk"))

  # AD bulk filter
  ad_bulk_filter <- vcf_table_annovar %>%
    filter(str_detect(sample, !!bulk_tumor),
           ALT_AD < min_AD_bulk) %>%
    distinct(id) %>%
    pull(id)

  # DP bulk filter
  dp_bulk_filter <- vcf_table_annovar %>%
    filter(str_detect(sample, !!bulk_tumor),
           DP < min_DP_bulk) %>%
    distinct(id) %>%
    pull(id)

  # removing high allele frequency in normal
  high_normal_af <- vcf_table_annovar %>%
    filter(sample == !!bulk_normal,
           AF > 0.05) %>%
    pull(Gene.refGene)

  # filter for AF in the tumor
  af_bulk_filter <- vcf_table_annovar %>%
    filter(str_detect(sample, !!bulk_tumor),
           AF < min_AF_bulk) %>%
    distinct(id) %>%
    pull(id)

  # removing calls too close to one another
  gene_closepos_calls <- vcf_table_annovar %>%
    filter(pos_dist < 1000 & pos_dist != 0) %>%
    pull(Gene.refGene)

  # removing genes according to filtering
  vcf_table_annovar <- vcf_table_annovar %>%
    filter(Gene.refGene %!in% gene_closepos_calls,
           Gene.refGene %!in% high_normal_af,
           id %!in% ad_bulk_filter,
           id %!in% dp_bulk_filter,
           id %!in% af_bulk_filter)

  vcf_table_annovar <- vcf_table_annovar %>%
    filter(ALT_AD >= min_AD)

  if (write_vcf == TRUE) {

    mutect_vcf_file <- paste0(mutect_vcf_path, "/tumor.filtered.pass.vcf")

    skip_lines <- system(paste0("grep '#' ", mutect_vcf_file, " | wc -l"), intern = T)

    mutect_vcf_f <- read_delim(mutect_vcf_file,
                               delim = "\t",
                               skip = as.numeric(skip_lines) - 1)


    mutect_vcf_filt <- mutect_vcf_f %>%
      filter(POS %in% vcf_table_annovar$POS)

    write.table(mutect_vcf_filt, paste0(mutect_vcf_path, "tumor.filtered.pass", min_AF_bulk,".filtered.noheader.vcf"),
                col.names = F, quote = FALSE, sep = "\t", row.names = F)


  }

  return(vcf_table_annovar)

}
