# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Tue Nov 24 17:06:17 2020
# setup file
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Tue Nov 24 17:06:25 2020

# Set up most of the things
# - Load packages
# - helpful functions
# - gene annotation for hg19

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Tue Nov 24 17:07:10 2020
# Load packages
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Tue Nov 24 17:07:15 2020
library(paletteer)
library(here)
library(rstatix)
library(parallel)
library(patchwork)
library(umap)
library(ggbeeswarm)
library(phangorn)
library(uwot)
library(tidyverse)
library(ComplexHeatmap)
library(amap)
library(cowplot)
library(ggtree)
library(ape)
library(future.apply)
library(tximport)
library(GenomicFeatures)
library(AnnotationHub)
library(jsonlite)
library(readr)
library(DESeq2)
library(Homo.sapiens)
library(conflicted)
library(vcfR)
library(jcolors)
library(janitor)
library(boot)
library(biomaRt)
library(jsonlite)
library(readr)
library(ComplexHeatmap)
library(DESeq2)
library(patchwork)
library(fgsea)
library(AnnotationDbi)
library(org.Hs.eg.db)
theme_set(theme_cowplot())
plan(multiprocess, workers = 20)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Tue Nov 24 17:17:44 2020
# resolving most frequent conflicts
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Tue Nov 24 17:17:52 2020
conflicted::conflict_prefer("which", "Matrix")
conflicted::conflict_prefer("desc", "dplyr")
conflicted::conflict_prefer("expand", "tidyr")
conflicted::conflict_prefer("layout", "plotly")
conflicted::conflict_prefer("ladderize", "ape")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("count", "dplyr")
conflicted::conflict_prefer("group_by", "dplyr")
conflicted::conflict_prefer("rename", "dplyr")
conflicted::conflict_prefer("select", "dplyr")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Tue Nov 24 17:09:02 2020
# set up helpful functions and preferred ggplot theme
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Tue Nov 24 17:09:09 2020

'%!in%' <- function(x,y)!('%in%'(x,y))
head2 <- function(x) head(x)[,1:5]
theme_set(theme_cowplot())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Sat Dec 12 14:03:04 2020
# cna pipeline lib
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Sat Dec 12 14:03:12 2020

bins_in_cna_pipeline <- read.table("extdata/lib/bins_in_cna_pipeline.txt", header = T)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Tue Nov 24 17:10:51 2020
# hg19 genes from Homo Sapiens package
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Tue Nov 24 17:11:11 2020

txdb <- Homo.sapiens
hg19_genes <- GenomicFeatures::genes(txdb, columns = "SYMBOL")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Wed Dec  9 09:55:57 2020
# colors setup
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Wed Dec  9 09:56:04 2020

get_colors <- function() {

  # This function is an unfortunate inheritance from the beginning of this project. Everything uses it
  # so it is was kept.

   # subclones
  hues1 <- paletteer::paletteer_d("ggsci::default_igv",
                                  n = 35) %>% unclass()


  colors_vec1 <- setNames(hues1, paste0("c", 1:35))

  #superclones

  colors_g_cl <-
    structure(unclass(paletteer::paletteer_d("yarrr::info", n = 9)),
              names = paste0("s", 1:9))
  colors_g_cl[1] <- "#E7A79BFF"
  colors_g_cl[2] <- "#90A8C0FF"
  colors_g_cl[6] <- "#B4DCF5FF"
  colors_g_cl[7] <- "#F2C695FF"

  colors_list <- list(subclones = colors_vec1,
                      superclones = colors_g_cl)

  return(colors_list)

}

colors_vector <- get_colors()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Dec 11 13:31:13 2020
# Ploidy Scale
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Dec 11 13:31:20 2020

# this function ploidy scale and rounds value to integer numbers
ploidy_scale <- function(ploidy_VAL, df, round = TRUE) {

  # population segmenter returns log ratios, transforming back to ratios values
  popseg_long_exp <- as.data.frame(2^df)

  # correcting to ploidy
  popseg_ploidy_cor <- as.data.frame(popseg_long_exp * ploidy_VAL)

  if (round == TRUE) {
    # rounding to the nearest integer values
    popseg_round <- as.data.frame(round(popseg_ploidy_cor, 0))
    return(popseg_round)
  } else return(popseg_ploidy_cor)

}

