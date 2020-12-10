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
library(Rphenograph)
library(aplot)
library(parallel)
library(umap)
library(phangorn)
library(uwot)
library(tidyverse)
library(ComplexHeatmap)
library(cowplot)
library(ggtree)
library(ape)
library(future.apply)
library(tximport)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(jsonlite)
library(readr)
library(DESeq2)
library(Homo.sapiens)
library(conflicted)

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

