# Figure XX: Feature coverage distributions.
library(dplyr)
library(plotly)
library(purrr)
library(tidyr)
library(zonator)

source("R/00_lib/utils.R")

# Load the project, creates japan_zproject. WARNING: loaded object is cached
# using memoise-package. If results have change, don't load the cached version.
# Also if the cache needs to updated, you will have to do this manually.
zproject_japan <- .load_zproject("zsetup", cache = TRUE, debug = TRUE)

breaks <- seq(0, 1, 0.1)

variant66 <- zonator::get_variant(zproject_japan, 66)
variant66_curves <- zonator::curves(variant66) %>% 
  dplyr::select(-cost, -min_pr, -ave_pr, -w_pr, -ext1, -ext2) %>% 
  dplyr::mutate(bin = ntile(pr_lost, 10)) %>% 
  tidyr::gather(feature, pr_rem, -pr_lost, -bin) 

bin_vectors <- function(id, dat) {
  dat <- dplyr::filter(dat, bin == id)
  freqs <- cut(dat$pr_rem, breaks = breaks, include.lowest = TRUE)
  return(data.frame(pr_lost = dat$pr_lost, pr_lost_bin = dat$bin, bin = freqs))
}  

get_frequency <- function(id, dat) {
  dat <- dat[[id]]
  freq_df <- data.frame(pr_lost_bin = id, freq = table(dat$bin))
  names(freq_df) <- c("x", "y", "z")
  return(freq_df)
}

bins <- unique(variant66_curves$bin)

dists_df <- purrr::map(bins, bin_vectors, dat = variant66_curves)
dists <- purrr::map(bins, get_frequency, dat = dists_df) %>% 
  dplyr::bind_rows() %>% 
  dplyr::mutate(x = x / 10, y = as.numeric(y) / 10)

dists %>%
  group_by(x) %>% 
  plot_ly(x = dists$y, y = dists$x, z = dists$z, type = "scatter3d", mode = "lines") 
