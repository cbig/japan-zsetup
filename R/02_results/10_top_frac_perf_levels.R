library(ggplot2)
library(hrbrthemes)
library(tidyverse)
library(zonator)

source("R/00_lib/utils.R")

# Helper functions --------------------------------------------------------

create_boxplot <- function(x, title = "") {
  
  p1 <- ggplot(x, aes(x = type, y = pr_rem, fill = group)) + 
    geom_boxplot(color = "darkgrey", outlier.alpha = 0.1, 
                 outlier.color = "lightgrey", outlier.size = 1) + 
    scale_y_continuous(breaks = seq(0, 1, 0.2),
                       labels = paste0(seq(0, 1, 0.2) * 100, "%"),
                       limits = c(0, 1)) +
    ylab("Fraction of species ranges covered\n") + xlab("") +
    theme_ipsum_rc() + ggtitle(title) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right", legend.direction = "vertical") + 
    #scale_fill_brewer("", type = "qual", palette = 1, direction = 1) 
    scale_fill_viridis("", begin = 0.0, discrete = TRUE) + 
    coord_flip()
  return(p1)
}

get_perf_level <- function(zproject, variant_id, fraction, group_label) {
  crvs_data <- zonator::get_variant(zproject, variant_id) %>% 
    zonator::results() %>%
    zonator::performance(pr.lost = 1 - top_fraction) %>% 
    tidyr::gather(feature, pr_rem, -pr_lost) %>% 
    dplyr::mutate(group = group_label) 
  return(crvs_data)
}

get_grp_perf_levels <- function(zproject, variant_id, fraction) {
  variant <- zonator::get_variant(zproject_japan, variant_id) %>% 
    regroup_to_taxa()
  grp_perf_levels <- variant  %>% 
    zonator::results() %>% 
    zonator::performance(pr.lost = 1 - top_fraction) %>% 
    tidyr::gather(feature, pr_rem, -pr_lost)
  
  name_groups <- tibble::tibble(feature = zonator::featurenames(variant),
                                group = variant@groups$name) 
  
  grp_perf_levels <- dplyr::left_join(grp_perf_levels, 
                                      name_groups, 
                                      by = c("feature" = "feature")) 
  return(grp_perf_levels)
}


# Read in data ------------------------------------------------------------

# Load the project, creates japan_zproject. WARNING: loaded object is cached
# using memoise-package. If results have change, don't load the cached version.
# Also if the cache needs to updated, you will have to do this manually.
# > key <- list("zsetup")
# > saveCache(zproject_japan, key = key)
zproject_japan <- .load_zproject("zsetup", cache = TRUE, debug = TRUE)

# Extract variant data ----------------------------------------------------

#top_fraction <- 0.0394
#top_fraction <- 0.17
top_fraction <- 1 - 0.906
# NO CONDITION

# Following variants are used:
#
# Taxon-specific (ts):
#
# 04_amp_abf_wgt
# 14_bir_abf_wgt
# 24_frf_abf_wgt
# 34_mam_abf_wgt
# 44_pla_abf_wgt
# 54_rep_abf_wgt
variant_ids <- c(4, 14, 24, 34, 44, 54)
variant_names <- c("Amphibians", "Birds", "Freshwater_fish", "Mammals", 
                   "Plants", "Reptiles")
ts_perf_levels <- purrr::map2(variant_ids, variant_names, get_perf_level,
                              zproject = zproject_japan, fraction = top_fraction) %>% 
  dplyr::bind_rows()
ts_perf_levels$type <- "Taxon-specific"

# Inclusive (joint):
# 
# 64_abf_wgt
incl_perf_levels <- get_grp_perf_levels(zproject_japan, 64, top_fraction)
incl_perf_levels$type <- "Inclusive"

# Pre-loaded based on 3.94% hierarchical mask
# 
# 71_abf_wgt_hm_st_top
pl_perf_levels <- get_grp_perf_levels(zproject_japan, 71, top_fraction)
pl_perf_levels$type <- "Taxon-specific combined"

# Bind everything together
all_dat <- dplyr::bind_rows(ts_perf_levels, incl_perf_levels, 
                            pl_perf_levels)
all_dat$type <- factor(all_dat$type, ordered = TRUE,
                       levels = c("Inclusive", 
                                  "Taxon-specific combined", 
                                  "Taxon-specific"))
all_dat$group <- factor(all_dat$group, ordered = TRUE,
                        levels = variant_names)

# Print out means for types 
all_dat %>% 
  dplyr::group_by(type) %>% 
  summarise(
    avg_rep = mean(pr_rem),
    median_rep = median(pr_rem),
    min_rep = min(pr_rem),
    max_rep = max(pr_rem)
  )

# Print out means for the different taxa
all_dat %>% 
  dplyr::filter(type == "Inclusive") %>% 
  dplyr::group_by(group) %>% 
  summarise(
    avg_rep = mean(pr_rem),
    median_rep = median(pr_rem),
    min_rep = min(pr_rem),
    max_rep = max(pr_rem)
  )

# Plot --------------------------------------------------------------------

p1 <- create_boxplot(all_dat, title = "Top 17%")

img_height <- 18
img_width <- 18

#ggsave("figs/figure_05_caz.png", p1, width = img_width, height = img_height, 
#       units = "cm")
ggsave("figs/ts_top_00394_boxplot.png", p1, width = img_width, height = img_height, 
       units = "cm")

