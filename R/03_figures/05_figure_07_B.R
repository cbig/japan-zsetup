# Figure 7B: Boxplots for the inclusive analysis
library(dplyr)
library(ggplot2)
library(gridExtra)
library(hrbrthemes)
library(tibble)
library(tidyr)
library(viridis)
library(zonator)

source("R/00_lib/utils.R")

# Helper functions --------------------------------------------------------

create_boxplot <- function(x) {
  
  boost_label <- function(string) {
    return(paste0("Area protected: ", string))
  }
  
  p1 <- ggplot(x, aes(x = group, y = pr_rem, fill = label)) + 
    geom_boxplot(outlier.alpha = 0.1, outlier.color = "lightgrey",
                 outlier.size = 1, notch = TRUE) + 
    facet_wrap(~fraction, labeller = labeller(fraction = boost_label),
               nrow = 1, ncol = 3) + 
    scale_y_continuous(breaks = seq(0, 1, 0.2),
                       labels = paste0(seq(0, 1, 0.2) * 100, "%"),
                       limits = c(0, 1)) +
    ylab("Fraction of species ranges covered\n") + xlab("") +
    theme_ipsum_rc() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top", legend.direction = "horizontal") + 
    scale_fill_brewer("", type = "qual", palette = 2, direction = 1)
    #scale_fill_viridis("", begin = 0.5, discrete = TRUE)
  return(p1)
}

get_fractions <- function(dat, fractions, label, method, name_groups) {
  frac_list <- list()
  for (fraction in fractions) {
    frac_data <- dat %>% 
      dplyr::select(-cost, -min_pr, -ave_pr, -w_pr, -ext1, -ext2) %>% 
      dplyr::filter(pr_lost >= fraction) %>% 
      tidyr::gather(feature, pr_rem, -pr_lost) %>% 
      dplyr::left_join(., name_groups, by = c("feature" = "feature")) %>% 
      dplyr::mutate(fraction = fraction, label = label, method = method)
    frac_list[[as.character(fraction)]] <- frac_data
  }
  all_dat <- dplyr::bind_rows(frac_list)
  unique_levels <- sort(unique(all_dat$fraction), decreasing = TRUE)
  unique_labels <- paste0(round(100 * (1 - unique_levels), 0), "%")
  all_dat$fraction <- factor(all_dat$fraction, ordered = TRUE,
                             levels = unique_levels,
                             labels = unique_labels)
  return(all_dat)
}

# Read in data ------------------------------------------------------------

# Load the project, creates japan_zproject. WARNING: loaded object is cached
# using memoise-package. If results have change, don't load the cached version.
# Also if the cache needs to updated, you will have to do this manually.
# > key <- list("zsetup")
# > saveCache(zproject_japan, key = key)
zproject_japan <- .load_zproject("zsetup", cache = TRUE, debug = TRUE)

# Extract variant data ----------------------------------------------------

# Reference variants; no masking, but still with condition
#caz_nohm_variant <- zonator::get_variant(zproject_japan, 65)
abf_nohm_variant <- zonator::get_variant(zproject_japan, 66)

#caz_nohm_variant <- regroup_to_taxa(caz_nohm_variant)
abf_nohm_variant <- regroup_to_taxa(abf_nohm_variant)

#caz_nohm_crvs <- zonator::curves(caz_nohm_variant, groups = FALSE)
abf_nohm_crvs <- zonator::curves(abf_nohm_variant, groups = FALSE)

# Hierarchical mask variant for HM3:
# 69_caz_wgt_con_hm3
# 70_abf_wgt_con_hm3
# 
#caz_hm3_variant <- zonator::get_variant(zproject_japan, 69)
abf_hm3_variant <- zonator::get_variant(zproject_japan, 70)

#caz_hm3_variant <- regroup_to_taxa(caz_hm3_variant)
abf_hm3_variant <- regroup_to_taxa(abf_hm3_variant)

#caz_hm3_crvs <- zonator::curves(caz_hm3_variant, groups = FALSE)
abf_hm3_crvs <- zonator::curves(abf_hm3_variant, groups = FALSE)

name_groups <- tibble::tibble(feature = zonator::featurenames(abf_hm3_variant),
                              group = abf_hm3_variant@groups$name) 

# Get representation levels -----------------------------------------------

fractions <- c(0.978, 0.906, 0.83)

# CAZ

#caz_nohm_all_dat <- get_fractions(caz_nohm_crvs, fractions = fractions,
#                                  label = "Without PAN", method = "CAZ", 
#                                  name_groups = name_groups)

#caz_hm3_all_dat <- get_fractions(caz_hm3_crvs, fractions = fractions,
#                                 label = "With PAN", method = "CAZ",
#                                 name_groups = name_groups)

#caz_all_dat <- dplyr::bind_rows(caz_nohm_all_dat, caz_hm3_all_dat)
#caz_all_dat$label <- factor(caz_all_dat$label, ordered = TRUE,
#                            levels = c("Without PAN", "With PAN"))

# ABF

abf_nohm_all_dat <- get_fractions(abf_nohm_crvs, fractions = fractions,
                                  label = "Without PAN", method = "ABF",
                                  name_groups = name_groups)

abf_hm3_all_dat <- get_fractions(abf_hm3_crvs, fractions = fractions,
                                 label = "With PAN", method = "ABF",
                                 name_groups = name_groups)

abf_all_dat <- dplyr::bind_rows(abf_nohm_all_dat, abf_hm3_all_dat)
abf_all_dat$label <- factor(abf_all_dat$label, ordered = TRUE,
                            levels = c("Without PAN", "With PAN"))

# Print out means for the different taxa
abf_all_dat %>% 
  dplyr::filter(fraction == "9%") %>% 
  dplyr::group_by(group) %>% 
  summarise(
    avg_rep = mean(pr_rem),
    median_rep = median(pr_rem),
    min_rep = min(pr_rem),
    max_rep = max(pr_rem)
  )
abf_all_dat %>% 
  dplyr::filter(fraction == "17%") %>% 
  dplyr::group_by(group) %>% 
  summarise(
    avg_rep = mean(pr_rem),
    median_rep = median(pr_rem),
    min_rep = min(pr_rem),
    max_rep = max(pr_rem)
  )

# Put everything together
#all_dat <- dplyr::bind_rows(caz_all_dat, abf_all_dat)

# Plot --------------------------------------------------------------------

#p1 <- create_boxplot(caz_all_dat)
p2 <- create_boxplot(abf_all_dat)
  
img_height <- 10
img_width <- 18

ggsave("figs/figure_07/figure_07_B.png", p2, width = img_width, height = img_height, 
       units = "cm", dpi = 600)
