# Figure 5: Performance curves for the per-taxon analyses
library(dplyr)
library(ggplot2)

source("R/00_lib/utils.R")

# Load the project, creates japan_zproject. WARNING: loaded object is cached
# using memoise-package. If results have change, don't load the cached version.
# Also if the cache needs to updated, you will have to do this manually.
# > key <- list("zsetup")
# > saveCache(zproject_japan, key = key)
zproject_japan <- .load_zproject("zsetup", cache = TRUE, debug = TRUE)

pa_breaks <- c(0.022, 0.094, 0.17, 0.5, 0.75, 1.0)

# Parametric (mean + SD)
crvs <- get_stat_curves(zproject_japan, variant_ids = c(4, 14, 24, 34, 44, 54)) %>% 
  dplyr::mutate(step = factor("taxon-specific", levels = c("taxon-specific", "inclusive")))
p1 <- crvs %>% 
  plot_curves(title = "", non_param = FALSE, invert_x = TRUE,
              highlights = pa_breaks,
              labels = c("Amphibians", "Birds", "Freshwater fish", "Mammals",
                         "Plants", "Reptiles"))

v64 <- zonator::get_variant(zproject_japan, 64)
v64 <- regroup_to_taxa(v64)
v64_crvs <- calc_row_stats(v64, groups = TRUE) %>% 
  dplyr::mutate(variant = factor(variant, levels = unique(variant),
                                 ordered = TRUE)) %>% 
  dplyr::mutate(step = factor("inclusive", levels = c("taxon-specific", "inclusive")))

p2 <- v64_crvs %>% 
  plot_curves(title = "", non_param = FALSE, invert_x = TRUE,
              highlights = pa_breaks,
              labels = c("Amphibians", "Birds", "Freshwater fish", "Mammals",
                         "Plants", "Reptiles"))

# Extra -------------------------------------------------------------------

all_crvs <- dplyr::bind_rows(crvs, v64_crvs)

p3 <- ggplot(all_crvs, aes(x = pr_lost, y = mean), linetype = step) +
  geom_line() + facet_wrap(~ variant)

# Save images ---------------------------------------------------------

img_height <- 10
img_width <- 18

ggsave("figs/figure_05.png", p1, width = img_width, height = img_height, 
       units = "cm", dpi = 600)
