# Figure 7: Performance curves for the per-taxon analyses
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
variant_ids <- c(66, 70)

# Get exact figures
perf_levels <- get_perf_levels(zproject_japan, variant_ids, pa_breaks[1:3])

# Plot
p1 <- get_stat_curves(zproject_japan, variant_ids) %>% 
  plot_curves(title = "", non_param = FALSE, invert_x = TRUE,
              highlights = pa_breaks, max_y = 0.6,
              labels = c("Without PAN", "With PAN"))

# Save images ---------------------------------------------------------

img_height <- 10
img_width <- 18

ggsave("figs/figure_07/figure_07_A.png", p1, width = img_width, height = img_height, 
       units = "cm", dpi = 600)
