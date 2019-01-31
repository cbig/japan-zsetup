# Figure 3: Performance curves for the per-taxon analyses
library(dplyr)
library(ggplot2)
library(gridExtra)

source("R/00_lib/utils.R")

# Load the project, creates japan_zproject. WARNING: loaded object is cached
# using memoise-package. If results have change, don't load the cached version.
# Also if the cache needs to updated, you will have to do this manually.
# > key <- list("zsetup")
# > saveCache(zproject_japan, key = key)
zproject_japan <- .load_zproject("zsetup", cache = TRUE, debug = TRUE)

pa_breaks <- c(0.022, 0.094, 0.17, 0.5, 0.75, 1.0)

labels <- c("Variant 1: Baseline", "Variant 2: Weights", 
            "Variant 3: Condition", "Variant 4: With current PAN")

# ABF
p1 <- get_stat_curves(zproject_japan, variant_ids = c(62, 64, 66, 70)) %>% 
  plot_curves(title = "ABF", non_param = FALSE, invert_x = TRUE,
              highlights = pa_breaks, max_y = 1.0,
              labels = labels)

# CAZ
#p2 <- get_stat_curves(zproject_japan, variant_ids = c(61, 63, 65, 69)) %>% 
#  plot_curves(title = "CAZ", non_param = FALSE, invert_x = TRUE,
#              highlights = pa_breaks, max_y = 1.0,
#              labels = labels)

#p3 <- grid.arrange(p1, p2, nrow = 1, ncol = 2)

# Save images ---------------------------------------------------------

ggsave("figs/figure_S03.png", p1, width = 18, height = 12, units = "cm")
