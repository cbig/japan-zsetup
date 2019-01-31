# Table S4: Japanese prefectures
library(dplyr)
library(ggplot2)
library(hrbrthemes)
library(zonator)

source("R/00_lib/utils.R")

# Load the project, creates japan_zproject. WARNING: loaded object is cached
# using memoise-package. If results have change, don't load the cached version.
# Also if the cache needs to updated, you will have to do this manually.
# > key <- list("zsetup")
# > saveCache(zproject_japan, key = key)
zproject_japan <- .load_zproject("zsetup", cache = TRUE, debug = TRUE)

# Get the PPA data for variant 66_abf_wgt_con
lsm_data66 <- get_ppa_data(zproject_japan, 
                           prefectures_file = "../Data.150928/prefectures.csv",
                           66, "abf_wgt_con")

# Plot
p1 <- ggplot(lsm_data66, aes(x = mrank, y = log(sdsum), label = prefecture)) +
  geom_point(aes(size = area, alpha = .2)) +
  geom_text(hjust = 0.5, size = 3, nudge_y = 0.1, check_overlap = FALSE) +
  scale_size(range = c(1, 10)) + xlim(c(0.2, 1)) +
  ylab("log(species distribution sum)") + xlab("Mean rank") +
  theme_ipsum_rc() + theme(legend.position = "none")

ggsave("figs/figure_S04.png", p1, width = 18, height = 18, units = "cm")
