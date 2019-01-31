# Figure 4: Performance curves for the inclusice analysis
library(dplyr)
library(ggplot2)

source("R/00_lib/utils.R")

# Load the project, creates japan_zproject. WARNING: loaded object is cached
# using memoise-package. If results have change, don't load the cached version.
# Also if the cache needs to updated, you will have to do this manually.
zproject_japan <- .load_zproject("zsetup", cache = FALSE, debug = TRUE)

pa_breaks <- c(0.02, 0.09, 0.17)

# Parametric (mean + SD)
p1 <- get_stat_curves(zproject_japan, variant_ids = c(61, 63, 65)) #%>% 
plot_curves(p1, title = "CAZ", non_param = FALSE, invert_x = TRUE,
              highlights = pa_breaks,
              labels = c("Baseline", "Weights", "Condition"))

p2 <- get_stat_curves(zproject_japan, variant_ids = c(62, 64)) %>% 
  plot_curves(title = "ABF", non_param = FALSE, invert_x = TRUE,
              nrow = 1, ncol = 2, highlights = pa_breaks,
              labels = c("Baseline", "Weights"))

p3 <- get_stat_curves(zproject_japan, variant_ids = c(65, 69)) %>% 
  plot_curves(title = "", non_param = FALSE, invert_x = TRUE,
              nrow = 1, ncol = 2, highlights = pa_breaks)

p4 <- get_stat_curves(zproject_japan, variant_ids = c(66, 70)) %>% 
  plot_curves(title = "", non_param = FALSE, invert_x = TRUE,
              plot_min = TRUE, plot_max = TRUE,
              nrow = 1, ncol = 2, highlights = pa_breaks,
              labels = c("Without PAN", "With PAN"))

# Non-parametric (median + quartiles)
p5 <- get_stat_curves(zproject_japan, variant_ids = c(61, 63)) %>% 
  plot_curves(title = "CAZ", non_param = TRUE, invert_x = TRUE,
              nrow = 1, ncol = 2, highlights = pa_breaks,
              labels = c("Baseline", "Weights"))

p6 <- get_stat_curves(zproject_japan, variant_ids = c(62, 64)) %>% 
  plot_curves(title = "ABF", non_param = TRUE, invert_x = TRUE,
              nrow = 1, ncol = 2, highlights = pa_breaks,
              labels = c("Baseline", "Weights"))

p7 <- get_stat_curves(zproject_japan, variant_ids = c(65, 69)) %>% 
  plot_curves(title = "", non_param = TRUE, invert_x = TRUE,
              nrow = 1, ncol = 2, highlights = pa_breaks)

p8 <- get_stat_curves(zproject_japan, variant_ids = c(66, 70)) %>% 
  plot_curves(title = "", non_param = TRUE, invert_x = TRUE,
              nrow = 1, ncol = 2, highlights = pa_breaks, 
              labels = c("Without PAN", "With PAN"))

p9 <- get_stat_curves(zproject_japan, variant_ids = c(62, 61, 64, 63, 66, 65)) %>% 
  plot_curves(title = "", non_param = TRUE, invert_x = TRUE,
              nrow = 3, ncol = 2, highlights = pa_breaks,
              labels = c("ABF", "CAZ", "ABF / Weights", "CAZ / Weights",
                         "ABF / Weights / HII", "CAZ / weights / HII"))

# Save images ---------------------------------------------------------

img_height <- 11
img_width <- 18

ggsave("figs/figure_04_1_caz_mean.png", p1, width = img_width, height = img_height, 
       units = "cm")
ggsave("figs/figure_04_2_abf_mean.png", p2, width = img_width, height = img_height, 
       units = "cm")
ggsave("figs/figure_04_3_caz_hii_mean.png", p3, width = img_width, height = img_height, 
       units = "cm")
ggsave("figs/figure_04_4_abf_hii_mean.png", p4, width = img_width, height = img_height, 
       units = "cm")
ggsave("figs/figure_04_5_caz_median.png", p5, width = img_width, height = img_height, 
       units = "cm")
ggsave("figs/figure_04_6_abf_median.png", p6, width = img_width, height = img_height, 
       units = "cm")
ggsave("figs/figure_04_7_caz_hii_median.png", p7, width = img_width, height = img_height, 
       units = "cm")
ggsave("figs/figure_04_8_abf_hii_median.png", p8, width = img_width, height = img_height, 
       units = "cm")

ggsave("figs/figure_S04.png", p9, width = img_width, height = 28, 
       units = "cm")
