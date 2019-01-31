library(ggplot2)
library(zonator)

japan_project <- load_zproject(root = "zsetup/", debug = TRUE)

all_variants <- variants(japan_project)

# TEMP: not all variants are run yet, so let's do some subsetting
ready_variants <- c(1, 2, 6, 7, 11, 12, 16, 17, 21, 22, 26, 27)
all_variants <- all_variants[ready_variants]

# Let's rename the groups back to the IUCN categories
iucn_names <- c("1" = "LC", "2" = "NT", "3" = "VU", "4" = "EN", "5" = "CR", 
                "6" = "EW", "7" = "EX", "8" = "DD")

# Process individual variants
all_plots <- list()
for (i in 1:length(all_variants)) {
  variant_name <- all_variants[[i]]@name
  message("Working with ", variant_name)
  # Set groupnames
  groupnames(all_variants[[i]]) <- iucn_names
  # Plot group curves
  group_curves <- curves(results(all_variants[[i]]), groups = TRUE)
  p <- plot(group_curves, main = variant_name, invert.x = TRUE)
  all_plots[[i]] <- p
  print(p)
}

library(gridExtra)
marrangeGrob(all_plots, nrow = 6, ncol = 2)
