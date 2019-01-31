library(ggplot2)
library(gridExtra)
library(raster)
library(rasterVis)
library(zonator)

source("R/00_lib/utils.R")

# Load the project, creates japan_zproject. WARNING: loaded object is cached
# using memoise-package. If results have change, don't load the cached version.
zproject_japan <- .load_zproject("zsetup", cache = TRUE, debug = TRUE)

# Read in the PA mask files
low_pa <- raster("../Data.150928/environment/low_PA.tif")
low_pa[low_pa == 0] <- NA
low_pa[low_pa == -1] <- NA
med_pa <- raster("../Data.150928/environment/medium_PA.tif")
med_pa[med_pa == 0] <- NA
med_pa[med_pa == -1] <- NA
high_pa <- raster("../Data.150928/environment/high_PA.tif")
high_pa[high_pa == 0] <- NA
high_pa[high_pa == -1] <- NA


plot_priority_in_pa <- function(variant_no) {
  
  message("Working with variant: ", variant_no, "...")
  
  # Read in results
  variant_results <- results(get_variant(zproject_japan, variant_no))
  variant_rank <- rank_raster(variant_results)
  
  p1 <- plot_hist(variant_rank, high_pa, add.median = TRUE, 
                  add.mean = FALSE, binwidth = 0.02, 
                  title = paste("Variant", variant_no, "for high PAs"))
  
  p2 <- plot_hist(variant_rank, med_pa, add.median = TRUE, 
                  add.mean = FALSE, binwidth = 0.02, 
                  title = paste("Variant", variant_no, "for medium PAs"))
  
  p3 <- plot_hist(variant_rank, low_pa, add.median = TRUE, 
                  add.mean = FALSE, binwidth = 0.02, 
                  title = paste("Variant", variant_no, "for low PAs"))
  
  return(grid.arrange(p1, p2, p3, nrow = 1, ncol = 3))
}

for (i in c(4, 9, 14, 19, 24, 29, 34)) {
  out_file <- file.path("figs", paste0("variant_", i, "_PA_rank_histogram.png"))
  ggsave(out_file, plot_priority_in_pa(i), width = 20, height = 7.5, 
         units = "cm")
}
