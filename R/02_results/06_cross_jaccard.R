library(raster)
library(zonator)

source("R/00_lib/utils.R")

# Helper functions --------------------------------------------------------

calculate_jaccards <- function(rank_stack, x.min, x.max, y.min, y.max) {
  
  variant_names <- c("Amhibians", "Birds", "Freswater fish", "Mammals", 
                     "Plants", "Reptiles", "All")
  
  jaccards <- matrix(nrow = nlayers(rank_stack), ncol = nlayers(rank_stack))
  for (i in 1:nrow(jaccards)) {
    for (j in 1:ncol(jaccards)) {
      if (i == j) {
        jaccards[i, j] <- 1
      }
      else {
        if (is.na(jaccards[j, i])) {
          message(paste0("Calculating Jaccard index between ", 
                         names(rank_stack[[i]]), 
                         " and ", names(rank_stack[[j]])))
          jaccards[i, j] <- jaccard(rank_stack[[i]], rank_stack[[j]], 
                                    x.min = x.min, x.max = x.max, 
                                    y.min = y.min, y.max = y.max)
        }
        else {
          jaccards[i, j] <- NA
        }
      }
    }
  }
  jaccards <- as.data.frame(jaccards)
  colnames(jaccards) <- variant_names
  rownames(jaccards) <- variant_names
  return(jaccards)
}

# Calculate ---------------------------------------------------------------

# Load the project, creates japan_zproject. WARNING: loaded object is cached
# using memoise-package. If results have change, don't load the cached version.
# > key <- list("zsetup")
# > saveCache(zproject_japan, key = key)
zproject_japan <- .load_zproject("zsetup", cache = TRUE, debug = TRUE)

# Get rank rasters for selected variants and create a RasterStack object. NOTE:
# here we don't check for the existence of results and the following will fail
# if any of the target rasters are missing.
# NOTE: the NoPAN variants do not include condition (HII) becuause 
# Fig 2 doesn't either.
variants_caz_nopan <- c(3, 13, 23, 33, 43, 53, 63)
variants_caz_pan <- c(9, 19, 29, 39, 49, 59, 69)
variants_abf_nopan <- c(4, 14, 24, 34, 44, 54, 64)
variants_abf_pan <- c(10, 20, 30, 40, 50, 60, 70)

# Transform the list of rasters into a RasterStack. NOTE: we are using 
# quick = TRUE beacause we know all priority rank rasters have the same 
# dimensions.
rank_stack_caz_nopan <- stack(lapply(variants_caz_nopan, 
                               function(x) rank_raster(get_variant(zproject_japan, x))), 
                        quick = TRUE)
rank_stack_caz_pan <- stack(lapply(variants_caz_pan, 
                                     function(x) rank_raster(get_variant(zproject_japan, x))), 
                              quick = TRUE)
rank_stack_abf_nopan <- stack(lapply(variants_abf_nopan, 
                               function(x) rank_raster(get_variant(zproject_japan, x))), 
                        quick = TRUE)
rank_stack_abf_pan <- stack(lapply(variants_abf_pan, 
                                     function(x) rank_raster(get_variant(zproject_japan, x))), 
                              quick = TRUE)

# Top 17% (no PAN)
jaccards_caz_nopan <- calculate_jaccards(rank_stack_caz_nopan, 
                                         x.min = 0.83, x.max = 1.0, 
                                         y.min = 0.83, y.max = 1.0)
jaccards_abf_nopan <- calculate_jaccards(rank_stack_abf_nopan, 
                                         x.min = 0.83, x.max = 1.0, 
                                         y.min = 0.83, y.max = 1.0)

# Top expansion areas (17%-9%, PAN included)
jaccards_caz_pan <- calculate_jaccards(rank_stack_caz_pan, 
                                       x.min = 0.83, x.max = 0.91, 
                                       y.min = 0.83, y.max = 0.91)
jaccards_abf_pan <- calculate_jaccards(rank_stack_abf_pan, 
                                       x.min = 0.83, x.max = 0.91, 
                                       y.min = 0.83, y.max = 0.91) 

write.csv(t(jaccards_abf_nopan), file = "tables/jaccards_abf_nopan.csv", na = "-", 
          row.names = TRUE)
write.csv(t(jaccards_abf_pan), file = "tables/jaccards_abf_pan.csv", na = "-", 
          row.names = TRUE)
write.csv(t(jaccards_caz_nopan), file = "tables/jaccards_caz_nopan.csv", na = "-", 
          row.names = TRUE)
write.csv(t(jaccards_caz_pan), file = "tables/jaccards_caz_pan.csv", na = "-", 
          row.names = TRUE)
