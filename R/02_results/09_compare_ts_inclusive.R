library(ggplot2)
library(hrbrthemes)
library(raster)
library(rasterVis)
library(tidyverse)
library(zonator)

source("R/00_lib/utils.R")

# Helper functions --------------------------------------------------------

calc_common_top_fraction <- function(rasters, type, top_fraction) {
  
  if (!is(rasters, "RasterBrick")) {
    stop("'rasters' object not a RasterBrick")
  }
  
  
  if (type == "intersect") {
    message("Calculating fraction of the landscape needed to cover the top ", 
            top_fraction * 100, "% in all taxon-secific priority ranks (intersect)...")  
  } else if (type == "union") {
    message("Calculating fraction of the landscape needed to cover the top ", 
            top_fraction * 100, "% with occurrence from at least one taxon-secific ",
            "priority rank (union)...")  
  } else {
    stop("Unknown type: ", type)
  }
  
  # Total number of layers
  n_layers <- raster::nlayers(rasters)
  # Total area
  area_mask <- sum(!is.na(rasters))
  total_area <- sum(raster::getValues(area_mask > 0))
  
  # Get individual top fractions
  top_fracs <- rasters >= (1.0 - top_fraction)
  # Sum the brickunio
  agg_top_fracs <- sum(top_fracs)
  # agg_area_mask <- sum(!is.na(top_fracs))
  # agg_total_area <- sum(raster::getValues(agg_area_mask > 0))
  # Calculate the fraction of cells which are in top fraction of 
  # A) all layers (type = "intersect")
  # B) any of the layers (type = "union")
  if (type == "intersect") {
    agg_top_area <- sum(raster::getValues(agg_top_fracs == n_layers), 
                        na.rm = TRUE)
  } else if (type == "union") {
    agg_top_area <- sum(raster::getValues(agg_top_fracs > 0), 
                        na.rm = TRUE)
  }
  if (agg_top_area == 0) {
    stop("Layers have no common cells in the top fraction ", top_fraction)
  } else {
    # Return the fraction of common cells in the top fraction
    return(data.frame(target_fraction = top_fraction,
                      area_fraction_covered = agg_top_area / total_area))
  }
}

create_mask <- function(rasters, top_fraction, type) {
  if (!is(rasters, "RasterBrick")) {
    stop("'rasters' object not a RasterBrick")
  }
  
  if (type %in% c("intersect", "union")) {
    message("Creating ", type, " mask for ", top_fraction * 100, "% top fraction...") 
  } else {
    stop("Unknown type: ", type)
  }
  
  # Get individual top fractions
  top_fracs <- rasters >= (1.0 - top_fraction)
  # Sum the brickunio
  agg_top_fracs <- sum(top_fracs)
 
  if (type == "intersect") {
    res_raster <- agg_top_fracs == n_layers
  } else if (type == "union") {
    res_raster <- agg_top_fracs > 0
  }
  
  # Reclassify
  m <- c(0, 1, 1,  2)
  rclmat <- matrix(m, ncol = 2, byrow = TRUE)
  res_raster <- raster::reclassify(res_raster, rclmat)
  
  return(res_raster)
}

plot_fractions <- function(rasters, top_fraction, labels = NULL, ...) {
  top_fracs <- rasters >= (1.0 - top_fraction)
  
  process_layer <- function(x) {
    rat <- data.frame(ID = c(0, 1),
                      status = c("Not Selected", "Selected"))
    xx <- raster::ratify(x)
    suppressWarnings({levels(xx) <- rat})
    return(xx)
  }
  
  top_fracs <- raster::stack(lapply(1:raster::nlayers(top_fracs), 
                                    function(x) {
                                      return(process_layer(top_fracs[[x]]))
                                    }))
  
  colours <- c("grey40", "#4DAF4A")
  
  if (!is.null(labels)) {
    names(top_fracs) <- labels
  }
  
  p <- rasterVis::levelplot(top_fracs, scales = list(draw = FALSE),
                            col.regions = colours, xlab = "", ylab = "", 
                            colorkey = list(space = "bottom", height = 1), ...)
  return(p)
}

# Read in data ------------------------------------------------------------

# Load the project, creates japan_zproject. WARNING: loaded object is cached
# using memoise-package. If results have change, don't load the cached version.
# Also if the cache needs to updated, you will have to do this manually.
zproject_japan <- .load_zproject("zsetup", cache = TRUE, debug = TRUE)

# Read in pixel-based rank data -------------------------------------------

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
# 
# Inclusive (joint):
# 
# 64_abf_wgt
# 
ts_variant_ids <- c(4, 14, 24, 34, 44, 54)
incl_variant_id <- 64

ts_rank_rasters <- ts_variant_ids %>% 
  purrr::map(zonator::get_variant, x = zproject_japan) %>% 
  purrr::map(zonator::rank_raster) %>% 
  raster::brick()

incl_rank_raster <- zonator::rank_raster(zonator::get_variant(zproject_japan, 
                                                              incl_variant_id))

# CHECK 1: intersection of top 17% fractions ------------------------------

# Calculate area fraction needed with 5% intervals
areas_needed <- purrr::map_df(seq(0.1, 0.9, 0.05), calc_common_top_fraction,
                              rasters = ts_rank_rasters, type = "intersect")

# Hone in to a more specific range
detailed_areas_needed <- purrr::map_df(seq(0.63, 0.635, 0.001), 
                                       calc_common_top_fraction,
                                       rasters = ts_rank_rasters, 
                                       type = "intersect")

# CHECK 2: union of 17% ---------------------------------------------------

union_areas_needed <- purrr::map_df(seq(0.01, 0.1, 0.01), 
                                    calc_common_top_fraction,
                                    rasters = ts_rank_rasters,
                                    type = "union")

detailed_union_areas_needed <- purrr::map_df(seq(0.039, 0.04, 0.0002), 
                                             calc_common_top_fraction,
                                             rasters = ts_rank_rasters,
                                             type = "union")


# Create mask -------------------------------------------------------------

union_mask <- create_mask(ts_rank_rasters, 0.0394, type = "union")

# Test, this should be close to 17%
# sum(raster::getValues(union_mask == 2), na.rm = TRUE) / sum(raster::getValues(union_mask >= 1), na.rm = TRUE)

raster::writeRaster(union_mask, "../Data.150928/environment/ts_top_mask_00394.tif",
                    overwrite = TRUE)

# Plot figures ------------------------------------------------------------

p1 <- ggplot(areas_needed, aes(x = target_fraction, 
                                y = area_fraction_covered)) +
  geom_vline(xintercept = 0.633, linetype = "dashed", color = "darkgrey") +
  geom_hline(yintercept = 0.17, linetype = "dashed", color = "darkgrey") +
  geom_line() +
  scale_y_continuous(breaks = seq(0, 0.7, 0.1),
                     labels = paste0(seq(0, 0.7, 0.1) * 100)) +
  scale_x_continuous(breaks = seq(0, 1.0, 0.1),
                     labels = paste0(seq(0, 1.0, 0.1) * 100)) +
  ylab("Taxon-specific top fraction (%)") +
  xlab("Land area needed (%)\n") +
  theme_ipsum_rc()

p2 <- ggplot(union_areas_needed, aes(x = target_fraction, 
                               y = area_fraction_covered)) +
  geom_vline(xintercept = 0.0394, linetype = "dashed", color = "darkgrey") +
  geom_hline(yintercept = 0.17, linetype = "dashed", color = "darkgrey") +
  geom_line() +
  scale_y_continuous(breaks = seq(0, 0.4, 0.1),
                     labels = paste0(seq(0, 0.4, 0.1) * 100)) +
  scale_x_continuous(breaks = seq(0, 0.1, 0.01),
                     labels = paste0(seq(0, 0.1, 0.01) * 100)) +
  xlab("Taxon-specific top fraction (%)") +
  ylab("Land area needed (%)\n") +
  theme_ipsum_rc()

ggsave("figs/ts_land_area_intersect.png", p1, width = 7, height = 7)
ggsave("figs/ts_land_area_union.png", p2, width = 7, height = 7)

# Rasters

layer_names <- c("Amphibians", "Birds", "Freshwater fish", "Mammals", 
                 "Plants", "Reptiles")

p3 <- plot_fractions(ts_rank_rasters, 0.0394, labels = layer_names, 
                     main = "Top 3.94%")

jaccards <- zonator::cross_jaccard(ts_rank_rasters, c(0.9606))
colnames(jaccards[[1]]) <- layer_names
rownames(jaccards[[1]]) <- layer_names

png("figs/ts_top_00394_fractions.png", width = 1200, height = 800)
p3
dev.off()
