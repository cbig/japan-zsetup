# Figure 4: Rank priority maps for the different taxa

library(dplyr)
library(grid)
library(lazyeval)
library(magick)
library(maptools)
library(purrr)
library(raster)
library(RColorBrewer)
library(raster)
library(rgdal)
library(sp)
library(tmap)

source("R/00_lib/utils.R")

# Load the project, creates japan_zproject. WARNING: loaded object is cached
# using memoise-package. If results have change, don't load the cached version.
# Also if the cache needs to updated, you will have to do this 
# > key <- list("zsetup")
# > saveCache(zproject_japan, key = key)
zproject_japan <- .load_zproject("zsetup", cache = TRUE, debug = TRUE)

# Get boundaries for Japan
sp_japan <- raster::getData('GADM', country = 'JPN', level = 0)

# Set parameters ----------------------------------------------------------

inner_margins <- c(0.02, 0.02, 0.02, 0)
bbox_japan <- matrix(c(-1409425, 2649200, 
                       1135015, 5132704),
                     nrow = 2, ncol = 2,
                     dimnames = list(c("x", "y"), 
                                     c("min", "max")))

# Define raster color scheme
breaks <- c(0, 0.25, 0.50, 0.83, 0.91, 0.98, 1)
colors <- rev(RColorBrewer::brewer.pal(length(breaks) - 1, "RdYlBu"))
labels <- (100 - breaks * 100)
labels <- cbind(labels[1:(length(labels) - 1)], labels[2:length(labels)])
labels[,2] <- paste(labels[,2], "%")
labels[6,2] <- ""
labels <- apply(labels, 1, paste, collapse = " - ")
labels[6] <- gsub(" - ", " %", labels[6])

# Manually add "lowest" and "highest" to the legend
labels[1] <- paste(labels[1], "(lowest)")
labels[6] <- paste(labels[6], "(highest)")

# Make plots --------------------------------------------------------------

construct_legend <- function(rank_raster) {
  rank_map_legend <- tmap::tm_shape(rank_raster, bbox = bbox_japan, is.master = TRUE) +
    tmap::tm_raster(title = "Best X% of \nthe solution", palette = colors,
                    labels = labels, breaks = breaks,
                    auto.palette.mapping = FALSE,
                    legend.show = TRUE) +
    tmap::tm_style("white", legend.only = TRUE, legend.position = c("left", "center"))
  return(rank_map_legend)
}

construct_map <- function(rank_raster, title, legend.show = TRUE) {
  rank_map <- tmap::tm_shape(rank_raster, bbox = bbox_japan, is.master = TRUE) +
    tmap::tm_raster(title = "Best X% of \nthe solution", 
                    palette = colors, labels = labels,
                    breaks = breaks, auto.palette.mapping = FALSE,
                    legend.show = legend.show) +
    tmap::tm_shape(sp_japan, bbox = bbox_japan) +
    tmap::tm_borders(col = "black", lwd = 0.05) +
    tmap::tm_style("white", title = title, title.size = 1.2, 
                   inner.margins = inner_margins, frame = FALSE)
  return(rank_map)
}

# Read in pixel-based rank data -------------------------------------------

# Following variants are used:
# 03_amp_caz_wgt
# 04_amp_abf_wgt
# 13_bir_caz_wgt
# 14_bir_abf_wgt
# 23_frf_caz_wgt
# 24_frf_abf_wgt
# 33_mam_caz_wgt
# 34_mam_abf_wgt
# 43_pla_caz_wgt
# 44_pla_abf_wgt
# 53_rep_caz_wgt
# 54_rep_abf_wgt

process_batch <- function(variant_ids, variant_labels) {
  
  if (length(variant_ids) != 6) {
    stop("6 variant IDs needed")
  } 
  
  # Define CRS
  crs_utm54n <- sp::CRS("+init=epsg:32654")
  
  variant_rank_rasters <- variant_ids %>% 
    purrr::map(zonator::get_variant, x = zproject_japan) %>% 
    purrr::map(zonator::rank_raster) %>% 
    purrr::map(raster::projectRaster, crs = crs_utm54n)
  
  variant_rank_maps <- purrr::map2(variant_rank_rasters,
                                   variant_labels, construct_map, legend.show = FALSE)

  # Use a single rank map to construct the legend
  rank_map_legend <- construct_legend(variant_rank_rasters[[1]])
  
  return(list(maps = variant_rank_maps, legend = rank_map_legend))
}

save_rank_maps <- function(rank_map_list, main_file_name, 
                           legend_file_name = NULL) {
  
  rank_map_objects <- rank_map_list[["maps"]]
  if (!is.null(legend_file_name)) {
    legend_object <- rank_map_list[["legend"]]
    tmap_save(legend_object, legend_file_name, 
              width = 900, height = 900, dpi = 600)
    # Crop legend
    img_legend <- magick::image_read(legend_file_name)
    img_legend <- magick::image_crop(img_legend, geometry = "900x800+20+80")
    magick::image_write(img_legend, path = legend_file_name)
    message("Saved cropped legend")
  }
  
  message("Saving main composite image...")
  
  # Where to save the sub figures
  component_files <- paste0("figs/figure_04/", "figure_04_0", 1:6, ".tiff")
  h_map <- 1600
  w_map <- 1600
  dpi <- 600
  
  res <- purrr::map2(rank_map_objects, component_files, 
                     tmap_save, width = w_map, height = h_map, dpi = dpi)
  # Read images back in and stack
  sub_figs_stack_1 <- magick::image_read(component_files[c(1, 3, 5)]) %>% 
    magick::image_append(stack = TRUE)
  sub_figs_stack_2 <- magick::image_read(component_files[c(2, 4, 6)]) %>% 
    magick::image_append(stack = TRUE)
  # Compose both stacks
  fig <- magick::image_append(c(sub_figs_stack_1, sub_figs_stack_2))
  # Write image
  magick::image_write(fig, path = main_file_name)
     
  message("Saved composite rank maps")
}

variant_ids_abf <- c(4, 14, 24, 34, 44, 54)
variant_ids_caz <- c(3, 13, 23, 33, 43, 53)
variant_labels <- c("Amphibians", "Birds", "Freshwater fish", "Mammals",
                    "Plants","Reptiles")

rank_maps_abf <- process_batch(variant_ids_abf, variant_labels)
#rank_maps_caz <- process_batch(variant_ids_caz, variant_labels)

# Save plots --------------------------------------------------------------

save_rank_maps(rank_maps_abf, 
               main_file_name = "figs/figure_04_main_abf.tiff", 
               legend_file_name = "figs/figure_04_legend.tiff")
#save_rank_maps(rank_maps_caz, 
#               main_file_name = "figs/figure_04_main_caz.png")
