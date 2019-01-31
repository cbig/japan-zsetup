# Figure 6: Priority areas maps without and with
# the existing protected areas as mask layers
library(dplyr)
library(grid)
library(lazyeval)
library(magick)
library(maptools)
library(raster)
library(RColorBrewer)
library(rgdal)
library(sp)
library(tmap)

source("R/00_lib/utils.R")

# Helper functions --------------------------------------------------------

get_map_params <- function(type, legend.reversed = FALSE) {
  
  if (type == "full") {  
    breaks <- c(0, 0.25, 0.50, 0.83, 0.91, 0.98, 1)
    colors <- rev(RColorBrewer::brewer.pal(length(breaks) - 1, "RdYlBu"))
    labels <- (100 - breaks * 100)
    labels <- cbind(labels[1:(length(labels) - 1)], labels[2:length(labels)])
    labels[,2] <- paste(labels[,2], "%")
    labels[6,2] <- ""
    labels <- apply(labels, 1, paste, collapse = " - ")
    labels[6] <- gsub(" - ", " %", labels[6])
  } else if (type == "3level") {
    # Define raster color scheme
    breaks <- c(0, 0.83, 0.91, 1)
    #colors <- c("lightgrey", "yellow", "blue", "red")
    colors <- c("lightgrey", rev(RColorBrewer::brewer.pal(6, "RdYlBu")[c(1, 6)]))
    labels <- (100 - breaks[2:(length(breaks) - 1)] * 100)
    labels <- c("Non-protected", "Aichi target (top 17%)", 
                                 "High/med rank PA (top 9%)")
  } else if (type == "4level") {
    # Define raster color scheme
    breaks <- c(0, 0.83, 0.91, 0.98, 1)
    #colors <- c("lightgrey", "yellow", "blue", "red")
    colors <- c("lightgrey", RColorBrewer::brewer.pal(6, "RdYlBu")[c(1,3,6)])
    labels <- (100 - breaks[2:(length(breaks) - 1)] * 100)
    labels <- c("Non-protected", "Aichi target (top 17%)", 
                "Med rank PA (top 9%)", "High rank PA (top 2%)")
  } else {
    stop("Unknown type")
  }
  
  params <- list()
  
  params$inner_margins <- c(0.02, 0.02, 0.02, 0)
  params$bbox <- matrix(c(-1409425, 2649200, 
                          1135015, 5132704),
                        nrow = 2, ncol = 2,
                        dimnames = list(c("x", "y"), 
                                        c("min", "max")))
  
  if (legend.reversed) {
    params$breaks <- rev(breaks)
    params$colors <- rev(colors)
    params$labels <- rev(labels)
  } else {
    params$breaks <- breaks
    params$colors <- colors
    params$labels <- labels
  }
  return(params)
}

create_raster_levels <- function(raster, params) {
  # Create a RasterLayer with a RAT
  rat_raster <- raster::ratify(raster)
  rat <- levels(rat_raster)[[1]]
  rat$priorities_cat <- cut(rat$ID, breaks = params$breaks)
  rat$priorities_cat <- factor(rat$priorities_cat,
                               levels = levels(rat$priorities_cat))
  levels(rat_raster) <- rat
  return(rat_raster)
}

construct_map <- function(raster, title, type, legend.show = TRUE) {
  
  params <- get_map_params(type = type, legend.reversed = FALSE)
  rank_map <- tmap::tm_shape(raster, bbox = params$bbox, is.master = TRUE) +
    tmap::tm_raster(title = "", 
                    palette = params$colors, labels = params$labels,
                    breaks = params$breaks, auto.palette.mapping = FALSE,
                    legend.show = TRUE) +
    tmap::tm_legend(legend.title.size = 3) +
    tmap::tm_style("white", title = title, title.size = 1, 
                   inner.margins = params$inner_margins, frame = FALSE)
  return(rank_map)
}

# Load the project, creates japan_zproject. WARNING: loaded object is cached
# using memoise-package. If results have change, don't load the cached version.
# Also if the cache needs to updated, you will have to do this manually.
# > key <- list("zsetup")
# > saveCache(zproject_japan, key = key)
zproject_japan <- .load_zproject("zsetup", cache = TRUE, debug = TRUE)

# Get boundaries for Japan
sp_japan <- raster::getData('GADM', country = 'JPN', level = 0)

# PA data -----------------------------------------------------------------

pa_raster <- raster::raster("../Data.150928/environment/all_PA.tif")
pa_raster[pa_raster == 2] <- 1
# Define CRS (WGS84), http://spatialreference.org/ref/epsg/wgs-84/
#crs_wgs84 <- sp::CRS("+init=epsg:4326")
crs_utm54n <- sp::CRS("+init=epsg:32654")
#crs_nsdic <- sp::CRS("+init=epsg:6933")

# Read in pixel-based rank data -------------------------------------------

# Following variants are used:
# 66_abf_wgt_con
# 70_abf_wgt_con_hm3

nopa_variant <- zonator::get_variant(zproject_japan, 66)
nopa_rank <- zonator::rank_raster(nopa_variant)
nopa_rank <- raster::projectRaster(nopa_rank, crs = crs_utm54n)

pa_variant <- zonator::get_variant(zproject_japan, 70)
pa_rank <- zonator::rank_raster(pa_variant)
pa_rank <- raster::projectRaster(pa_rank, crs = crs_utm54n)

# Make plots --------------------------------------------------------------

maps <- list()

maps[[1]] <- construct_map(pa_rank, title = "(a) with current PAN", 
                             type = "full", legend.show = TRUE)
maps[[2]] <- construct_map(pa_rank, title = "(b) with current PAN", 
                            type = "3level", legend.show = TRUE)

maps[[3]] <- construct_map(nopa_rank, title = "(c) without current PAN", 
                               type = "full", legend.show = TRUE)
maps[[4]] <- construct_map(nopa_rank, title = "(d) without current PAN", 
                              type = "3level", legend.show = TRUE)

# Save plots --------------------------------------------------------------

component_files <- paste0("figs/figure_06/", "figure_06_0", 1:4, ".png")
h_map <- 2000
w_map <- 2000
dpi <- 600

res <- purrr::map2(maps, component_files, 
                   tmap_save, width = w_map, height = h_map, dpi = dpi)
# Read images back in and stack
sub_figs_stack_1 <- magick::image_read(component_files[c(1, 3)]) %>% 
  magick::image_append(stack = TRUE)
sub_figs_stack_2 <- magick::image_read(component_files[c(2, 4)]) %>% 
  magick::image_append(stack = TRUE)
# Compose both stacks
fig <- magick::image_append(c(sub_figs_stack_1, sub_figs_stack_2))
# Write image
magick::image_write(fig, path = "figs/figure_06.png")
