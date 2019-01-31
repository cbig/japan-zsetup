# Figure XX: Hex-grid priority pattern maps
library(dplyr)
library(tidyr)
library(sp)
library(raster)
library(rgeos)
library(rgbif)
library(viridis)
library(gridExtra)
library(rasterVis)
set.seed(1)

study_area <- getData("GADM", country = "LKA", level = 0, 
                      path = "tmp/hexagonal-grids/") %>% 
  disaggregate %>% 
  geometry

study_area <- sapply(study_area@polygons, slot, "area") %>% 
{which(. == max(.))} %>% 
  study_area[.]

plot(study_area, col = "grey50", bg = "light blue", axes = TRUE)
text(81.5, 9.5, "Study Area:\nSri Lanka")

size <- 0.5
hex_points <- spsample(study_area, type = "hexagonal", cellsize = size)
hex_grid <- HexPoints2SpatialPolygons(hex_points, dx = size)
plot(study_area, col = "grey50", bg = "light blue", axes = TRUE)
plot(hex_points, col = "black", pch = 20, cex = 0.5, add = T)
plot(hex_grid, border = "orange", add = T)

make_grid <- function(x, type, cell_width, cell_area, clip = FALSE) {
  if (!type %in% c("square", "hexagonal")) {
    stop("Type must be either 'square' or 'hexagonal'")
  }
  
  if (missing(cell_width)) {
    if (missing(cell_area)) {
      stop("Must provide cell_width or cell_area")
    } else {
      if (type == "square") {
        cell_width <- sqrt(cell_area)
      } else if (type == "hexagonal") {
        cell_width <- sqrt(2 * cell_area / sqrt(3))
      }
    }
  }
  # buffered extent of study area to define cells over
  ext <- as(extent(x) + cell_width, "SpatialPolygons")
  projection(ext) <- projection(x)
  # generate grid
  if (type == "square") {
    g <- raster(ext, resolution = cell_width)
    g <- as(g, "SpatialPolygons")
  } else if (type == "hexagonal") {
    # generate array of hexagon centers
    g <- spsample(ext, type = "hexagonal", cellsize = cell_width, offset = c(0, 0))
    # convert center points to hexagons
    g <- HexPoints2SpatialPolygons(g, dx = cell_width)
  }
  
  # clip to boundary of study area
  if (clip) {
    g <- gIntersection(g, x, byid = TRUE)
  } else {
    g <- g[x, ]
  }
  # clean up feature IDs
  row.names(g) <- as.character(1:length(g))
  return(g)
}

study_area_utm <- CRS("+proj=utm +zone=44 +datum=WGS84 +units=km +no_defs") %>% 
  spTransform(study_area, .)
# without clipping
hex_grid <- make_grid(study_area_utm, cell_area = 625, clip = FALSE)
plot(study_area_utm, col = "grey50", bg = "light blue", axes = FALSE)
plot(hex_grid, border = "orange", add = TRUE)
box()
# with clipping
hex_grid <- make_grid(study_area_utm, cell_area = 625, clip = TRUE)
plot(study_area_utm, col = "grey50", bg = "light blue", axes = FALSE)
plot(hex_grid, border = "orange", add = TRUE)
box()

japan <- getData(name = "GADM", country = "JPN", level = 0, 
                   path = "tmp/hexagonal-grids/") %>% 
  disaggregate %>% 
  geometry

japan_utm <- CRS("+proj=utm +zone=54 +datum=WGS84 +units=km +no_defs") %>% 
  spTransform(japan, .)

hex_jpn <- make_grid(japan_utm, type = "hexagonal", cell_width = 100, clip = FALSE)
plot(japan_utm, col = "grey50", bg = "light blue", axes = FALSE)
plot(hex_jpn, border = "orange", add = TRUE)
box()

amp_rank <- raster::raster("zsetup/amphibians/05_amp_caz_wgt_con/05_amp_caz_wgt_con_out/05_amp_caz_wgt_con.CAZ_DE.rank.compressed.tif") %>% 
  projectRaster(t_crop, to = raster(hex_jpn, res = 1)) %>% 
  setNames('rank')
hex_amp_rank <- extract(amp_rank, hex_jpn, fun = mean, na.rm = TRUE, sp = TRUE)

breaks <- c(0, 0.25, 0.50, 0.83, 0.91, 0.98, 1)
colors <- rev(RColorBrewer::brewer.pal(length(breaks) - 1, "RdYlBu"))
labels <- (100 - breaks * 100)
labels <- cbind(labels[1:(length(labels) - 1)], labels[2:length(labels)])
labels[,2] <- paste(labels[,2], "%")
labels[6,2] <- ""
labels <- apply(labels, 1, paste, collapse = " - ")
labels[6] <- gsub(" - ", " %", labels[6])

p2 <- spplot(hex_amp_rank,
             col.regions = colors,
             at = breaks,
             colorkey = list(
               labels = list(at = breaks, 
                             labels = labels)
             )
)
plot(japan_utm, col = "grey50", bg = "light blue", axes = FALSE)
plot(p2)
