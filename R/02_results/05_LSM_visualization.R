library(RColorBrewer)
library(htmlwidgets)
library(leaflet)
library(raster)
library(rgdal)
library(zonator)

source("R/00_lib/utils.R")

# Load the project, creates japan_zproject. WARNING: loaded object is cached
# using memoise-package. If results have change, don't load the cached version.
zproject_japan <- .load_zproject("zsetup", cache = TRUE, debug = TRUE)

output_dir <- "zsetup/taxa_all/34_caz_wgt_con/34_caz_wgt_con_out"
dsn_geojson <- file.path(output_dir, "34_caz_wgt_con_PPA_LSM_1.geojson")
dsn_rank_tif <- file.path(output_dir, "34_caz_wgt_con.CAZ_DEA.rank.compressed.tif")
variant34_LSM <- readOGR(dsn_geojson, layer = ogrListLayers(dsn_geojson))
variant34_rank <- raster(dsn_rank_tif)

zpal <- colorBin(
  palette = zlegend("spectral")$colors,
  domain = c(0, 1)
)

prefecture_popup <- paste0("<strong>Prefecture: </strong>", 
                           variant34_LSM$prefecture, 
                           "<br><strong>Mean rank: </strong>", 
                           variant34_LSM$Mean_rank,
                           "<br><strong>Spp dist. sum: </strong>", 
                           variant34_LSM$Spp_distribution_sum,
                           "<br><strong>Spp dist. > 10%:    </strong>",
                           variant34_LSM$Plus_10,
                           "<br><strong>Spp dist. > 1%: </strong>",
                           variant34_LSM$Plus_1,
                           "<br><strong>Spp dist. > 0.1%: </strong>",
                           variant34_LSM$Plus_01)

base_map <- leaflet(variant34_LSM) %>% addTiles()

poly_map <- base_map %>%         
              addPolygons(
                stroke = FALSE, fillOpacity = 0.8, smoothFactor = 0.5,
                color = ~zpal(Mean_rank), popup = prefecture_popup
              ) %>% 
            addLegend("bottomright", pal = zpal, values = c(0, 1),
                      title = "Mean rank",
                      #labFormat = labelFormat(prefix = "$"),
                      opacity = 1
            )

html_output <- gsub(".geojson", ".html", normalizePath(dsn_geojson))
saveWidget(poly_map, file = html_output)

raster_map <- base_map %>% 
  addRasterImage(variant34_rank, colors = zpal, opacity = 0.8) %>%
  addLegend(pal = zpal, values = c(0, 1),
            title = "Rank priority")
