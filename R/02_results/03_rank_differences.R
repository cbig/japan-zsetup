library(raster)
library(zonator)

japan_project <- create_zproject(root = "zsetup/", debug = TRUE)

result_rasters <- rank_rasters(japan_project)
names(result_rasters) <- gsub("do_", "", names(result_rasters))

# Substractions -----------------------------------------------------------

output_dir <- "zsetup/output/substractions"
dir.create(output_dir)

writeRaster(result_rasters[["X02_caz_wgt"]] - result_rasters[["X01_caz"]],
            filename = file.path(output_dir, "substraction_01_from_02.tif"))

writeRaster(result_rasters[["X03_caz_wgt_cond"]] - result_rasters[["X02_caz_wgt"]] ,
            filename = file.path(output_dir, "substraction_02_from_03.tif"))

writeRaster(result_rasters[["X04_abf"]] - result_rasters[["X01_caz"]],
            filename = file.path(output_dir, "substraction_01_from_04.tif"))

writeRaster(result_rasters[["X06_caz_wgt_msk"]] - result_rasters[["X02_caz_wgt"]] ,
            filename = file.path(output_dir, "substraction_02_from_06.tif"))
