# NOTE: This script is used to generate Zonation setups for public release.
# See XXXXXX for a more detailed explanation.
#
# NOTE: you will need the latest version for this to work
# zonator > 0.5.0
# devtools::install_github("cbig/zonator")
library(dplyr)
library(raster)
library(readr)
library(readxl)
library(zonator)

source("R/00_lib/utils.R")

# Generate variants for all taxa ------------------------------------------

# Define names for variants. "[ID]" is a placeholder for running id, "[TX]" is
# for taxon codes.
variant_templates <- c("[TX]_caz",
                       "[TX]_abf",
                       "[TX]_caz_wgt",
                       "[TX]_abf_wgt",
                       "[TX]_caz_wgt_con",
                       "[TX]_abf_wgt_con",
                       "[TX]_caz_wgt_con_hm3",
                       "[TX]_abf_wgt_con_hm3")

zsetup_root <- "zsetup_release"

# Make a list to hold information on the different taxa names, short codes
# and Excel sheet names (sheet names not actually needed).
spp_data <- list("amphibians" = list("code" = "amp", "sheet" = "amphibians"),
                 "birds" = list("code" = "bir", "sheet" = "birds"),
                 "freshwater_fish" = list("code" = "frf",
                                           "sheet" = "freshwater_fish"),
                 "mammals" = list("code" = "mam", "sheet" = "mammals"),
                 "plants" = list("code" = "pla", "sheet" = "plants"),
                 "reptiles" = list("code" = "rep", "sheet" = "reptiles"))

aux_data_file <- "../Zendo/input_data/biodiv_features.csv"
# Read in the weights file containting information on the family, AUC,
# Red-list status and WC for taxa in turn
#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !! w_field defines which weights are used
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
w_field <- "weight"

dat <- as.data.frame(readr::read_csv(aux_data_file))

# Read in aggregate weights
agg_weights_file <- "../Zendo/input_data/aggregate_weights.csv"
agg_weights <- readr::read_csv(agg_weights_file) %>%
  dplyr::filter(agg_category == "Aggregation") %>%
  dplyr::select(-agg_category) %>%
  tidyr::gather(taxon, agg_weight)

# Auxiliary data files

ppa_raster_file <- "../Zendo/input_data/prefecture.tif"
ppa_config_file <- "ppa_config.txt"

# FIXME: freswater fish may need a different condition file
condition_raster_file <- "../Zendo/input_data/hii_rescaled.tif"
condition_config_file <- "condition_config.txt"

# Hierarchical masks
hm3_raster_file <- "../Zendo/input_data/PA_3_levels.tif"

# Set up the project

japan_zproject <- initiate_zproject(zsetup_root = zsetup_root,
                                    variant_templates = variant_templates,
                                    spp_data = spp_data,
                                    data_dir = "../Zenodo/input_data",
                                    prefix_spp_paths = "../..")

# Set run configuration parameters ----------------------------------------

# In principle, the preprocessing of the same type of variant for each taxon
# (e.g. 01_amp_caz, 06_bir_caz, 11_frf_caz etc.) is the same. Thus  let's
# automate it a bit.

# Use a counting variable to identify variant IDs.
variant_id <- 1
for (taxon in spp_data) {

  variant <- get_variant(japan_zproject, variant_id)

  # Generate groups
  # Parse groups based on the Red-list status
  groups <- lookup_rl_group(dat, featurenames(variant))

  # [XX]_[TX]_caz -------------------------------------------------------------

  message("Editing variant: ", variant@name)
  # Set groups
  groups(variant) <- as.vector(groups)
  # Set groups use and groups file
  variant <- set_dat_param(variant, "use groups", 1)
  # Note that groups file location is always relative to the bat file
  groups_file <- file.path(variant@name, paste0(variant@name, "_groups.txt"))
  variant <- set_dat_param(variant, "groups file", groups_file)

  # Set post-processing (LSM). First, let's create the file itself (zonator
  # can't handle this yet). The file needs to be created only once per raxon
  # since all the variants can use the same file.
  ppa_file_name <- file.path(zsetup_root, taxon$sheet, ppa_config_file)
  ppa_cmd_string <- paste(c("LSM", ppa_raster_file, 0, -1, 0), collapse = " ")
  write(ppa_cmd_string, ppa_file_name)
  # Need to define ppa_config.txt relative to the bat-file (same dir)-
  variant <- set_dat_param(variant, "post-processing list file",
                           ppa_config_file)

  # Save variant
  save_zvariant(variant, dir = file.path(zsetup_root, taxon$sheet),
                overwrite = TRUE, debug_msg = FALSE)
  create_sh_file(variant)

  # [XX]_[TX]_abf -------------------------------------------------------------

  variant_id <- variant_id + 1
  variant <- get_variant(japan_zproject, variant_id)
  message("Editing variant: ", variant@name)

  # Manipulate the spp data
  variant_spp_data <- sppdata(variant)
  sppdata(variant) <- variant_spp_data

  # Set groups, use the same as with caz. NOTE: This must be done after setting
  # the spp data
  groups(variant) <- as.vector(groups)

  # Set cell removal rule
  variant <- set_dat_param(variant, "removal rule", 2)

  # Set groups use and groups file
  variant <- set_dat_param(variant, "use groups", 1)
  # Note that groups file location is always relative to the bat file
  groups_file <- file.path(variant@name, paste0(variant@name, "_groups.txt"))
  variant <- set_dat_param(variant, "groups file", groups_file)

  # Set PPA file
  variant <- set_dat_param(variant, "post-processing list file",
                           ppa_config_file)

  # Save variant
  save_zvariant(variant, dir = file.path(zsetup_root, taxon$sheet),
                overwrite = TRUE, debug_msg = FALSE)
  create_sh_file(variant)

  # [XX]_[TX]_caz_wgt ---------------------------------------------------------

  variant_id <- variant_id + 1
  variant <- get_variant(japan_zproject, variant_id)
  message("Editing variant: ", variant@name)

  # Set weights
  wgts <- lookup_weight(dat, featurenames(variant), w_field = w_field)
  # sppweights()<- not implemented yet.
  variant_spp_data <- sppdata(variant)
  variant_spp_data$weight <- wgts
  sppdata(variant) <- variant_spp_data

  # Set groups, use the same as with previous
  groups(variant) <- as.vector(groups)

  # Set cell removal rule
  variant <- set_dat_param(variant, "removal rule", 1)

  # Set groups use and groups file
  variant <- set_dat_param(variant, "use groups", 1)
  # Note that groups file location is always relative to the bat file
  groups_file <- file.path(variant@name, paste0(variant@name, "_groups.txt"))
  variant <- set_dat_param(variant, "groups file", groups_file)

  # Set PPA file
  variant <- set_dat_param(variant, "post-processing list file",
                           ppa_config_file)

  # Save variant
  save_zvariant(variant, dir = file.path(zsetup_root, taxon$sheet),
                overwrite = TRUE, debug_msg = FALSE)
  create_sh_file(variant)

  # [XX]_[TX]_abf_wgt ---------------------------------------------------------

  variant_id <- variant_id + 1
  variant <- get_variant(japan_zproject, variant_id)
  message("Editing variant: ", variant@name)

  # Set weights
  wgts <- lookup_weight(dat, featurenames(variant), w_field = w_field)
  # sppweights()<- not implemented yet.
  variant_spp_data <- sppdata(variant)
  variant_spp_data$weight <- wgts
  sppdata(variant) <- variant_spp_data

  # Set groups, use the same as with previous
  groups(variant) <- as.vector(groups)

  # Set cell removal rule
  variant <- set_dat_param(variant, "removal rule", 2)

  # Set groups use and groups file
  variant <- set_dat_param(variant, "use groups", 1)
  # Note that groups file location is always relative to the bat file
  groups_file <- file.path(variant@name, paste0(variant@name, "_groups.txt"))
  variant <- set_dat_param(variant, "groups file", groups_file)

  # Set PPA file
  variant <- set_dat_param(variant, "post-processing list file",
                           ppa_config_file)

  # Save variant
  save_zvariant(variant, dir = file.path(zsetup_root, taxon$sheet),
                overwrite = TRUE, debug_msg = FALSE)
  create_sh_file(variant)

  # [XX]_[TX]_caz_wgt_con -----------------------------------------------------

  variant_id <- variant_id + 1
  variant <- get_variant(japan_zproject, variant_id)
  message("Editing variant: ", variant@name)

  # Set weights
  wgts <- lookup_weight(dat, featurenames(variant), w_field = w_field)
  # sppweights()<- not implemented yet.
  variant_spp_data <- sppdata(variant)
  variant_spp_data$weight <- wgts
  sppdata(variant) <- variant_spp_data

  # Set groups, use the same as with previous
  groups(variant) <- as.vector(groups)

  # Set cell removal rule
  variant <- set_dat_param(variant, "removal rule", 1)

  # Set groups use and groups file
  variant <- set_dat_param(variant, "use groups", 1)
  # Note that groups file location is always relative to the bat file
  groups_file <- file.path(variant@name, paste0(variant@name, "_groups.txt"))
  variant <- set_dat_param(variant, "groups file", groups_file)

  # Set PPA file
  variant <- set_dat_param(variant, "post-processing list file",
                           ppa_config_file)

  # Set condition. First, let's create the file itself (zonator
  # can't handle this yet). The file needs to be created only once per raxon
  # since all the variants can use the same file.
  condition_file_name <- file.path(zsetup_root, taxon$sheet,
                                   condition_config_file)
  condition_cmd_string <- paste(c(1, condition_raster_file), collapse = " ")
  write(condition_cmd_string, condition_file_name)
  # Need to define ppa_config.txt relative to the bat-file (same dir)-
  variant <- set_dat_param(variant, "use condition layer", 1)
  variant <- set_dat_param(variant, "condition file", condition_config_file)
  # We also need to (manually) set condition reference in the groups data
  variant@groups$condition.group <- 1

  # Save variant
  save_zvariant(variant, dir = file.path(zsetup_root, taxon$sheet),
                overwrite = TRUE, debug_msg = FALSE)
  create_sh_file(variant)

  # [XX]_[TX]_abf_wgt_con -----------------------------------------------------

  variant_id <- variant_id + 1
  variant <- get_variant(japan_zproject, variant_id)
  message("Editing variant: ", variant@name)

  # Set weights
  wgts <- lookup_weight(dat, featurenames(variant), w_field = w_field)

  # sppweights()<- not implemented yet.
  variant_spp_data <- sppdata(variant)
  variant_spp_data$weight <- wgts
  sppdata(variant) <- variant_spp_data

  # Set groups, use the same as with previous
  groups(variant) <- as.vector(groups)

  # Set cell removal rule
  variant <- set_dat_param(variant, "removal rule", 2)

  # Set groups use and groups file
  variant <- set_dat_param(variant, "use groups", 1)
  # Note that groups file location is always relative to the bat file
  groups_file <- file.path(variant@name, paste0(variant@name, "_groups.txt"))
  variant <- set_dat_param(variant, "groups file", groups_file)

  # Set PPA file
  variant <- set_dat_param(variant, "post-processing list file",
                           ppa_config_file)

  # Set condition. First, let's create the file itself (zonator
  # can't handle this yet). The file needs to be created only once per raxon
  # since all the variants can use the same file.
  condition_file_name <- file.path(zsetup_root, taxon$sheet,
                                   condition_config_file)
  condition_cmd_string <- paste(c(1, condition_raster_file), collapse = " ")
  write(condition_cmd_string, condition_file_name)
  # Need to define ppa_config.txt relative to the bat-file (same dir)-
  variant <- set_dat_param(variant, "use condition layer", 1)
  variant <- set_dat_param(variant, "condition file", condition_config_file)
  # We also need to (manually) set condition reference in the groups data
  variant@groups$condition.group <- 1

  # Save variant
  save_zvariant(variant, dir = file.path(zsetup_root, taxon$sheet),
                overwrite = TRUE, debug_msg = FALSE)
  create_sh_file(variant)

  # [XX]_[TX]_caz_wgt_con_hm3 -------------------------------------------------

  variant_id <- variant_id + 1
  variant <- get_variant(japan_zproject, variant_id)
  message("Editing variant: ", variant@name)

  # Set weights
  wgts <- lookup_weight(dat, featurenames(variant), w_field = w_field)
  # sppweights()<- not implemented yet.
  variant_spp_data <- sppdata(variant)
  variant_spp_data$weight <- wgts
  sppdata(variant) <- variant_spp_data

  # Set groups, use the same as with previous
  groups(variant) <- as.vector(groups)

  # Set cell removal rule
  variant <- set_dat_param(variant, "removal rule", 1)

  # Set groups use and groups file
  variant <- set_dat_param(variant, "use groups", 1)
  # Note that groups file location is always relative to the bat file
  groups_file <- file.path(variant@name, paste0(variant@name, "_groups.txt"))
  variant <- set_dat_param(variant, "groups file", groups_file)

  # Set PPA file
  variant <- set_dat_param(variant, "post-processing list file",
                           ppa_config_file)

  # Set condition. First, let's create the file itself (zonator
  # can't handle this yet). The file needs to be created only once per raxon
  # since all the variants can use the same file.
  condition_file_name <- file.path(zsetup_root, taxon$sheet,
                                   condition_config_file)
  condition_cmd_string <- paste(c(1, condition_raster_file), collapse = " ")
  write(condition_cmd_string, condition_file_name)
  # Need to define ppa_config.txt relative to the bat-file (same dir)-
  variant <- set_dat_param(variant, "use condition layer", 1)
  variant <- set_dat_param(variant, "condition file", condition_config_file)
  # We also need to (manually) set condition reference in the groups data
  variant@groups$condition.group <- 1

  # Set hierarchical mask hm2_raster_file
  variant <- set_dat_param(variant, "use mask", 1)
  variant <- set_dat_param(variant, "mask file", hm3_raster_file)

  # Save variant
  save_zvariant(variant, dir = file.path(zsetup_root, taxon$sheet),
                overwrite = TRUE, debug_msg = FALSE)
  create_sh_file(variant)

  # [XX]_[TX]_abf_wgt_con_hm3 -------------------------------------------------

  variant_id <- variant_id + 1
  variant <- get_variant(japan_zproject, variant_id)
  message("Editing variant: ", variant@name)

  # Set weights
  wgts <- lookup_weight(dat, featurenames(variant), w_field = w_field)
  # sppweights()<- not implemented yet.
  variant_spp_data <- sppdata(variant)
  variant_spp_data$weight <- wgts
  sppdata(variant) <- variant_spp_data

  # Set groups, use the same as with previous
  groups(variant) <- as.vector(groups)

  # Set cell removal rule
  variant <- set_dat_param(variant, "removal rule", 2)

  # Set groups use and groups file
  variant <- set_dat_param(variant, "use groups", 1)
  # Note that groups file location is always relative to the bat file
  groups_file <- file.path(variant@name, paste0(variant@name, "_groups.txt"))
  variant <- set_dat_param(variant, "groups file", groups_file)

  # Set PPA file
  variant <- set_dat_param(variant, "post-processing list file",
                           ppa_config_file)

  # Set condition. First, let's create the file itself (zonator
  # can't handle this yet). The file needs to be created only once per raxon
  # since all the variants can use the same file.
  condition_file_name <- file.path(zsetup_root, taxon$sheet,
                                   condition_config_file)
  condition_cmd_string <- paste(c(1, condition_raster_file), collapse = " ")
  write(condition_cmd_string, condition_file_name)
  # Need to define ppa_config.txt relative to the bat-file (same dir)-
  variant <- set_dat_param(variant, "use condition layer", 1)
  variant <- set_dat_param(variant, "condition file", condition_config_file)
  # We also need to (manually) set condition reference in the groups data
  variant@groups$condition.group <- 1

  # Set hierarchical mask hm2_raster_file
  variant <- set_dat_param(variant, "use mask", 1)
  variant <- set_dat_param(variant, "mask file", hm3_raster_file)

  # Save variant
  save_zvariant(variant, dir = file.path(zsetup_root, taxon$sheet),
                overwrite = TRUE, debug_msg = FALSE)
  create_sh_file(variant)

  variant_id <- variant_id + 1
}

taxon <- "taxa_all"

# Generate variants for all species together -------------------------------

variant <- get_variant(japan_zproject, variant_id)

# Generate groups
# Parse groups based on the Red-list status
groups <- lookup_rl_group(dat, featurenames(variant))
# We need to manually fix few things in groups
groups["Salamandrella_keyserlingii_comp"] <- groups["Salamandrella_keyserlingii"]

# [XX]_caz

message("Editing variant: ", variant@name)
# Set groups
groups(variant) <- as.vector(groups)
# Set groups use and groups file
variant <- set_dat_param(variant, "use groups", 1)
# Note that groups file location is always relative to the bat file
groups_file <- file.path(variant@name, paste0(variant@name, "_groups.txt"))
variant <- set_dat_param(variant, "groups file", groups_file)

# Set post-processing (LSM). First, let's create the file itself (zonator
# can't handle this yet). The file needs to be created only once per raxon
# since all the variants can use the same file.
ppa_file_name <- file.path(zsetup_root, taxon, ppa_config_file)
ppa_cmd_string <- paste(c("LSM", ppa_raster_file, 0, -1, 0), collapse = " ")
write(ppa_cmd_string, ppa_file_name)
# Need to define ppa_config.txt relative to the bat-file (same dir)-
variant <- set_dat_param(variant, "post-processing list file",
                         ppa_config_file)

# Save variant
save_zvariant(variant, dir = file.path(zsetup_root, taxon), overwrite = TRUE,
              debug_msg = FALSE)
create_sh_file(variant)

# [XX]_abf

variant_id <- variant_id + 1
variant <- get_variant(japan_zproject, variant_id)
message("Editing variant: ", variant@name)

# Set z-value on taxon level. First, we need to know how many species are in
# what taxon. HACK: this is not the safest way of doing the z-value allocation
# but will do.

spp_rasters <- list.files("../Data.150928/",
                          pattern = "[A-Z][a-z]+_.+\\.(tif|img)$",
                          recursive = TRUE)
# Split by path separator; subfolder is the taxon.
spp_rasters <- as.data.frame(do.call("rbind",
                                     strsplit(spp_rasters, .Platform$file.sep)))
names(spp_rasters) <- c("taxon_name", "species")

# Manipulate the spp data
variant_spp_data <- sppdata(variant)
sppdata(variant) <- variant_spp_data

# Set groups, use the same as with caz
groups(variant) <- as.vector(groups)

# Set cell removal rule
variant <- set_dat_param(variant, "removal rule", 2)

# Set groups use and groups file
variant <- set_dat_param(variant, "use groups", 1)
# Note that groups file location is always relative to the bat file
groups_file <- file.path(variant@name, paste0(variant@name, "_groups.txt"))
variant <- set_dat_param(variant, "groups file", groups_file)

# Set PPA file
variant <- set_dat_param(variant, "post-processing list file",
                         ppa_config_file)

# Save variant
save_zvariant(variant, dir = file.path(zsetup_root, taxon), overwrite = TRUE,
              debug_msg = FALSE)
create_sh_file(variant)

# [XX]_caz_wgt
variant_id <- variant_id + 1
variant <- get_variant(japan_zproject, variant_id)
message("Editing variant: ", variant@name)

dat_sk <- dat[dat$st.species == "Salamandrella_keyserlingii",]
dat_sk$st.species <- "Salamandrella_keyserlingii_comp"
dat <- rbind(dat, dat_sk)

# sppweights()<- not implemented yet.
variant_spp_data <- sppdata(variant)

variant_spp_data$weight <- dat[[w_field]]
sppdata(variant) <- variant_spp_data

# Set groups, use the same as with previous
groups(variant) <- as.vector(groups)

# Set cell removal rule
variant <- set_dat_param(variant, "removal rule", 1)

# Set groups use and groups file
variant <- set_dat_param(variant, "use groups", 1)
# Note that groups file location is always relative to the bat file
groups_file <- file.path(variant@name, paste0(variant@name, "_groups.txt"))
variant <- set_dat_param(variant, "groups file", groups_file)

# Set PPA file
variant <- set_dat_param(variant, "post-processing list file",
                         ppa_config_file)

# Save variant
save_zvariant(variant, dir = file.path(zsetup_root, taxon), overwrite = TRUE,
              debug_msg = FALSE)
create_sh_file(variant)

# [XX]_abf_wgt
variant_id <- variant_id + 1
variant <- get_variant(japan_zproject, variant_id)
message("Editing variant: ", variant@name)

# Manipulate sppdata
variant_spp_data <- sppdata(variant)
variant_spp_data$weight <- dat[[w_field]]
sppdata(variant) <- variant_spp_data

# Set groups, use the same as with previous
groups(variant) <- as.vector(groups)

# Set cell removal rule
variant <- set_dat_param(variant, "removal rule", 2)

# Set groups use and groups file
variant <- set_dat_param(variant, "use groups", 1)
# Note that groups file location is always relative to the bat file
groups_file <- file.path(variant@name, paste0(variant@name, "_groups.txt"))
variant <- set_dat_param(variant, "groups file", groups_file)

# Set PPA file
variant <- set_dat_param(variant, "post-processing list file",
                         ppa_config_file)

# Save variant
save_zvariant(variant, dir = file.path(zsetup_root, taxon), overwrite = TRUE,
              debug_msg = FALSE)
create_sh_file(variant)

# [XX]_[TX]_caz_wgt_con
variant_id <- variant_id + 1
variant <- get_variant(japan_zproject, variant_id)
message("Editing variant: ", variant@name)

# Manipulate sppdata
variant_spp_data <- sppdata(variant)
variant_spp_data$weight <- dat[[w_field]]
sppdata(variant) <- variant_spp_data

# Set groups, use the same as with previous
groups(variant) <- as.vector(groups)

# Set cell removal rule
variant <- set_dat_param(variant, "removal rule", 1)

# Set groups use and groups file
variant <- set_dat_param(variant, "use groups", 1)
# Note that groups file location is always relative to the bat file
groups_file <- file.path(variant@name, paste0(variant@name, "_groups.txt"))
variant <- set_dat_param(variant, "groups file", groups_file)

# Set PPA file
variant <- set_dat_param(variant, "post-processing list file",
                         ppa_config_file)

# Set condition. First, let's create the file itself (zonator
# can't handle this yet). The file needs to be created only once per raxon
# since all the variants can use the same file.
condition_file_name <- file.path(zsetup_root, taxon, condition_config_file)
condition_cmd_string <- paste(c(1, condition_raster_file), collapse = " ")
write(condition_cmd_string, condition_file_name)
# Need to define ppa_config.txt relative to the bat-file (same dir)-
variant <- set_dat_param(variant, "use condition layer", 1)
variant <- set_dat_param(variant, "condition file", condition_config_file)
# We also need to (manually) set condition reference in the groups data
variant@groups$condition.group <- 1

# Save variant
save_zvariant(variant, dir = file.path(zsetup_root, taxon), overwrite = TRUE,
              debug_msg = FALSE)
create_sh_file(variant)

# [XX]_[TX]_abf_wgt_con
variant_id <- variant_id + 1
variant <- get_variant(japan_zproject, variant_id)
message("Editing variant: ", variant@name)

# Manipulate sppdata
variant_spp_data <- sppdata(variant)
variant_spp_data$weight <- dat[[w_field]]
sppdata(variant) <- variant_spp_data

# Set groups, use the same as with previous
groups(variant) <- as.vector(groups)

# Set cell removal rule
variant <- set_dat_param(variant, "removal rule", 2)

# Set groups use and groups file
variant <- set_dat_param(variant, "use groups", 1)
# Note that groups file location is always relative to the bat file
groups_file <- file.path(variant@name, paste0(variant@name, "_groups.txt"))
variant <- set_dat_param(variant, "groups file", groups_file)

# Set PPA file
variant <- set_dat_param(variant, "post-processing list file",
                         ppa_config_file)

# Set condition. First, let's create the file itself (zonator
# can't handle this yet). The file needs to be created only once per raxon
# since all the variants can use the same file.
condition_file_name <- file.path(zsetup_root, taxon, condition_config_file)
condition_cmd_string <- paste(c(1, condition_raster_file), collapse = " ")
write(condition_cmd_string, condition_file_name)
# Need to define ppa_config.txt relative to the bat-file (same dir)-
variant <- set_dat_param(variant, "use condition layer", 1)
variant <- set_dat_param(variant, "condition file", condition_config_file)
# We also need to (manually) set condition reference in the groups data
variant@groups$condition.group <- 1

# Save variant
save_zvariant(variant, dir = file.path(zsetup_root, taxon), overwrite = TRUE,
              debug_msg = FALSE)
create_sh_file(variant)

# [XX]_[TX]_caz_wgt_con_hm2
variant_id <- variant_id + 1
variant <- get_variant(japan_zproject, variant_id)
message("Editing variant: ", variant@name)

# Manipulate sppdata
variant_spp_data <- sppdata(variant)
variant_spp_data$weight <- dat[[w_field]]
sppdata(variant) <- variant_spp_data

# Set groups, use the same as with previous
groups(variant) <- as.vector(groups)

# Set cell removal rule
variant <- set_dat_param(variant, "removal rule", 1)

# Set groups use and groups file
variant <- set_dat_param(variant, "use groups", 1)
# Note that groups file location is always relative to the bat file
groups_file <- file.path(variant@name, paste0(variant@name, "_groups.txt"))
variant <- set_dat_param(variant, "groups file", groups_file)

# Set PPA file
variant <- set_dat_param(variant, "post-processing list file",
                         ppa_config_file)

# Set condition. First, let's create the file itself (zonator
# can't handle this yet). The file needs to be created only once per raxon
# since all the variants can use the same file.
condition_file_name <- file.path(zsetup_root, taxon, condition_config_file)
condition_cmd_string <- paste(c(1, condition_raster_file), collapse = " ")
write(condition_cmd_string, condition_file_name)
# Need to define ppa_config.txt relative to the bat-file (same dir)-
variant <- set_dat_param(variant, "use condition layer", 1)
variant <- set_dat_param(variant, "condition file", condition_config_file)
# We also need to (manually) set condition reference in the groups data
variant@groups$condition.group <- 1

# Set hierarchical mask hm2_raster_file
variant <- set_dat_param(variant, "use mask", 1)
variant <- set_dat_param(variant, "mask file", hm2_raster_file)

# Save variant
save_zvariant(variant, dir = file.path(zsetup_root, taxon), overwrite = TRUE,
              debug_msg = FALSE)
create_sh_file(variant)

# [XX]_[TX]_abf_wgt_con_hm2
variant_id <- variant_id + 1
variant <- get_variant(japan_zproject, variant_id)
message("Editing variant: ", variant@name)

# Manipulate sppdata
variant_spp_data <- sppdata(variant)
variant_spp_data$weight <- dat[[w_field]]
sppdata(variant) <- variant_spp_data

# Set groups, use the same as with previous
groups(variant) <- as.vector(groups)

# Set cell removal rule
variant <- set_dat_param(variant, "removal rule", 2)

# Set groups use and groups file
variant <- set_dat_param(variant, "use groups", 1)
# Note that groups file location is always relative to the bat file
groups_file <- file.path(variant@name, paste0(variant@name, "_groups.txt"))
variant <- set_dat_param(variant, "groups file", groups_file)

# Set PPA file
variant <- set_dat_param(variant, "post-processing list file",
                         ppa_config_file)

# Set condition. First, let's create the file itself (zonator
# can't handle this yet). The file needs to be created only once per raxon
# since all the variants can use the same file.
condition_file_name <- file.path(zsetup_root, taxon, condition_config_file)
condition_cmd_string <- paste(c(1, condition_raster_file), collapse = " ")
write(condition_cmd_string, condition_file_name)
# Need to define ppa_config.txt relative to the bat-file (same dir)-
variant <- set_dat_param(variant, "use condition layer", 1)
variant <- set_dat_param(variant, "condition file", condition_config_file)
# We also need to (manually) set condition reference in the groups data
variant@groups$condition.group <- 1

# Set hierarchical mask hm2_raster_file
variant <- set_dat_param(variant, "use mask", 1)
variant <- set_dat_param(variant, "mask file", hm2_raster_file)

# Save variant
save_zvariant(variant, dir = file.path(zsetup_root, taxon), overwrite = TRUE,
              debug_msg = FALSE)
create_sh_file(variant)

# [XX]_[TX]_caz_wgt_con_hm3
variant_id <- variant_id + 1
variant <- get_variant(japan_zproject, variant_id)
message("Editing variant: ", variant@name)

# Manipulate sppdata
variant_spp_data <- sppdata(variant)
variant_spp_data$weight <- dat[[w_field]]
sppdata(variant) <- variant_spp_data

# Set groups, use the same as with previous
groups(variant) <- as.vector(groups)

# Set cell removal rule
variant <- set_dat_param(variant, "removal rule", 1)

# Set groups use and groups file
variant <- set_dat_param(variant, "use groups", 1)
# Note that groups file location is always relative to the bat file
groups_file <- file.path(variant@name, paste0(variant@name, "_groups.txt"))
variant <- set_dat_param(variant, "groups file", groups_file)

# Set PPA file
variant <- set_dat_param(variant, "post-processing list file",
                         ppa_config_file)

# Set condition. First, let's create the file itself (zonator
# can't handle this yet). The file needs to be created only once per raxon
# since all the variants can use the same file.
condition_file_name <- file.path(zsetup_root, taxon, condition_config_file)
condition_cmd_string <- paste(c(1, condition_raster_file), collapse = " ")
write(condition_cmd_string, condition_file_name)
# Need to define ppa_config.txt relative to the bat-file (same dir)-
variant <- set_dat_param(variant, "use condition layer", 1)
variant <- set_dat_param(variant, "condition file", condition_config_file)
# We also need to (manually) set condition reference in the groups data
variant@groups$condition.group <- 1

# Set hierarchical mask hm2_raster_file
variant <- set_dat_param(variant, "use mask", 1)
variant <- set_dat_param(variant, "mask file", hm3_raster_file)

# Save variant
save_zvariant(variant, dir = file.path(zsetup_root, taxon), overwrite = TRUE,
              debug_msg = FALSE)
create_sh_file(variant)

# [XX]_[TX]_abf_wgt_con_hm3
variant_id <- variant_id + 1
variant <- get_variant(japan_zproject, variant_id)
message("Editing variant: ", variant@name)

# Manipulate sppdata
variant_spp_data <- sppdata(variant)
variant_spp_data$weight <- dat[[w_field]]
sppdata(variant) <- variant_spp_data

# Set groups, use the same as with previous
groups(variant) <- as.vector(groups)

# Set cell removal rule
variant <- set_dat_param(variant, "removal rule", 2)

# Set groups use and groups file
variant <- set_dat_param(variant, "use groups", 1)
# Note that groups file location is always relative to the bat file
groups_file <- file.path(variant@name, paste0(variant@name, "_groups.txt"))
variant <- set_dat_param(variant, "groups file", groups_file)

# Set PPA file
variant <- set_dat_param(variant, "post-processing list file",
                         ppa_config_file)

# Set condition. First, let's create the file itself (zonator
# can't handle this yet). The file needs to be created only once per raxon
# since all the variants can use the same file.
condition_file_name <- file.path(zsetup_root, taxon, condition_config_file)
condition_cmd_string <- paste(c(1, condition_raster_file), collapse = " ")
write(condition_cmd_string, condition_file_name)
# Need to define ppa_config.txt relative to the bat-file (same dir)-
variant <- set_dat_param(variant, "use condition layer", 1)
variant <- set_dat_param(variant, "condition file", condition_config_file)
# We also need to (manually) set condition reference in the groups data
variant@groups$condition.group <- 1

# Set hierarchical mask hm2_raster_file
variant <- set_dat_param(variant, "use mask", 1)
variant <- set_dat_param(variant, "mask file", hm3_raster_file)

# Save variant
save_zvariant(variant, dir = file.path(zsetup_root, taxon), overwrite = TRUE,
              debug_msg = FALSE)
create_sh_file(variant)
