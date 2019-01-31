library(dplyr)
library(Hmisc)
library(readr)
library(readxl)
library(zonator)

source("R/00_lib/utils.R")


# Read in data ------------------------------------------------------------

# Load the project, creates japan_zproject. WARNING: loaded object is cached
# using memoise-package. If results have change, don't load the cached version.
# Also if the cache needs to updated, you will have to do this manually.
zproject_japan <- .load_zproject("zsetup", cache = TRUE, debug = TRUE)

# Read in the auxiliary data for different taxa
read_taxon <- function(taxon_name) {
  # Read in taxon data
  dat <- readxl::read_excel("../Data.150928/AUC&RL.xlsx", taxon_name)
  # Add taxon name as a column
  dat$taxon <- taxon_name
  return(dat)
}

all_taxa <- lapply(c("plants", "mammals", "birds", "reptiles", "amphibians", "freshwater_fish"),
                   read_taxon)
all_taxa <- dplyr::bind_rows(all_taxa)
all_taxa$name <- Hmisc::capitalize(gsub(" ", "_", all_taxa$st.species))

all_taxa_weights <- readr::read_csv("../Data.150928/all_taxa_weights.csv")
all_taxa_weights <- all_taxa_weights %>% 
  dplyr::select(st.species, AGG = agg_value, WGT = weight)

all_taxa <- dplyr::left_join(all_taxa, all_taxa_weights, by = c("st.species" = "st.species"))

# Check variables ---------------------------------------------------------

# panel.smooth function is built in.
# panel.cor puts correlation in upper panels, size proportional to correlation
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste(prefix, txt, sep = "")
  if (missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}


plot_pairs <- function(dat, taxon_name = NULL) {

  if (is.null(taxon_name)) {
    title <- "All taxa"
    dat <- dat %>% 
      dplyr::select(ED, END, THS, AGG, WGT, AUC)
  } else {
    title <- taxon_name
    dat <- dat %>%
      dplyr::filter(taxon == taxon_name) %>% 
      dplyr::select(ED, END, THS, WGT, AUC)
  }
  
  pairs(dat,lower.panel = panel.smooth, upper.panel = panel.cor, 
        pch = 20, main = title)  

}

plot_pairs(all_taxa)
plot_pairs(all_taxa, "amphibians")
plot_pairs(all_taxa, "birds")
plot_pairs(all_taxa, "freshwater_fish")
plot_pairs(all_taxa, "mammals")
plot_pairs(all_taxa, "plants")
plot_pairs(all_taxa, "reptiles")

# Regroup -----------------------------------------------------------------

plot_variant <- function(variant) {
  grp_curves <- zonator::curves(variant, groups = TRUE)
  plot(grp_curves, invert.x = TRUE, min = TRUE, mean = TRUE, max = TRUE)
}

regroup_endemic <- function(variant) {
  # Combine data
  spp_data <- zonator::sppdata(variant) %>% 
    dplyr::left_join(., all_taxa, by = c("name" = "name"))
  
  # Amphibians need to be manually patched for Salamandrella_keyserlingii_comp
  spp_data[spp_data$name == "Salamandrella_keyserlingii_comp",]$END <- spp_data[spp_data$name == "Salamandrella_keyserlingii",]$END 
  # Renumber
  spp_data$END <- ifelse(spp_data$END == 4, 1, 2)
  
  zonator::groups(variant) <- spp_data$END
  zonator::groupnames(variant) <- c("1" = "Endemic", "2" = "Non-endemic")
  return(variant)
}

regroup_ed <- function(variant) {
  # Combine data
  spp_data <- zonator::sppdata(variant) %>% 
    dplyr::left_join(., all_taxa, by = c("name" = "name"))
  
  # Amphibians need to be manually patched for Salamandrella_keyserlingii_comp
  spp_data[spp_data$name == "Salamandrella_keyserlingii_comp",]$ED <- spp_data[spp_data$name == "Salamandrella_keyserlingii",]$ED 
  # Renumber
  spp_data$ED <- ifelse(spp_data$ED >= quantile(spp_data$ED)[3], 1, 2)
  
  zonator::groups(variant) <- spp_data$ED
  zonator::groupnames(variant) <- c("1" = "Evolutionary distinct", 
                                    "2" = "Evolutionary indistinct")
  return(variant)
}


# Endemicity --------------------------------------------------------------


caz <- zonator::get_variant(zproject_japan, 61)
caz <- regroup_endemic(caz)
plot_variant(caz)

caz_wgt <- zonator::get_variant(zproject_japan, 63)
caz_wgt <- regroup_endemic(caz_wgt)
plot_variant(caz_wgt)

caz_wgt_con <- zonator::get_variant(zproject_japan, 65)
caz_wgt_con <- regroup_endemic(caz_wgt_con)
plot_variant(caz_wgt_con)

caz_wgt_con_hm2 <- zonator::get_variant(zproject_japan, 67)
caz_wgt_con_hm2 <- regroup_endemic(caz_wgt_con_hm2)
plot_variant(caz_wgt_con_hm2)


# Evolutionary distinctiviness --------------------------------------------

caz <- zonator::get_variant(zproject_japan, 61)
caz <- regroup_ed(caz)
plot_variant(caz)

caz_wgt <- zonator::get_variant(zproject_japan, 63)
caz_wgt <- regroup_ed(caz_wgt)
plot_variant(caz_wgt)

caz_wgt_con <- zonator::get_variant(zproject_japan, 65)
caz_wgt_con <- regroup_ed(caz_wgt_con)
plot_variant(caz_wgt_con)

caz_wgt_con_hm2 <- zonator::get_variant(zproject_japan, 67)
caz_wgt_con_hm2 <- regroup_endemic(caz_wgt_con_hm2)
plot_variant(caz_wgt_con_hm2)
