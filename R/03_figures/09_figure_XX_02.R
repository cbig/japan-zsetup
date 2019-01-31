# Figure XX: Surrogacy analysis from pre-loaded variants.
library(dplyr)
library(ggplot2)
library(hrbrthemes)
library(purrr)
library(tidyr)
library(viridis)
library(zonator)

source("R/00_lib/utils.R")


# Helper functions --------------------------------------------------------

build_grp_data <- function(variant_list) {
  grp_curve_dfs <- lapply(X = variant_list, FUN = zonator::curves, 
                          groups = TRUE)
  grp_curve_dfs <- lapply(seq_along(grp_curve_dfs), 
                          function(x) {
                            x_name <- names(grp_curve_dfs[x])
                            main_taxon <- unlist(strsplit(x_name, "_"))[2]
                            dat <- grp_curve_dfs[[x]] 
                            dat$main_taxon <- main_taxon
                            if (grepl("_con$", x_name)) {
                              dat$variant <- "Condition"
                            } else {
                              dat$variant <- "No condition"
                            }
                            return(dat)
                          })
  all_dat <- dplyr::bind_rows(grp_curve_dfs) %>% 
    dplyr::select(pr_lost, variant, main_taxon, starts_with("mean.")) %>% 
    tidyr::gather(taxon, pr_rem, -pr_lost, -variant, -main_taxon) %>% 
    dplyr::mutate(taxon = gsub("mean\\.", "", taxon))
  return(all_dat)
}

# Load the project, creates japan_zproject. WARNING: loaded object is cached
# using memoise-package. If results have change, don't load the cached version.
# Also if the cache needs to updated, you will have to do this manually.kus
# > key <- list("zsetup")
# > saveCache(zproject_japan, key = key)
zproject_japan <- .load_zproject("zsetup", debug = TRUE, cache = TRUE)

taxa_group_names <- c("1" = "Amphibians", "2" = "Birds", 
                      "3" = "Freshwater fish", "4" = "Mammals", 
                      "5" = "Plants", "6" = "Reptiles")

# Select a subset of variants
variant_names <- c("83_amp_all_caz_wgt", "85_amp_all_caz_wgt_con",
                   "93_bir_all_caz_wgt", "95_bir_all_caz_wgt_con",
                   "103_frf_all_caz_wgt", "105_frf_all_caz_wgt_con",
                   "113_mam_all_caz_wgt", "115_mam_all_caz_wgt_con",
                   "123_pla_all_caz_wgt", "125_pla_all_caz_wgt_con",
                   "133_rep_all_caz_wgt", "135_rep_all_caz_wgt_con")
# Match the numeric variant IDs.
variant_ids <- match(variant_names, names(zproject_japan))

# Subset the variants and set group names
sel_variants <- variant_ids %>%
  purrr::map(zonator::get_variant, x = zproject_japan) %>% 
  purrr::map(zonator::`groupnames<-`, taxa_group_names)
  
curves_data <- build_grp_data(sel_variants)
curves_data$variant <- factor(curves_data$variant, 
                              levels = c("No condition", "Condition"),
                              ordered = TRUE)

# For invert x-axis
curves_data$pr_prot <- 1.0 - curves_data$pr_lost

# FOR NOW, only select the no-condition variants. Keep the above code
# in case it's needed later.
curves_data <- dplyr::filter(curves_data, variant == "No condition") %>% 
  dplyr::select(-variant)

# Plot --------------------------------------------------------------------

x_lab <- "\nFraction of the landscape protected"
x_scale <- scale_x_continuous(breaks = seq(0, 1, 0.2),
                              labels = paste0(seq(0, 1, 0.2) * 100, "%"))
y_scale <- scale_y_continuous(breaks = seq(0, 1, 0.2),
                              labels = paste(100 * seq(0, 1, 0.2), "%"))

boost_label <- function(string) {
  lookup <- list("amp" = "Amphibians", "bir" = "Birds", 
                 "frf" = "Freshwater fish", "mam" = "Mammals",
                 "pla" = "Plants", "rep" = "Reptiles")
  return(sapply(string, function(x) lookup[[x]], USE.NAMES = FALSE))
}

p1 <- ggplot2::ggplot(curves_data, aes(x = pr_prot, y = pr_rem, color = taxon)) +
  geom_line(size = 0.5) + xlab(x_lab) + 
  facet_wrap(~main_taxon, labeller = labeller(main_taxon = boost_label),
             nrow = 2, ncol = 3) +
  scale_color_manual("", values =  rev(viridis(6))) +
  x_scale + y_scale + ylab("Average feature distribution covered\n") +
  theme_ipsum_rc() + theme(legend.position = "top")

ggsave("figs/figure_S03.png", units = "cm", width = 18, height = 16)
