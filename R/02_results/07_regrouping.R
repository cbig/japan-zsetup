library(zonator)

source("R/00_lib/utils.R")

iucn_names <- c("1" = "LC", "2" = "NT", "3" = "VU", "4" = "EN", "5" = "CR", 
                "6" = "EW", "7" = "EX", "8" = "DD")

zproject_japan <- .load_zproject("zsetup", debug = TRUE, cache = TRUE)

variant65 <- get_variant(zproject_japan, 65)

taxa_groups <- c(rep(1, 72), rep(2, 337), rep(3, 212), rep(4, 104), rep(5, 5565),
                 rep(6, 73))
taxa_groups_name <- c("1" = "Amp", "2" = "Bir", "3" = "FWF", "4" = "Mam", 
                      "5" = "Pla", "6" = "rep")

groups(variant65) <- taxa_groups
groupnames(variant65) <- taxa_groups_name

grp_curves_65 <- curves(variant65, groups = TRUE)
plot(grp_curves_65)


# EXTRA -------------------------------------------------------------------

# Write the groups data out
write.table(variant65@groups[,1:5], 
            file = "zsetup/cross_loading/71_load_caz_amp_all/71_load_caz_amp_all_groups.txt", 
            col.names = FALSE, row.names = FALSE, sep = " ")
