library(dplyr)
library(ggplot2)
library(tidyr)
library(zonator)

japan_project <- create_zproject(root = "zsetup/", debug = TRUE)

# We'll use variant 01 as a template
caz_01 <- get_variant(japan_zproject, 1)
caz_wgt_msk_06 <- get_variant(japan_project, 7)

# Let's rename the groups back to the IUCN categories
iucn_names <- c("1"="1-LC", "2"="2-NT", "3"="3-VU", "4"="4-EN", "5"="5-CR", 
                "6"="6-EW", "7"="7-EX", "8"="8-DD")

spps <- sppdata(caz_01) %>% 
  select(name, group)

spps$group_name <- iucn_names[spps$group]

# 01_caz ------------------------------------------------------------------

caz_01_curves <- curves(caz_01)

caz_01_curves_data <- as_data_frame(caz_01_curves) %>% 
  select(-cost, -min_pr, -ave_pr, -w_pr, -ext1, -ext2) %>%
  filter(pr_lost >= 0.9) %>% 
  gather(species, dist_rem, -pr_lost) %>% 
  mutate(group_name = spps[species, "group_name"])
  
p <- ggplot(caz_01_curves_data, aes(factor(group_name), dist_rem))
p + geom_boxplot() + ylab("Distribution remaining") + xlab("Group")

# 06_caz_wgt_msk --------------------------------------------------------------

caz_wgt_msk_06_curves <- curves(caz_wgt_msk_06)

caz_wgt_msk_06_curves_data <- as_data_frame(caz_wgt_msk_06_curves) %>% 
  select(-cost, -min_pr, -ave_pr, -w_pr, -ext1, -ext2) %>%
  filter(pr_lost >= 0.9) %>% 
  gather(species, dist_rem, -pr_lost) %>% 
  mutate(group_name = spps[species, "group_name"])

p <- ggplot(caz_wgt_msk_06_curves_data, aes(factor(group_name), dist_rem))
p + geom_boxplot() + ylab("Distribution remaining") + xlab("Group")
