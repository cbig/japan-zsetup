library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(moments)
library(R.cache)
library(RColorBrewer)
library(readr)
library(purrr)
library(tidyr)
library(viridis)
library(zonator)

# Common data -------------------------------------------------------------

# Red-list status used as a basis for grouping
rl_groups <- list(LC = 1, NT = 2, VU = 3, EN = 4, CR = 5, EW = 6, EX = 7, DD = 8)

# Common functions ---------------------------------------------------------

calc_stats <- function(x, mult = 1) {
  quant <- quantile(x)
  x <- na.omit(x)
  min <- min(x)
  q25 <- quant[2]
  sd <- sd(x)
  median <- median(x)
  mean <- mean(x)
  q75 <- quant[4]
  max <- max(x)
  skewness <- moments::skewness(x)
  kurtosis <- moments::kurtosis(x)
  zeros <- sum(x == 0.0)
  return(data.frame(min = min, q25 = q25, sdmin = mean - sd,
                    mean = mean, median = median,
                    sdmax = mean + sd, q75 = q75, max = max,
                    skewness = skewness, kurtosis = kurtosis,
                    zeros = zeros))
}

calc_row_stats <- function(x, groups) {
  # Use variant name as variant label
  variant <- x@name
  x_curves <- zonator::curves(x)
  pr_lost <- x_curves$pr_lost
  x_curves <- x_curves %>%
    dplyr::select(-pr_lost, -cost, -min_pr, -ave_pr, -w_pr, -ext1, -ext2)
  # Use groups if provided
  if (groups) {
    x_groups <- zonator::groups(x)
    x_group_stats <- list()
    x_groupnames <- zonator::groupnames(x)
    for (group in unique(x_groups)) {
      group_name <- paste0(variant, "_", group)
      x_data <- x_curves[,x_groups == group]
      if (!is.null(dim(x_data))) {
        x_data <- cbind(pr_lost = pr_lost,
                        dplyr::bind_rows(apply(x_data, 1, calc_stats)))
        x_data$variant <- group_name
        x_group_stats[[group]] <- x_data
      }
    }
    x_group_stats <- dplyr::bind_rows(x_group_stats)
  } else {
    x_group_stats <- cbind(pr_lost = pr_lost,
                           dplyr::bind_rows(apply(x_curves, 1, calc_stats)))
    x_group_stats$variant <- variant
  }
  if (min(x_group_stats$sdmin) < 0.0) {
    x_group_stats[x_group_stats$sdmin < 0.0,]$sdmin <- 0.0
  }
  if (max(x_group_stats$sdmax) > 1.0) {
    x_group_stats[x_group_stats$sdmax > 1.0,]$sdmax <- 1.0
  }
  return(x_group_stats)
}

create_sh_file <- function(x) {
  if (class(x) == "Zvariant") {
    bat_file <- x@bat.file
  } else {
    bat_file <- x
  }
  
  sh_file <- gsub("\\.bat", "\\.sh", bat_file)
  
  cmd_lines <- readLines(bat_file)
  new_cmd_lines <- c("#!/bin/sh")
  
  for (line in cmd_lines) {
    line <- gsub("call ", "", line)
    line <- gsub("\\.exe", "", line)
    new_cmd_lines <- c(new_cmd_lines, line)
  }
  
  file_con <- file(sh_file)
  writeLines(new_cmd_lines, file_con)
  close(file_con)
  Sys.chmod(sh_file)
  return(invisible(TRUE))
}

get_stat_curves <- function(zproject, variant_ids, groups = FALSE) {

  curves_data <- variant_ids %>%
    purrr::map(zonator::get_variant, x = zproject) %>%
    purrr::map(calc_row_stats, groups = groups) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(variant = factor(variant, levels = unique(variant),
                                   ordered = TRUE))
  return(curves_data)
}

get_perf_levels <- function(zproject, variant_ids, perf_levels) {
  get_level <- function(plevel, dat) {
    plevel_pr_lost <- 1 - plevel
    dat <- dat %>%
      dplyr::group_by(variant) %>%
      dplyr::filter(pr_lost >= plevel_pr_lost) %>%
      dplyr::slice(1) %>%
      dplyr::mutate(perf_level = plevel) %>%
      dplyr::ungroup()
    return(dat)
  }
  curves_dat <- get_stat_curves(zproject, variant_ids = variant_ids)
  curves_dat <- purrr::map(perf_levels, get_level, dat = curves_dat) %>%
    dplyr::bind_rows()
  return(curves_dat)
}

#' Get PPA data per planning unit based on existing variants.
#'
#' Zonation variants must have been run and results must exist and match
#' with the auxiliary PLU data provided.
#'
#' @param zproject Zproject object containing the variants.
#' @param prefectures_file Character path to the auxiliary PLU data.
#' @param variant_id Numeric identifying the variant used.
#' @param label Character label to be used.
#'
#' @return data.frame
#'
get_ppa_data <- function(zproject, prefectures_file,
                         variant_id, label) {

  # Read in the prefectures data
  prefectures <- readr::read_csv(prefectures_file)

  variant <- zonator::get_variant(zproject, variant_id)
  res <- zonator::results(variant)
  ppa_data <- zonator::ppa_lsm(res) %>%
    dplyr::left_join(., prefectures, by = c("Unit" = "no")) %>%
    dplyr::select(id = Unit, prefecture, area = Area, mrank = Mean_rank,
                  sdsum = Spp_distribution_sum, plus10 = Plus_10,
                  plus_1 = Plus_1, plus_01 = Plus_01, plus_0001 = Plus_0001) %>%
    # Calculate ranks for prefectures based on mean rank
    dplyr::mutate(norm_sdsum = (sdsum / area),
                  prefecture_rank = dplyr::dense_rank(-mrank), label = label)
  return(ppa_data)
}

#' Create a set of variant on file system.
#'
#' @param zsetup_root Character path to the parent dir of variants.
#' @param variant_templates Character vector containing the templates for variants
#'   to be created.
#' @param spp_data List holding auxiliary data for taxa.
#' @param data_dir Character path to data directory.
#' @param prefix_spp_paths Prefix path inserted into spp files.
#' @param dat_template_file Character path to template dat file used.
#'
#' @return Zproject object
#'
initiate_zproject <- function(zsetup_root = "zsetup", variant_templates,
                              spp_data, data_dir = "../Data.150928",
                              prefix_spp_paths = "../..",
                              dat_template_file = "templates/template.dat") {
  variant_id <- 1

  for (taxon in names(spp_data)) {
    variant_names <- variant_templates
    sub_project_dir <- file.path(zsetup_root, taxon)
    message("Creating sub-project in ", sub_project_dir)
    
    # Generate ids dynamically
    ids <- sprintf("%02d_",
                   variant_id:(variant_id + length(variant_names) - 1))

    # Generate actual variant names
    variant_names <- paste0(ids,
                            gsub(pattern = "\\[TX\\]",
                                 replacement = spp_data[[taxon]]$code,
                                 variant_names))
    
    # NOTE: set the data path manually
    create_zproject(name = sub_project_dir, dir = ".",
                    variants = variant_names,
                    dat_template_file = dat_template_file,
                    spp_template_dir = file.path(data_dir, spp_data[[taxon]]$sheet),
                    override_path = file.path(prefix_spp_paths, data_dir,
                                              spp_data[[taxon]]$sheet),
                    overwrite = TRUE, debug = TRUE)

    # Update variant ID for the next taxon
    variant_id <- variant_id + length(variant_names)
    message(" Created variants: \n  ", paste(variant_names, collapse = "\n  "))
  }

  # Finally, create variants with all species together
  sub_project_dir <- file.path(zsetup_root, "taxa_all")
  message("Creating sub-project in ", sub_project_dir)
  
  # No need for the taxa names anymore
  variant_names <- gsub("\\[TX\\]_", "", variant_templates)
  ids <- sprintf("%02d_", variant_id:(variant_id + length(variant_names) - 1))
  variant_names <- paste0(ids, variant_names)
  # Use "recursive = TRUE" with specific spp name template to get all the
  # rasters from different input folders
  create_zproject(name = sub_project_dir,
                  dir = ".",
                  variants = variant_names,
                  dat_template_file = dat_template_file,
                  spp_template_dir = data_dir,
                  override_path = file.path(prefix_spp_paths, data_dir),
                  overwrite = TRUE,
                  debug = FALSE,
                  recursive = TRUE,
                  spp_file_pattern = "[a-z]{3}_.+\\.(tif|img)$")
  message(" Created variants: \n  ", paste(variant_names, collapse = "\n  "))

  return(load_zproject(zsetup_root))
}

#' Load Zonation results.
#'
#' @param zsetup_root Character path to the parent dir of variants.
#' @param cache
#' @param ... other arguments passed on to \code{load_zproject()}.
#'
#' @return Zproject object
#'
.load_zproject <- function(zsetup_root = "zsetup", cache = TRUE, ...) {


  if (cache) {
    # Create the key
    key <- list(zsetup_root)
    zproject <- loadCache(key)
    if (!is.null(zproject)) {
      message("Loaded cached data")
      return(zproject)
    } else {
      warning("Couldn't find cached data, reloading...")
      zproject <- zonator::load_zproject(zsetup_root, ...)
      saveCache(zproject, key = key)
    }
  } else {
    zproject <- zonator::load_zproject(zsetup_root, ...)
  }
  return(zproject)
}

# Lookup a numeric value for red-list group for a species.
# @param table Dataframe table containing the data for a given taxon.
# @param species String character defining species.
lookup_rl_group <-
  Vectorize(function(table, species) {

    if (!"st.species" %in% names(table) | !"category" %in% names(table)) {
      stop("Table must have columns 'st.species' and 'category'")
    }
    if (!species %in% table$st.species) {
      warning(paste0("Species provided (", species, ") not found in the table."))
      return(NA)
    }
    if (!all(table$category %in% names(rl_groups))) {
      stop("All category values in table must be in: ", paste(names(rl_groups),
                                                              collapse = ", "))
    }
    
    rl_cat <- table[which(table$st.species == species), ]$category
    rl_cat <- unlist(rl_groups[rl_cat])
    # Remove vector names
    rl_cat <- as.vector(rl_cat)
    return(rl_cat)
  },
  c("species"), USE.NAMES = TRUE)

# Lookup a numeric value for weight for a species.
# @param table Dataframe table containing the data for a given taxon.
# @param species String character defining species.
# @param w_field String character name of the weight field used.
lookup_weight <-
  Vectorize(function(table, species, w_field) {

    if (!"st.species" %in% names(table) | !"weight" %in% names(table)) {
      stop("Table must have columns 'st.species' and 'weight'")
    }
    if (!species %in% table$st.species) {
      warning(paste0("Species provided (", species, ") not found in the table."))
      return(NA)
    }
    if (!species %in% table$st.species) {
      warning(paste0("Species provided (", species, ") not found in the table."))
      return(NA)
    }
    
    weight <- as.numeric(table[which(table$st.species == species), ][w_field])
    if (length(weight) > 1) {
      weight <- weight[1]
      warning("Multiple weights found for a single species, using the first: ",
              weight)
    }
    return(weight)
  },
  c("species"), USE.NAMES = TRUE)

plot_curves <- function(stats, title = NULL, non_param = FALSE,
                        invert_x = FALSE, labels = NULL, plot_min = FALSE,
                        plot_max = FALSE, wrap=FALSE, ncol=NULL, nrow=NULL, 
                        highlights = NULL, max_y = NULL,
                        ylab = "Mean fraction of species' ranges covered\n", 
                        xlab = "\nFraction of the landscape protected")    {

  if (non_param) {
    dat <- data.frame(x = stats$pr_lost,
                      ymin = stats$q25,
                      ymax = stats$q75,
                      ystat = stats$median,
                      min = stats$min,
                      max = stats$max)
  } else {
    dat <- data.frame(x = stats$pr_lost,
                      ymin = stats$sdmin,
                      ymax = stats$sdmax,
                      ystat = stats$mean,
                      min = stats$min,
                      max = stats$max)
  }

  dat$variant <- stats$variant
  # Set max if not provided
  if (is.null(max_y)) {
    max_y <- 1.0
  }
  
  # Set labels
  if (!is.null(labels)) {
    if (length(labels) != length(levels(dat$variant))) {
      warning("Number of labels not the same as variant levels, ignoring")
    } else {
      dat$variant <- factor(dat$variant, levels = levels(dat$variant),
                            labels = labels, ordered = TRUE)
    }
  }


  if (invert_x) {
    dat$x <- 1.0 - dat$x
    x_axis <- scale_x_continuous(breaks = seq(0, 1, 0.2),
                                 labels = paste0(seq(0, 1, 0.2) * 100, "%"))
  } else {
    x_axis <- scale_x_continuous(breaks = seq(0, 1, 0.2),
                                 labels = paste0(seq(1, 0, -0.2) * 100, "%"))
  }

  lookup <- function(labels) {
    labels <- sapply(1:length(labels),
                     function(x) paste0(LETTERS[x], ": ", labels[x]))
    return(labels)
  }

  p <- ggplot(dat, aes(x = x))

  if (!is.null(highlights)) {

    line_alpha <- 0.6
    linetype <- 3
    label_y <- 1
    label_size <- 3

    #colors <- ggplot2::scale_color_brewer(length(highlights))
    colors <- RColorBrewer::brewer.pal(6, name = "RdYlBu")[1:length(highlights)]
    
    for (i in 1:length(highlights)) {
      color <- colors[i]
      if (i == 1) {
        xmin <- 0
      } else {
        xmin <- highlights[i - 1]

      }
      xmax <- highlights[i]
      p <- p + annotate("rect", fill = color, alpha = 0.3,
                        xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf)
    }

    for (i in 1:length(highlights)) {
      p <- p + geom_vline(xintercept = highlights[i], alpha = line_alpha,
                          linetype = linetype)
    }

    #for (i in 1:length(highlights)) {
    #  label <- paste0(round(highlights[i] * 100, 0), "%")
    #  p <- p + annotate("text", x = highlights[i], y = label_y, label = label,
    #                    size = label_size)
    #}
  }
  
  if (wrap) {
    p <- p + geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "grey60", alpha = 0.6) +
      geom_line(aes(y = ystat))
  } else {
    p <- p + geom_line(aes(y = ystat, color = variant), size = 0.6)
  }

  if (plot_min) {
    p <- p + geom_line(aes(y = min), linetype = 2, size = 0.4, color = "black")
  }
  if (plot_max) {
    p <- p + geom_line(aes(y = max), linetype = 4, size = 0.4, color = "black")
  }

  p <- p + scale_y_continuous(breaks = seq(0, 1, 0.2),
                              labels = paste0(seq(0, 1, 0.2) * 100, "%"),
                              limits = c(0, max_y)) +
    x_axis + ggtitle(title) +
    ylab(ylab) +
    xlab(xlab)
  if (wrap) {
    if (is.null(nrow) | is.null(ncol)) {
      stop("nrow and ncol must be provided when wrapping")
    }
    p <- p + facet_wrap(~variant, ncol = ncol, nrow = nrow, 
                        labeller = labeller(variant = lookup))
  } else {
    if (length(levels(dat$variant)) > 1) {
      p <- p + scale_colour_brewer("", type = "qual", palette = 2, direction = 1)
    } else {
      p <- p + scale_color_manual("", values = "black")
    }
  }
    
  p <- p + theme_ipsum_rc() + theme(strip.text = element_text(hjust = 0.5))
  return(p)
}

#' Regroup features to taxa from RL categories.
#'
#' All analysis use groups based on the RL category (1-8). Regroup
#' features to taxon
#'
#' @note Only works for taxa_all variants.
#'
#' @param x Zvariant object (with groups).
#'
#' @return Regroupped Zvariant object.
#'
#' @author Joona Lehtomaki \email{joona.lehtomaki@@gmail.com}
#'
regroup_to_taxa <- function(x) {

  taxa_groups <- c(rep(1, 72), rep(2, 337), rep(3, 210), rep(4, 101), rep(5, 5533),
                   rep(6, 73))
  if (length(taxa_groups) != nrow(sppdata(x))) {
    stop("Number of taxa groups and rows in sppdata differ")
  }

  taxa_group_names <- c("1" = "Amphibians", "2" = "Birds",
                        "3" = "Freshwater_fish", "4" = "Mammals",
                        "5" = "Plants", "6" = "Reptiles")

  groups(x) <- taxa_groups
  groupnames(x) <- taxa_group_names
  return(x)
}

#' Write Zonation groups file based on a Zvariant object.
#'
#' Function accesses \code{groups} slot of a \code{Zvariant} object and writes
#' the content into a file (if groups are assigned). Accessing Zvariant object's
#' slost directly not a generally good idea. Done here until zonator actually
#' supports writing out files properly (does not yet in 0.4.1).
#'
#' @note This functionality should eventually be incorporated into
#'  \code{zonator}.
#'
#' @param x Zvariant object (with groups).
#' @param filename String character file path to be written.
#' @parma overwrite Logical indicating if an existing file should be overwritten.
#'
#' @return Invisible NULL.
#'
#' @author Joona Lehtomaki \email{joona.lehtomaki@@gmail.com}
#'
write_groups <- function(x, filename, overwrite = FALSE) {
  if (!class(x) == "Zvariant") {
    stop("Object provided must be of Zvariant class")
  }
  if (all(dim(x@groups) == c(0, 0))) {
    stop("Zvariant must have groups.")
  }
  groups <- x@groups
  # Strip out "name" (last) column
  groups <- groups[, -ncol(groups)]
  if (!file.exists(filename) & !overwrite) {
    stop("File exists and overwrite is off.")
  }
  write.table(groups, file = filename, sep = "\t", row.names = FALSE,
              col.names = FALSE)
  message("Wrote groups file ", filename)
  return(invisible(NULL))
}
