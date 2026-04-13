#!/usr/bin/env Rscript

library(optparse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(scales)
library(stringr)
library(RColorBrewer)
library(grid)
library(patchwork)
library(svglite)

############################################################
# Arguments
############################################################

option_list <- list(
  make_option(c("--input"),       type="character"),
  make_option(c("--output"),      type="character"),
  make_option(c("--output2"),     type="character"),
  make_option(c("--protein"),     type="character"),
  make_option(c("--genome_file"), type="character", default=NULL,
              help="CSV with columns: genome, domain, phylum, class, order, family, genus"),
  make_option(c("--tax_level"),   type="character", default=NULL,
              help="Taxonomic level to split plots by. One of: domain, phylum, class, order, family, genus")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input))   stop("ERROR: --input not provided")
if (is.null(opt$output))  stop("ERROR: --output not provided")
if (is.null(opt$output2)) stop("ERROR: --output2 not provided")
if (is.null(opt$protein)) stop("ERROR: --protein not provided")

valid_tax_levels <- c("domain", "phylum", "class", "order", "family", "genus")
if (!is.null(opt$genome_file) && is.null(opt$tax_level)) {
  stop("ERROR: --tax_level must be provided when --genome_file is used")
}
if (!is.null(opt$tax_level) && !(opt$tax_level %in% valid_tax_levels)) {
  stop(paste("ERROR: --tax_level must be one of:", paste(valid_tax_levels, collapse=", ")))
}

# Derive color level as one step below the split tax level in the hierarchy
tax_hierarchy <- c("species", "genus", "family", "order", "class", "phylum", "domain")

if (!is.null(opt$tax_level)) {
  tax_level_idx <- which(tax_hierarchy == opt$tax_level)

  if (length(tax_level_idx) == 0) {
    stop(paste("ERROR: --tax_level not found in hierarchy:", opt$tax_level))
  }
  if (tax_level_idx == 1) {
    stop("ERROR: cannot color one level below species — no lower level defined")
  }

  color_level <- tax_hierarchy[tax_level_idx - 1]
} else {
  color_level <- "order"   # sensible default when no genome_file / tax_level given
}

protein_name <- opt$protein

############################################################
# Load data
############################################################

df <- read_csv(opt$input, show_col_types = FALSE)

required_cols <- c("annotation_cluster", "cluster_label", "genome",
                   "IPS_acc", "IPS_desc", "gene_offset", "distance_bp")
missing <- setdiff(required_cols, colnames(df))
if (length(missing) > 0) stop(paste("Missing columns:", paste(missing, collapse=", ")))

############################################################
# Build a distinguishable colour palette for a character
# vector of names (alphabetically pre-sorted by the caller).
#
# Uses a curated high-contrast set:
#   up to 20 values  -> hand-picked 20-colour qualitative set
#   beyond 20        -> interpolated extension of the same set
############################################################

DISTINCT_20 <- c(
  "#E63946", "#F4A261", "#2A9D8F", "#457B9D", "#6A0572",
  "#F1C453", "#264653", "#A8DADC", "#E76F51", "#8338EC",
  "#06D6A0", "#FFB703", "#219EBC", "#FB8500", "#023047",
  "#8ECAE6", "#C77DFF", "#D62828", "#52B788", "#B5838D"
)

make_color_palette <- function(names_vec) {
  n <- length(names_vec)
  cols <- if (n <= 20) {
    DISTINCT_20[seq_len(n)]
  } else {
    colorRampPalette(DISTINCT_20)(n)
  }
  setNames(cols, names_vec)
}

############################################################
# Join order-level taxonomy from genome_file
# Also joins tax_label for the per-taxon split if requested.
############################################################

if (!is.null(opt$genome_file)) {
  tax_df <- read_csv(opt$genome_file, show_col_types = FALSE)

  # genome_file uses "genome_file" as the genome column header
  needed <- unique(c("genome_file", color_level, opt$tax_level))
  missing_tax <- setdiff(needed, colnames(tax_df))
  if (length(missing_tax) > 0) {
    stop(paste("genome_file is missing columns:", paste(missing_tax, collapse=", ")))
  }

  # Build the global color palette from ALL values in the color_level column
  # (alphabetically sorted, before any merge losses) so colors are consistent
  # across every taxon plot produced in this run.
  all_color_values <- sort(unique(na.omit(tax_df[[color_level]])))
  all_color_values <- all_color_values[str_trim(all_color_values) != ""]
  all_color_values <- c(all_color_values, "Unknown")   # reserve a slot for unknowns
  global_order_palette <- make_color_palette(all_color_values)

  tax_df <- tax_df %>%
    select(genome_file, genome_order = all_of(color_level),
           tax_label = all_of(opt$tax_level)) %>%
    mutate(
      genome_order = if_else(is.na(genome_order) | str_trim(genome_order) == "",
                             "Unknown", as.character(genome_order)),
      tax_label    = if_else(is.na(tax_label)    | str_trim(tax_label)    == "",
                             "Unknown", as.character(tax_label))
    )

  # Merge on genome (df) = genome_file (tax_df)
  df <- df %>%
    left_join(tax_df, by = c("genome" = "genome_file")) %>%
    mutate(
      genome_order = if_else(is.na(genome_order), "Unknown", genome_order),
      tax_label    = if_else(is.na(tax_label),    "Unknown", tax_label)
    )

} else {
  global_order_palette <- make_color_palette("Unknown")
  df <- df %>% mutate(genome_order = "Unknown", tax_label = "All")
}

############################################################
# Build IPS_acc -> IPS_desc lookup
############################################################

ips_lookup <- df %>%
  filter(!is.na(IPS_acc), IPS_acc != "-") %>%
  select(IPS_acc, IPS_desc) %>%
  distinct() %>%
  mutate(
    acc_tokens  = str_split(IPS_acc,  "; "),
    desc_tokens = str_split(IPS_desc, "; ")
  ) %>%
  filter(lengths(acc_tokens) == lengths(desc_tokens)) %>%
  unnest(cols = c(acc_tokens, desc_tokens)) %>%
  filter(acc_tokens != "-") %>%
  select(acc_token = acc_tokens, desc_token = desc_tokens) %>%
  distinct(acc_token, .keep_all = TRUE)

############################################################
# Representative label resolver (IPS description only)
############################################################

resolve_label <- function(cluster_lbl) {
  if (is.na(cluster_lbl) || cluster_lbl == "-") return(cluster_lbl)

  tokens <- str_trim(str_split(cluster_lbl, ";")[[1]])
  tokens <- tokens[tokens != "-" & tokens != ""]
  if (length(tokens) == 0) return(cluster_lbl)

  exact <- ips_lookup %>% filter(acc_token == cluster_lbl)
  if (nrow(exact) > 0) return(exact$desc_token[1])

  counts <- ips_lookup %>%
    filter(acc_token %in% tokens) %>%
    mutate(n = sapply(acc_token, function(a) sum(tokens == a))) %>%
    arrange(desc(n))
  if (nrow(counts) > 0) return(counts$desc_token[1])

  return(cluster_lbl)
}

############################################################
# Shared theme / palettes
############################################################

shared_theme <- theme_bw(base_size = 16) +
  theme(
    strip.text        = element_text(face = "bold", size = 14),
    strip.background  = element_rect(fill = "#e8e8e8", color = NA),
    axis.text.y       = element_text(size = 10),
    legend.text       = element_text(size = 10),
    legend.title      = element_text(size = 12),
    legend.spacing.y  = unit(4, "mm"),
    legend.key.height = unit(6, "mm")
  )

no_y_theme <- theme(
  axis.text.y  = element_blank(),
  axis.ticks.y = element_blank(),
  axis.title.y = element_blank()
)

synteny_palette <- scale_fill_gradientn(
  colours = c("#008080", "#40E0D0", "#FFFACD", "#FF69B4", "#C2185B"),
  limits  = c(0, 1),
  name    = "Synteny\nScore",
  labels  = percent_format(accuracy = 1)
)

synteny_color_palette <- scale_color_gradientn(
  colours = c("#008080", "#40E0D0", "#FFFACD", "#FF69B4", "#C2185B"),
  limits  = c(0, 1),
  guide   = "none"
)



############################################################
# Helper: build the 1x3 combined plot
#
#   Panel 1 - Stacked genome count bar (coloured by order)
#             Secondary x-axis shows proportion.
#             Y-axis ticks = wrapped IPS description.
#   Panel 2 - |Gene Offset| boxplot + dotplot (synteny gradient)
#   Panel 3 - |Distance (bp)| boxplot + dotplot (synteny gradient)
#
# Panels 2 & 3 suppress y labels (aligned via patchwork).
############################################################

build_plot <- function(all_clusters, df_subset, total_genomes_subset,
                       order_palette, rank_by, title_suffix, taxon_label) {

  # --- Rank and factor clusters ---
  top_table <- all_clusters %>%
    arrange(desc(.data[[rank_by]])) %>%
    slice_head(n = 20) %>%
    mutate(
      annotation_cluster = factor(annotation_cluster, levels = rev(annotation_cluster))
    )

  top_ids     <- levels(top_table$annotation_cluster)
  synteny_map <- setNames(top_table$synteny_score, top_table$annotation_cluster)

  # Y-axis label: wrapped IPS description only (no cluster ID)
  desc_map <- setNames(
    str_wrap(top_table$rep_desc, width = 40),
    top_table$annotation_cluster
  )

  n_cl <- nrow(top_table)

  # --- Stacked bar data: genomes per (cluster x order) ---
  stacked_data <- df_subset %>%
    filter(annotation_cluster %in% top_ids) %>%
    mutate(annotation_cluster = factor(annotation_cluster, levels = top_ids)) %>%
    group_by(annotation_cluster, genome_order) %>%
    summarise(count = n_distinct(genome), .groups = "drop")

  # Full grid so missing order/cluster combos get 0 (clean stacking)
  all_combinations <- expand.grid(
    annotation_cluster = factor(top_ids, levels = top_ids),
    genome_order       = names(order_palette),
    stringsAsFactors   = FALSE
  ) %>%
    left_join(stacked_data, by = c("annotation_cluster", "genome_order")) %>%
    mutate(count = if_else(is.na(count), 0L, count))

  # Orders that actually appear (for legend breaks)
  active_orders <- all_combinations %>%
    filter(count > 0) %>%
    pull(genome_order) %>%
    unique()

  # --- Raw rows for synteny panels ---
  raw <- df_subset %>%
    filter(annotation_cluster %in% top_ids) %>%
    mutate(
      annotation_cluster = factor(annotation_cluster, levels = top_ids),
      synteny_score      = synteny_map[as.character(annotation_cluster)]
    )

  # ---- Panel 1: Stacked genome count bar ----
  p_stack <- ggplot(all_combinations,
                    aes(y    = annotation_cluster,
                        x    = count,
                        fill = genome_order)) +
    geom_col(width = 0.7, position = "stack") +
    scale_fill_manual(
      values = order_palette,
      name   = str_to_title(color_level),
      breaks = active_orders
    ) +
    scale_y_discrete(labels = desc_map) +
    scale_x_continuous(
      name     = "Genome Count",
      labels   = comma_format(accuracy = 1),
      sec.axis = sec_axis(
        ~ . / total_genomes_subset,
        name   = "Genome Proportion",
        labels = percent_format(accuracy = 1)
      )
    ) +
    labs(y = NULL, title = "Genome Count / Proportion by Order") +
    shared_theme

  # ---- Panel 2: |Gene Offset| ----
  p_offset <- ggplot(raw,
                     aes(y     = annotation_cluster,
                         x     = gene_offset,
                         fill  = synteny_score,
                         color = synteny_score)) +
    geom_boxplot(alpha = 0.5, width = 0.7, outlier.shape = 19) +
    # geom_dotplot(binaxis  = "x", dotsize  = 0.06, binwidth = 1, alpha = 0.7) +
    synteny_palette +
    synteny_color_palette +
    scale_x_continuous(labels = comma_format(accuracy = 1)) +
    scale_y_discrete(labels = desc_map) +
    labs(y = NULL, x = "Gene Offset", title = "Gene Offset") +
    shared_theme +
    no_y_theme

  # ---- Panel 3: |Distance (bp)| ----
  p_dist <- ggplot(raw,
                   aes(y     = annotation_cluster,
                       x     = distance_bp,
                       fill  = synteny_score,
                       color = synteny_score)) +
    geom_boxplot(alpha = 0.5, width = 0.7, outlier.shape = 19) +
    geom_dotplot(binaxis  = "x", stackdir = "center",
                 dotsize  = 0.06, binwidth = 1, alpha = 0.7) +
    synteny_palette +
    synteny_color_palette +
    scale_x_continuous(labels = comma_format(accuracy = 1)) +
    scale_y_discrete(labels = desc_map) +
    labs(y = NULL, x = "Distance (bp)", title = "Distance (bp)") +
    shared_theme +
    no_y_theme

  # ---- Combine 1x3 ----
  plot_height <- max(8, n_cl * 0.65)

  plot_title <- if (taxon_label == "All") {
    paste0("Top 20 neighborhood clusters for ", protein_name, " (", title_suffix, ")")
  } else {
    paste0("Top 20 neighborhood clusters for ", protein_name,
           " | ", opt$tax_level, ": ", taxon_label, " (", title_suffix, ")")
  }

  p_combined <- (p_stack | p_offset | p_dist) +
    plot_annotation(
      title = plot_title,
      theme = theme(plot.title = element_text(face = "bold", size = 20))
    ) +
    plot_layout(widths = c(1.8, 1, 1), guides = "collect") &
    theme(legend.position = "right")

  list(plot = p_combined, height = plot_height)
}

############################################################
# Helper: derive output path for a given taxon
############################################################

taxon_output_path <- function(base_path, taxon_label) {
  if (taxon_label == "All") return(base_path)

  # Sanitise taxon name for use as a directory component
  safe_taxon <- str_replace_all(taxon_label, "[^A-Za-z0-9._-]", "_")

  filename <- paste0(safe_taxon, ".clusters.svg")

  dir.create(base_path, recursive = TRUE, showWarnings = FALSE)

  file.path(base_path, filename)
}

############################################################
# Helper: run full pipeline for one taxon slice
############################################################

run_for_taxon <- function(df_slice, taxon_label) {

  total_genomes_subset <- df_slice %>% pull(genome) %>% unique() %>% length()
  cat(sprintf("  Taxon '%s': %d genomes\n", taxon_label, total_genomes_subset))

  all_clusters <- df_slice %>%
    filter(!is.na(annotation_cluster)) %>%
    group_by(annotation_cluster, cluster_label) %>%
    summarise(
      genome_count  = n_distinct(genome),
      synteny_score = {
        offset_counts <- table(gene_offset)
        max(offset_counts) / n_distinct(genome)
      },
      .groups = "drop"
    ) %>%
    mutate(
      proportion = genome_count / total_genomes_subset,
      rep_desc   = sapply(cluster_label, resolve_label)
    )

  if (nrow(all_clusters) == 0) {
    cat(sprintf("  Skipping '%s': no clusters found.\n", taxon_label))
    return(invisible(NULL))
  }

  # Use the global palette built from the full genome_file before any merge losses
  order_palette <- global_order_palette

  out_path1 <- taxon_output_path(opt$output,  taxon_label)
  out_path2 <- taxon_output_path(opt$output2, taxon_label)

  out1 <- build_plot(all_clusters, df_slice, total_genomes_subset,
                     order_palette, rank_by = "genome_count",
                     title_suffix = "by Genome Count", taxon_label = taxon_label)
  ggsave(out_path1, out1$plot, width = 30, height = out1$height, dpi = 300)
  cat("  Plot saved:", out_path1, "\n")

  out2 <- build_plot(all_clusters, df_slice, total_genomes_subset,
                     order_palette, rank_by = "synteny_score",
                     title_suffix = "by Synteny Score", taxon_label = taxon_label)
  ggsave(out_path2, out2$plot, width = 30, height = out2$height, dpi = 300)
  cat("  Plot saved:", out_path2, "\n")
}

############################################################
# Main: split by taxon and loop
############################################################

taxa_list <- df %>%
  group_by(tax_label) %>%
  summarise(n_genomes = n_distinct(genome), .groups = "drop") %>%
  arrange(desc(n_genomes))

cat(sprintf("\nFound %d taxa at level '%s'\n",
            nrow(taxa_list),
            if (!is.null(opt$tax_level)) opt$tax_level else "All"))

cat(sprintf("Producing plots for all %d taxa...\n\n", nrow(taxa_list)))

for (tl in taxa_list$tax_label) {
  df_slice <- df %>% filter(tax_label == tl)
  run_for_taxon(df_slice, tl)
}

cat("\nDone.\n")