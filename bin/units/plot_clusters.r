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
  make_option(c("--input"),   type="character"),
  make_option(c("--output"),  type="character"),
  make_option(c("--output2"), type="character"),
  make_option(c("--protein"), type="character")
)

opt <- parse_args(OptionParser(option_list = option_list))

if(is.null(opt$input))   stop("ERROR: --input not provided")
if(is.null(opt$output))  stop("ERROR: --output not provided")
if(is.null(opt$output2)) stop("ERROR: --output2 not provided")
if(is.null(opt$protein)) stop("ERROR: --protein not provided")

protein_name <- opt$protein

############################################################
# Load data
############################################################

df <- read_csv(opt$input, show_col_types = FALSE)

############################################################
# Validate columns
############################################################

required_cols <- c("annotation_cluster", "cluster_label", "genome", "IPS_acc", "IPS_desc",
                   "gene_offset", "distance_bp")

missing <- setdiff(required_cols, colnames(df))

if(length(missing) > 0){
  stop(paste("Missing columns:", paste(missing, collapse=", ")))
}

############################################################
# Build IPS_acc -> IPS_desc lookup (one desc per acc token)
############################################################

ips_lookup <- df %>%
  filter(!is.na(IPS_acc), IPS_acc != "-") %>%
  select(IPS_acc, IPS_desc) %>%
  distinct() %>%
  # Each row may have multiple acc tokens separated by "; "
  # Explode to one token per row and pair with its corresponding desc token
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
# Representative label resolver
# For a given cluster_label ("; "-separated IPS_acc string),
# return the IPS_desc of the best matching IPS_acc token.
############################################################

resolve_label <- function(cluster_lbl) {
  if (is.na(cluster_lbl) || cluster_lbl == "-") return(cluster_lbl)

  tokens <- str_trim(str_split(cluster_lbl, ";")[[1]])
  tokens <- tokens[tokens != "-" & tokens != ""]

  if (length(tokens) == 0) return(cluster_lbl)

  # 1. Exact match: cluster_label == a single IPS_acc token
  exact <- ips_lookup %>% filter(acc_token == cluster_lbl)
  if (nrow(exact) > 0) return(exact$desc_token[1])

  # 2. Best overlap: acc_token that appears most in the token list
  counts <- ips_lookup %>%
    filter(acc_token %in% tokens) %>%
    mutate(n = sapply(acc_token, function(a) sum(tokens == a))) %>%
    arrange(desc(n))

  if (nrow(counts) > 0) return(counts$desc_token[1])

  # 3. Fallback
  return(cluster_lbl)
}

############################################################
# Compute total genomes
############################################################

total_genomes <- df %>%
  pull(genome) %>%
  unique() %>%
  length()

############################################################
# Helper: build the 1x4 combined plot for a given top-20 slice
#   rank_by : "genome_count" | "synteny_score"
#   bar_mode: "cluster_colors" | "synteny"
############################################################

build_plot <- function(all_clusters, rank_by, title_suffix) {

  # --- Slice top 20 by chosen metric ---
  top_table <- all_clusters %>%
    arrange(desc(.data[[rank_by]])) %>%
    slice_head(n = 20) %>%
    mutate(
      annotation_cluster = factor(
        annotation_cluster,
        levels = rev(annotation_cluster)   # highest at top
      )
    )

  top_ids           <- levels(top_table$annotation_cluster)
  cluster_label_map <- setNames(top_table$legend_label,  top_table$annotation_cluster)
  synteny_map       <- setNames(top_table$synteny_score, top_table$annotation_cluster)

  n_cl <- nrow(top_table)

  # --- Cluster colors for bar panels (Set3 palette) ---
  cluster_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_cl)
  names(cluster_colors) <- top_ids

  # --- Pivot for bar panels ---
  freq_long <- top_table %>%
    tidyr::pivot_longer(
      cols      = c(genome_count, proportion),
      names_to  = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric = recode(metric,
        genome_count = "Genome Count",
        proportion   = "Genome Proportion"
      ),
      metric = factor(metric, levels = c("Genome Count", "Genome Proportion"))
    )

  # --- Raw rows for top 20 (with synteny score joined for panels 3+4) ---
  raw <- df %>%
    filter(annotation_cluster %in% top_ids) %>%
    mutate(
      annotation_cluster = factor(annotation_cluster, levels = top_ids),
      synteny_score      = synteny_map[as.character(annotation_cluster)]
    )

  # --- Panel 1: Genome Count (Set3 colors, full legend) ---
  p_count <- ggplot(
      freq_long %>% filter(metric == "Genome Count"),
      aes(y = annotation_cluster, x = value, fill = annotation_cluster)
    ) +
    geom_col(width = 0.7) +
    scale_fill_manual(
      values = cluster_colors,
      labels = cluster_label_map,
      name   = "Cluster"
    ) +
    scale_x_continuous(labels = comma_format(accuracy = 1)) +
    labs(y = "Cluster ID", x = "Genome Count", title = "Genome Count") +
    shared_theme

  # --- Panel 2: Genome Proportion (Set3 colors, secondary axis = genome count) ---
  p_prop <- ggplot(
      freq_long %>% filter(metric == "Genome Proportion"),
      aes(y = annotation_cluster, x = value, fill = annotation_cluster)
    ) +
    geom_col(width = 0.7) +
    scale_fill_manual(
      values = cluster_colors,
      labels = cluster_label_map,
      guide  = "none"
    ) +
    scale_x_continuous(
      labels   = percent_format(accuracy = 1),
      sec.axis = sec_axis(
        ~ . * total_genomes,
        name   = "Genome Count",
        labels = comma_format(accuracy = 1)
      )
    ) +
    labs(y = NULL, x = "Genome Proportion", title = "Genome Proportion") +
    shared_theme +
    no_y_theme

  # --- Panel 3: |Gene Offset| boxplot + dotplot (synteny gradient) ---
  p_offset <- ggplot(raw,
                     aes(y     = annotation_cluster,
                         x     = (gene_offset),
                         fill  = synteny_score,
                         color = synteny_score)) +
    geom_boxplot(alpha = 0.5, width = 0.7, outlier.shape = 19) +
    geom_dotplot(binaxis = "x", stackdir = "center",
                 dotsize = 0.06, binwidth = 1, alpha = 0.7) +
    synteny_palette +
    synteny_color_palette +
    scale_x_continuous(labels = comma_format(accuracy = 1)) +
    labs(y = NULL, x = "|Gene Offset|", title = "|Gene Offset|") +
    shared_theme +
    no_y_theme

  # --- Panel 4: |Distance (bp)| boxplot + dotplot (synteny gradient) ---
  p_dist <- ggplot(raw,
                   aes(y     = annotation_cluster,
                       x     = (distance_bp),
                       fill  = synteny_score,
                       color = synteny_score)) +
    geom_boxplot(alpha = 0.5, width = 0.7, outlier.shape = 19) +
    geom_dotplot(binaxis = "x", stackdir = "center",
                 dotsize = 0.06, binwidth = 1, alpha = 0.7) +
    synteny_palette +
    synteny_color_palette +
    scale_x_continuous(labels = comma_format(accuracy = 1)) +
    labs(y = NULL, x = "|Distance (bp)|", title = "|Distance (bp)|") +
    shared_theme +
    no_y_theme

  # --- Combine 1x4 ---
  plot_height <- max(8, n_cl * 0.6)

  p_combined <- (p_count | p_prop | p_offset | p_dist) +
    plot_annotation(
      title = paste0("Top 20 neighborhood clusters for ", protein_name,
                     " (", title_suffix, ")"),
      theme = theme(plot.title = element_text(face = "bold", size = 20))
    ) +
    plot_layout(widths = c(1.4, 1, 1, 1), guides = "collect") &
    theme(legend.position = "right")

  list(plot = p_combined, height = plot_height)
}

############################################################
# Compute full cluster table (all clusters, not sliced yet)
############################################################

all_clusters <- df %>%
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
    proportion   = genome_count / total_genomes,
    rep_desc     = sapply(cluster_label, resolve_label),
    legend_label = paste0(annotation_cluster, " - ", str_wrap(rep_desc, width = 50))
  )

############################################################
# Shared theme / palettes (used inside build_plot)
############################################################

shared_theme <- theme_bw(base_size = 16) +
  theme(
    strip.text        = element_text(face = "bold", size = 14),
    strip.background  = element_rect(fill = "#e8e8e8", color = NA),
    axis.text.y       = element_text(size = 11),
    legend.text       = element_text(size = 10),
    legend.title      = element_text(size = 12),
    legend.spacing.y  = unit(6, "mm"),
    legend.key.height = unit(8, "mm")
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
# Build & save: ranked by genome count
############################################################

out1 <- build_plot(all_clusters, rank_by = "genome_count",  title_suffix = "by Genome Count")
ggsave(opt$output,  out1$plot, width = 36, height = out1$height, dpi = 300)
cat("Plot saved:", opt$output, "\n")

############################################################
# Build & save: ranked by synteny score
############################################################

out2 <- build_plot(all_clusters, rank_by = "synteny_score", title_suffix = "by Synteny Score")
ggsave(opt$output2, out2$plot, width = 36, height = out2$height, dpi = 300)
cat("Plot saved:", opt$output2, "\n")