#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))

# -----------------------------
# Parse arguments
# -----------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop(paste(
    "Usage: Rscript coencoding_heatmap.R <csv> <genome_col> <protein_col> <outfile_abs> <outfile_cond>",
    "  --proteins <p1 p2 ...>",
    "  [--group_col <col>]",
    sep = "\n"
  ))
}

input_csv    <- args[1]
genome_col   <- args[2]
protein_col  <- args[3]
outfile_abs  <- args[4]   # absolute: P(A and B) / N
outfile_cond <- args[5]   # conditional: P(B | A) = P(A and B) / P(A)

# --proteins
proteins_idx <- which(args == "--proteins")
if (length(proteins_idx) == 0) stop("Must provide --proteins followed by protein names")
remaining    <- args[(proteins_idx + 1):length(args)]
flag_pos     <- which(startsWith(remaining, "--"))
protein_list <- if (length(flag_pos) == 0) remaining else remaining[seq_len(flag_pos[1] - 1)]
protein_list <- trimws(gsub("^['\"]|['\"]$", "", protein_list))

# --group_col (optional: facet by e.g. phylum / tax level)
group_idx <- which(args == "--group_col")
group_col <- if (length(group_idx) > 0) args[group_idx + 1] else NULL

# -----------------------------
# Load & filter
# -----------------------------

df <- read_csv(input_csv, show_col_types = FALSE)

required_cols <- c(genome_col, protein_col, group_col)
missing_cols  <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))

df <- df %>% filter(.data[[protein_col]] %in% protein_list)
cat("Rows after filtering:", nrow(df), "\n")

# -----------------------------
# Compute protein x protein co-encoding proportions
#
# cell[A, B] = #{genomes with both A and B} / #{total genomes in group}
# diagonal   = #{genomes with A} / #{total genomes}  (single-protein prevalence)
# Matrix is symmetric; both triangles are computed for plotting.
# -----------------------------

compute_coenc <- function(data, label) {
  # One row per genome x protein (ignore copy number)
  wide <- data %>%
    distinct(genome = .data[[genome_col]], protein = .data[[protein_col]]) %>%
    mutate(present = 1L) %>%
    pivot_wider(names_from = protein, values_from = present, values_fill = 0L)

  # Ensure every requested protein column exists
  for (p in protein_list) {
    if (!p %in% colnames(wide)) wide[[p]] <- 0L
  }

  n <- nrow(wide)

  expand_grid(protein_A = protein_list, protein_B = protein_list) %>%
    rowwise() %>%
    mutate(
      group        = label,
      n_genomes    = n,
      n_both       = sum(wide[[protein_A]] == 1L & wide[[protein_B]] == 1L),
      n_A          = sum(wide[[protein_A]] == 1L),
      prop_abs     = n_both / n,                              # P(A and B) / N
      prop_cond    = if_else(n_A == 0, NA_real_, n_both / n_A)  # P(B | A)
    ) %>%
    ungroup()
}

# -----------------------------
# Run — optionally per group
# -----------------------------

if (is.null(group_col)) {
  n_total  <- n_distinct(df[[genome_col]])
  coenc_df <- compute_coenc(df, paste0("All genomes (n=", n_total, ")"))
} else {
  df <- df %>% filter(!is.na(.data[[group_col]]))
  coenc_df <- df %>%
    group_by(group_val = .data[[group_col]]) %>%
    group_map(~ {
      n <- n_distinct(.x[[genome_col]])
      compute_coenc(.x, paste0(.y$group_val, " (n=", n, ")"))
    }) %>%
    bind_rows()
}

coenc_df <- coenc_df %>%
  mutate(
    protein_A = factor(protein_A, levels = protein_list),
    protein_B = factor(protein_B, levels = protein_list)
  )

cat("Absolute proportion range:",  round(range(coenc_df$prop_abs,  na.rm = TRUE), 3), "\n")
cat("Conditional proportion range:", round(range(coenc_df$prop_cond, na.rm = TRUE), 3), "\n")

# -----------------------------
# Cluster proteins by absolute co-encoding similarity (shared order for both plots)
# -----------------------------

mean_mat <- coenc_df %>%
  group_by(protein_A, protein_B) %>%
  summarise(prop_abs = mean(prop_abs, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = protein_B, values_from = prop_abs, values_fill = 0) %>%
  column_to_rownames("protein_A") %>%
  as.matrix()

ordered_proteins <- if (nrow(mean_mat) > 2) {
  rownames(mean_mat)[hclust(dist(mean_mat), method = "average")$order]
} else {
  protein_list
}

# Apply shared ordering and shared helpers
coenc_df <- coenc_df %>%
  mutate(
    protein_A   = factor(protein_A, levels = ordered_proteins),
    protein_B   = factor(protein_B, levels = ordered_proteins),
    is_diagonal = protein_A == protein_B
  )

abs_df  <- coenc_df %>% mutate(
  proportion = prop_abs,
  cell_label = if_else(proportion == 0, "", sprintf("%.0f%%", proportion * 100))
)

cond_df <- coenc_df %>% mutate(
  proportion = prop_cond,
  cell_label = if_else(is.na(proportion) | proportion == 0, "", sprintf("%.0f%%", proportion * 100))
)

# -----------------------------
# Plot helper — used for both absolute and conditional heatmaps
# -----------------------------

make_plot <- function(df, title, subtitle, caption, legend_title) {
  n_groups <- length(unique(df$group))
  n_prot   <- length(protein_list)

  ggplot(df, aes(x = protein_B, y = fct_rev(protein_A), fill = proportion)) +

    geom_tile(aes(color = is_diagonal), linewidth = 0.5) +

    geom_text(
      aes(label = cell_label,
          color = I(if_else(!is.na(proportion) & proportion > 0.55, "white", "grey25"))),
      size = if (n_prot <= 10) 3.2 else 2.4
    ) +

    scale_fill_gradientn(
      name     = legend_title,
      colors   = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
      limits   = c(0, 1),
      breaks   = c(0, 0.25, 0.5, 0.75, 1),
      labels   = scales::percent,
      na.value = "grey90"
    ) +

    scale_color_manual(values = c("TRUE" = "grey55", "FALSE" = "white"), guide = "none") +

    { if (n_groups > 1) facet_wrap(~group, ncol = min(3, ceiling(sqrt(n_groups)))) } +

    scale_x_discrete(position = "top") +
    scale_y_discrete() +
    coord_equal(clip = "off") +

    theme_minimal(base_size = 12) +
    theme(
      axis.title        = element_blank(),
      axis.text.x       = element_text(angle = 40, hjust = 0, size = 9),
      axis.text.y       = element_text(size = 9),
      strip.text        = element_text(size = 9, face = "bold"),
      panel.grid        = element_blank(),
      legend.position   = "right",
      legend.title      = element_text(size = 9),
      legend.text       = element_text(size = 8),
      legend.key.height = unit(1.2, "cm"),
      plot.title        = element_text(size = 14, face = "bold"),
      plot.subtitle     = element_text(size = 9, color = "grey45", margin = margin(b = 8)),
      plot.caption      = element_text(size = 7.5, color = "grey50", hjust = 0),
      plot.margin       = margin(t = 8, r = 8, b = 8, l = 8, unit = "mm")
    ) +

    labs(title = title, subtitle = subtitle, caption = caption)
}

# -----------------------------
# Build both plots
# -----------------------------

n_groups <- length(unique(coenc_df$group))
n_prot   <- length(protein_list)
group_caption <- if (!is.null(group_col)) paste0("grouped by ", group_col) else
                   paste0("n = ", n_distinct(df[[genome_col]]), " genomes total")

p_abs <- make_plot(
  df           = abs_df,
  title        = "Protein co-encoding heatmap — absolute",
  subtitle     = paste0(n_prot, " proteins  \u00b7  diagonal = single-protein prevalence  \u00b7  proteins ordered by hierarchical clustering"),
  caption      = paste0("Cell [A, B]: proportion of all genomes encoding both A and B  \u00b7  ", group_caption),
  legend_title = "P(A \u2229 B) / N"
)

p_cond <- make_plot(
  df           = cond_df,
  title        = "Protein co-encoding heatmap — conditional",
  subtitle     = paste0(n_prot, " proteins  \u00b7  row = conditioning protein  \u00b7  diagonal = 1 by definition"),
  caption      = paste0("Cell [A, B]: P(B | A) = proportion of genomes with A that also have B  \u00b7  ", group_caption),
  legend_title = "P(B | A)"
)

# -----------------------------
# Save both
# -----------------------------

dir.create(dirname(outfile_abs),  recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(outfile_cond), recursive = TRUE, showWarnings = FALSE)

tile_cm    <- if (n_prot <= 8) 1.4 else if (n_prot <= 15) 1.1 else 0.85
n_cols_f   <- if (n_groups > 1) min(3, ceiling(sqrt(n_groups))) else 1
n_rows_f   <- if (n_groups > 1) ceiling(n_groups / n_cols_f) else 1
panel_size <- n_prot * tile_cm

plot_width  <- max(8, n_cols_f * (panel_size + 2) + 4)
plot_height <- max(6, n_rows_f * (panel_size + 2.5) + 2)

ggsave(outfile_abs,  p_abs,  width = plot_width, height = plot_height, units = "cm", dpi = 300)
message("Saved: ", outfile_abs)

ggsave(outfile_cond, p_cond, width = plot_width, height = plot_height, units = "cm", dpi = 300)
message("Saved: ", outfile_cond)