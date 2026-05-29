#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(viridis)
  library(ggrepel)
})

# -----------------------------
# Parse arguments
# -----------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("Usage: Rscript copy_heatmap.R <csv> <tax_level> <genome_col> <protein_col> <outfile> --proteins <p1 p2 ...>")
}

input_csv   <- args[1]
tax_level   <- args[2]
genome_col  <- args[3]
protein_col <- args[4]
outfile     <- args[5]

proteins_idx <- which(args == "--proteins")

if (length(proteins_idx) == 0) {
  stop("Must provide --proteins followed by protein names")
}

protein_list <- args[(proteins_idx + 1):length(args)]
protein_list <- gsub("^['\"]|['\"]$", "", protein_list)
protein_list <- trimws(protein_list)

# -----------------------------
# Full taxonomic hierarchy (coarse -> fine)
# -----------------------------

TAX_HIERARCHY  <- c("phylum", "class", "order", "family", "genus")
tax_level_idx  <- match(tax_level, TAX_HIERARCHY)

if (is.na(tax_level_idx)) {
  stop(paste("tax_level must be one of:", paste(TAX_HIERARCHY, collapse = ", ")))
}

ancestor_ranks <- TAX_HIERARCHY[seq_len(tax_level_idx)]

# -----------------------------
# Load dataframe
# -----------------------------

df <- read_csv(input_csv, show_col_types = FALSE)

required_cols <- c(genome_col, protein_col, ancestor_ranks)
missing_cols  <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop(paste("Missing columns in CSV:", paste(missing_cols, collapse = ", ")))
}

df <- df %>%
  filter(.data[[protein_col]] %in% protein_list) %>%
  filter(!is.na(.data[[tax_level]]))

df[[protein_col]] <- factor(df[[protein_col]], levels = protein_list)

cat("Rows after filtering:", nrow(df), "\n")
cat("Unique proteins:", unique(df[[protein_col]]), "\n")
cat("Unique taxa:", unique(df[[tax_level]]), "\n")

# -----------------------------
# Hierarchical sort order
# -----------------------------

taxon_sort_df <- df %>%
  select(all_of(ancestor_ranks)) %>%
  distinct() %>%
  mutate(across(all_of(ancestor_ranks), ~ replace_na(as.character(.), ""))) %>%
  mutate(sort_key = do.call(paste, c(select(., all_of(ancestor_ranks)), sep = "|"))) %>%
  arrange(sort_key) %>%
  pull(.data[[tax_level]])

cat("Taxon order (hierarchical):\n")
cat(paste(taxon_sort_df, collapse = "\n"), "\n\n")

# -----------------------------
# Count copies per genome
# -----------------------------

copy_counts <- df %>%
  group_by(
    genome  = .data[[genome_col]],
    taxon   = .data[[tax_level]],
    protein = .data[[protein_col]]
  ) %>%
  summarise(copies = n(), .groups = "drop")

# -----------------------------
# Count genomes per taxon
# -----------------------------

genome_counts <- df %>%
  distinct(genome = .data[[genome_col]], taxon = .data[[tax_level]]) %>%
  count(taxon, name = "genomes")

all_taxa <- unique(df[[tax_level]])

# -----------------------------
# Ensure every protein x taxon combination appears
# -----------------------------

max_copy_df <- expand_grid(
  taxon   = all_taxa,
  protein = factor(protein_list, levels = protein_list)
) %>%
  left_join(
    copy_counts %>%
      group_by(taxon, protein) %>%
      summarise(max_copy = max(copies), .groups = "drop"),
    by = c("taxon", "protein")
  ) %>%
  mutate(max_copy = replace_na(max_copy, 0))

# -----------------------------
# Expand copy levels
# -----------------------------

prop_df <- max_copy_df %>%
  rowwise() %>%
  mutate(copy_level = list(if (max_copy > 0) 1:max_copy else NA_real_)) %>%
  unnest(copy_level, keep_empty = TRUE) %>%
  left_join(genome_counts, by = "taxon") %>%
  rowwise() %>%
  mutate(
    genomes_with_copy = if (is.na(copy_level) || max_copy == 0) {
      0L
    } else {
      sum(
        copy_counts$copies[
          copy_counts$taxon == taxon & copy_counts$protein == protein
        ] >= copy_level
      )
    },
    proportion = if (max_copy == 0 || is.na(copy_level)) 0 else genomes_with_copy / genomes
  ) %>%
  ungroup() %>%
  mutate(copy_level_plot = if_else(is.na(copy_level), 1, copy_level))

# -----------------------------
# Row labels - apply hierarchical order
# -----------------------------

prop_df <- prop_df %>%
  left_join(genome_counts, by = "taxon", suffix = c("", ".y")) %>%
  mutate(
    genomes     = coalesce(genomes, genomes.y),
    taxon_label = paste0(taxon, " (", genomes, ")")
  ) %>%
  select(-ends_with(".y"))

ordered_taxon_labels <- prop_df %>%
  distinct(taxon, taxon_label) %>%
  mutate(taxon = factor(taxon, levels = taxon_sort_df)) %>%
  arrange(taxon) %>%
  pull(taxon_label)

prop_df <- prop_df %>%
  mutate(taxon_label = factor(taxon_label, levels = rev(ordered_taxon_labels)))

prop_df$protein <- factor(prop_df$protein, levels = protein_list)

cat("Proportion summary:\n")
print(summary(prop_df$proportion))

# -----------------------------
# Rare protein detection
# -----------------------------

if (length(protein_list) > 1) {
  rare_proteins <- prop_df %>%
    group_by(protein) %>%
    summarise(max_prop = max(proportion), .groups = "drop") %>%
    filter(max_prop < 0.05) %>%
    pull(protein)
} else {
  rare_proteins <- character(0)
}

dot_df  <- prop_df %>% filter(!protein %in% rare_proteins)
text_df <- prop_df %>%
  filter(protein %in% rare_proteins, copy_level_plot == 1) %>%
  mutate(label = as.character(genomes_with_copy))

# -----------------------------
# Dot spacing
# -----------------------------

DOT_STEP <- 0.2

n_proteins        <- length(protein_list)
global_max_copies <- if (nrow(dot_df) > 0) max(dot_df$copy_level_plot) else 1

if (nrow(dot_df) > 0) {
  max_copies_per_protein <- dot_df %>%
    group_by(protein) %>%
    summarise(n_levels = max(copy_level_plot), .groups = "drop")

  dot_df <- dot_df %>%
    left_join(max_copies_per_protein, by = "protein") %>%
    mutate(
      protein_idx = match(as.character(protein), protein_list),
      x_offset = if_else(
        n_levels == 1,
        0,
        (copy_level_plot - 1) * DOT_STEP - (n_levels - 1) * DOT_STEP / 3
      ),
      x_pos = protein_idx + x_offset
    )
} else {
  dot_df <- dot_df %>%
    mutate(
      n_levels    = integer(0),
      protein_idx = integer(0),
      x_offset    = numeric(0),
      x_pos       = numeric(0)
    )
}

text_df <- text_df %>%
  mutate(
    protein_idx = match(as.character(protein), protein_list),
    x_pos = protein_idx
  )

dot_perfect  <- dot_df %>% filter(proportion == 1.0)
dot_gradient <- dot_df %>% filter(proportion > 0 & proportion < 1.0)
dot_empty    <- dot_df %>% filter(proportion == 0, copy_level_plot == 1)

x_expand_add <- (global_max_copies - 1) * DOT_STEP / 3 + 0.2

# -----------------------------
# Alternating row bands for table feel
# -----------------------------

n_taxa     <- length(levels(prop_df$taxon_label))
row_levels <- levels(prop_df$taxon_label)

row_bands <- tibble(
  taxon_label = factor(row_levels, levels = row_levels),
  row_num     = seq_along(row_levels),
  shade       = row_num %% 2 == 0
) %>%
  filter(shade)

x_lo <- 0.5
x_hi <- n_proteins + x_expand_add

# -----------------------------
# Constants
# — PALETTE_HIGH matches proportion == 1.0 after direction = -1
#   (viridis flipped: yellow=0, purple=1)
# -----------------------------

PALETTE_HIGH <- "#440154"   # viridis purple = top of flipped scale
BAND_FILL    <- "#F2EFE9"
EMPTY_COL    <- "#CCCCCC"
HEADER_LINE  <- "#2A4858"
TEXT_MAIN    <- "#1C2B33"
TEXT_MUTED   <- "#7A8C96"

# -----------------------------
# Custom theme
# -----------------------------

theme_table <- function(base_size = 15) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),

      # ── No vertical grid lines ──
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),

      axis.title     = element_blank(),
      axis.ticks     = element_blank(),
      axis.text.x.top = element_text(
        angle  = 0, hjust = 0.5, vjust = -0.5,
        colour = TEXT_MAIN,
        face   = "bold",
        size   = base_size,
        margin = margin(b = 4)
      ),
      axis.text.y = element_text(
        hjust  = 1,
        colour = TEXT_MAIN,
        size   = base_size,
        margin = margin(r = 6)
      ),
      legend.position    = "right",
      legend.direction   = "vertical",
      legend.title       = element_text(colour = TEXT_MAIN, face = "bold", size = base_size*0.82),
      legend.text        = element_text(colour = TEXT_MUTED, size  = base_size*0.72),
      legend.key.height  = unit(0.9, "lines"),
      legend.key.width   = unit(0.55, "lines"),
      legend.margin      = margin(l = 8),
      plot.title = element_text(
        colour = HEADER_LINE,
        face   = "bold",
        size   = base_size * 1.15,
        margin = margin(b = 12)
      ),
      plot.subtitle = element_text(
        colour = TEXT_MUTED,
        size   = base_size * 0.85,
        margin = margin(b = 10)
      ),
      plot.margin = margin(t = 8, r = 22, b = 14, l = 10, unit = "mm")
    )
}

# -----------------------------
# Build plot
# -----------------------------

p <- ggplot(prop_df, aes(y = taxon_label)) +

  # Absent dots
  geom_point(
    data   = dot_empty,
    aes(x  = x_pos),
    shape  = 21, size = 5,
    fill   = NA, colour = EMPTY_COL,
    stroke = 0.6
  ) +

  # Gradient dots (0 < proportion < 1)
  geom_point(
    data   = dot_gradient,
    aes(x  = x_pos, fill = proportion),
    shape  = 21, size = 5,
    colour = "white", stroke = 0.25
  ) +

  # Perfect dots (proportion == 1) — star shape
  geom_point(
    data   = dot_perfect,
    aes(x  = x_pos),
    shape  = 24,          # 5-pointed star (open, native ggplot)
    size   = 5,
    colour = "white", 
    fill = PALETTE_HIGH,
    stroke = 1.0
  ) +

  # Rare-protein count labels
  geom_text(
    data  = text_df,
    aes(x = x_pos, label = label),
    size  = 3.2, colour = TEXT_MUTED,
    fontface = "italic"
  ) +

  # ── Viridis fill scale, flipped ──
  scale_fill_viridis_c(
    name      = "Genome proportion",
    option    = "cividis",
    direction = -1,          # yellow = 0, purple = 1
    limits    = c(0, 1),
    na.value  = "grey80",
    guide     = guide_colorbar(barheight = unit(5, "lines"))
  ) +

  scale_x_continuous(
    position = "top",
    breaks   = seq_along(protein_list),
    labels   = protein_list,
    expand   = expansion(add = x_expand_add)
  ) +

  geom_hline(yintercept = n_taxa + 0.5, colour = HEADER_LINE, linewidth = 0.9) +
  geom_hline(yintercept = 0.5,          colour = "#B0A898",   linewidth = 0.5) +

  coord_cartesian(clip = "off") +

  theme_table(base_size = 13) +

  labs(
    title    = paste("Protein copy distribution across", tax_level),
   )

# -----------------------------
# Save
# -----------------------------

dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)

col_width   <- 0.75 + (global_max_copies - 1) * DOT_STEP
plot_width  <- max(8, 2.5 + col_width * n_proteins)
plot_height <- max(4, 1.5 + 0.45 * length(unique(prop_df$taxon_label)))

ggsave(outfile, p, width = plot_width, height = plot_height, dpi = 300)
message("Saved: ", outfile)