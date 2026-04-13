#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))

# -----------------------------
# Parse arguments
# -----------------------------

args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 5){
  stop("Usage: Rscript copy_heatmap.R <csv> <tax_level> <genome_col> <protein_col> <outfile> --proteins <p1 p2 ...>")
}

input_csv   <- args[1]
tax_level   <- args[2]
genome_col  <- args[3]
protein_col <- args[4]
outfile     <- args[5]

proteins_idx <- which(args == "--proteins")

if(length(proteins_idx) == 0){
  stop("Must provide --proteins followed by protein names")
}

protein_list <- args[(proteins_idx + 1):length(args)]

protein_list <- gsub("^['\"]|['\"]$", "", protein_list)
protein_list <- trimws(protein_list)

# -----------------------------
# Load dataframe
# -----------------------------

df <- read_csv(input_csv, show_col_types = FALSE)

df <- df %>%
  filter(.data[[protein_col]] %in% protein_list) %>%
  filter(!is.na(.data[[tax_level]]))

df[[protein_col]] <- factor(df[[protein_col]], levels = protein_list)

cat("Rows after filtering:", nrow(df), "\n")
cat("Unique proteins:", unique(df[[protein_col]]), "\n")
cat("Unique taxa:", unique(df[[tax_level]]), "\n")

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

cat("Copy counts rows:", nrow(copy_counts), "\n")

# -----------------------------
# Count genomes per taxon
# -----------------------------

genome_counts <- df %>%
  distinct(
    genome = .data[[genome_col]],
    taxon  = .data[[tax_level]]
  ) %>%
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
          copy_counts$taxon   == taxon &
          copy_counts$protein == protein
        ] >= copy_level
      )
    },
    proportion = if (max_copy == 0 || is.na(copy_level)) 0 else genomes_with_copy / genomes
  ) %>%
  ungroup() %>%
  mutate(copy_level_plot = if_else(is.na(copy_level), 1, copy_level))

# -----------------------------
# Row labels
# -----------------------------

prop_df <- prop_df %>%
  left_join(genome_counts, by = "taxon", suffix = c("", ".y")) %>%
  mutate(
    genomes     = coalesce(genomes, genomes.y),
    taxon_label = paste0(taxon, " (n=", genomes, ")")
  ) %>%
  select(-ends_with(".y"))

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
#
# DOT_STEP is fixed at 0.3 units — always comfortable to read.
# Plot width grows to accommodate however many copy levels exist,
# rather than squishing dots to fit a fixed canvas.
# -----------------------------

DOT_STEP <- 0.2   # fixed spacing between copy-level dots (data units)

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
      # centre the spread around the protein position
      x_offset = if_else(
        n_levels == 1,
        0,
        (copy_level_plot - 1) * DOT_STEP - (n_levels - 1) * DOT_STEP / 2
      ),
      x_pos = protein_idx + x_offset
    )

} else {
  dot_df <- dot_df %>%
    mutate(n_levels = integer(0), protein_idx = integer(0),
           x_offset = numeric(0), x_pos = numeric(0))
}

text_df <- text_df %>%
  mutate(
    protein_idx = match(as.character(protein), protein_list),
    x_pos = protein_idx
  )

dot_perfect  <- dot_df %>% filter(proportion == 1.0)
dot_gradient <- dot_df %>% filter(proportion > 0 & proportion < 1.0)
dot_empty    <- dot_df %>% filter(proportion == 0, copy_level_plot == 1)

# Padding: half a cell + half the max spread width on each side
x_expand_add <- 0.5 + (global_max_copies - 1) * DOT_STEP / 2 + 0.5

# -----------------------------
# Plot
# -----------------------------

p <- ggplot(prop_df, aes(y = taxon_label)) +

  geom_point(
    data  = dot_empty,
    aes(x = x_pos),
    shape = 21, size = 3, fill = NA, color = "grey70"
  ) +

  geom_point(
    data = dot_gradient,
    aes(x = x_pos, color = proportion),
    size = 3
  ) +

  geom_point(
    data   = dot_perfect,
    aes(x  = x_pos),
    shape  = 21, fill = "#FF4D6D", color = "black",
    size   = 3, stroke = 1.2
  ) +

  geom_text(
    data  = text_df,
    aes(x = x_pos, label = label),
    size  = 3.5, color = "grey40"
  ) +

  scale_color_gradient(
    name     = "Genome proportion",
    low      = "#00C9A7",
    high     = "#FF4D6D",
    limits   = c(0, 1),
    na.value = "grey80"
  ) +

  scale_x_continuous(
    position = "top",
    breaks   = seq_along(protein_list),
    labels   = protein_list,
    expand   = expansion(add = x_expand_add)
  ) +

  coord_cartesian(clip = "off") +

  theme_classic(base_size = 14) +

  theme(
    axis.title  = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 0),
    axis.text.y = element_text(hjust = 1),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    plot.margin = margin(t = 10, r = 20, b = 15, l = 10, unit = "mm")
  ) +

  labs(
    title   = paste("Protein copy distribution across", tax_level),
    caption = str_wrap(
      "Outlined dot = 100% of genomes | Gradient = partial presence | Open circle = absent | Grey numbers = count (<5% of genomes)",
      width = 80
    )
  )

# -----------------------------
# Save — width grows with both protein count AND max copy number
# -----------------------------

dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)

# Each protein column needs room for its copy-level spread:
# base width per protein = 1.5 inches + (max_copies - 1) * DOT_STEP inches
col_width   <- 1.5 + (global_max_copies - 1) * DOT_STEP
plot_width  <- max(8, 2.5 + col_width * n_proteins)
plot_height <- max(4, 1.5 + 0.45 * length(unique(prop_df$taxon_label)))

ggsave(outfile, p, width = plot_width, height = plot_height, dpi = 300)
message("Saved: ", outfile)