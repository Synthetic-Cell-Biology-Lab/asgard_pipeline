#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(scales)
  library(colorspace)
  library(ggtext)       # for element_markdown / richtext labels
  library(systemfonts)  # better font resolution
  library(tibble)
})

# ==============================
# ARGUMENTS
# ==============================
option_list <- list(
  make_option(c("-i", "--input"),  type="character", help="Input CSV"),
  make_option(c("-o", "--output"), type="character", help="Output SVG/PDF/PNG")
)

opt         <- parse_args(OptionParser(option_list=option_list))
INPUT_CSV   <- opt$input
OUTPUT_FILE <- opt$output

# ==============================
# LOAD & DEDUPLICATE
# ==============================
df <- read_csv(INPUT_CSV, show_col_types = FALSE)

if (!all(c("class", "order", "genome_file") %in% colnames(df))) {
  stop("CSV must contain 'class', 'order', and 'genome_file' columns")
}

df <- df %>%
  group_by(genome_file) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    class = ifelse(is.na(class) | class == "", "Unknown", class),
    order = ifelse(is.na(order) | order == "", "Unknown", order)
  )

# ==============================
# COUNT
# ==============================
counts <- df %>%
  count(class, order, name = "count")

class_totals <- counts %>%
  group_by(class) %>%
  summarise(total = sum(count), .groups = "drop") %>%
  arrange(desc(total))

# highest total at top → reverse for ggplot y-axis
counts$class <- factor(counts$class, levels = rev(class_totals$class))

# ==============================
# PUBLICATION COLOR PALETTE
# Hand-picked perceptually-distinct hues; desaturated tints for shading
# Uses HCL space for uniform luminance across classes
# ==============================
class_levels <- levels(counts$class)
n_classes    <- length(class_levels)

# Evenly-spaced HCL hues — perceptually uniform, print-safe
hcl_hues <- seq(15, 375 - (360 / n_classes), length.out = n_classes)
base_colors <- hcl(h = hcl_hues, c = 65, l = 58)
names(base_colors) <- class_levels

color_map <- c()

for (cls in class_levels) {
  sub <- counts %>% filter(class == cls) %>% arrange(order)
  n   <- nrow(sub)

  if (n == 1) {
    shades <- base_colors[[cls]]
  } else {
    # Ramp in HCL space: low chroma/high luminance → full color
    # avoids the muddy mid-tones of RGB interpolation
    start <- lighten(base_colors[[cls]], amount = 0.55, space = "HLS")
    shades <- colorRampPalette(c(start, base_colors[[cls]]))(n)
  }

  keys          <- paste(cls, sub$order, sep = "||")
  names(shades) <- keys
  color_map     <- c(color_map, shades)
}

counts <- counts %>%
  mutate(key = paste(class, order, sep = "||"))

counts$key <- factor(counts$key, levels = names(color_map))

# ==============================
# LEGEND LABELS & ORDER
# ==============================
legend_order <- counts %>%
  arrange(class, order) %>%
  pull(key) %>%
  as.character()

legend_labels <- counts %>%
  arrange(class, order) %>%
  mutate(label = paste0(class, " \u2014 ", order, "  (n\u202f=\u202f", count, ")")) %>%
  select(key, label) %>%
  deframe()

# ==============================
# BAR-END TOTAL LABELS
# ==============================
bar_labels <- counts %>%
  group_by(class) %>%
  summarise(total = sum(count), .groups = "drop")

# ==============================
# THEME CONSTANTS
# ==============================
FONT_MAIN  <- "Helvetica Neue"   # falls back gracefully
FONT_TITLE <- "Georgia"
BG_COLOR   <- "#FAFAF8"          # warm off-white — easier on eyes than pure white
GRID_COLOR <- "#E8E5DF"
TEXT_DARK  <- "#1A1A1A"
TEXT_MID   <- "#555550"
TEXT_LIGHT <- "#888880"

# ==============================
# PLOT
# ==============================
p <- ggplot(counts, aes(x = count, y = class, fill = key)) +

  # subtle reference lines before bars so bars sit on top
  geom_vline(
    xintercept = seq(0, max(class_totals$total) * 1.1, by = 50),
    color      = GRID_COLOR,
    linewidth  = 0.4
  ) +

  geom_bar(
    stat      = "identity",
    position  = "stack",
    color     = BG_COLOR,          # segment dividers match background — clean
    linewidth = 0.35,
    width     = 0.68
  ) +

  # bold total at bar end
  geom_text(
    data        = bar_labels,
    aes(x = total, y = class, label = total),
    inherit.aes = FALSE,
    hjust       = -0.35,
    vjust       = 0.38,
    size        = 3.2,
    fontface    = "bold",
    color       = TEXT_DARK,
    family      = FONT_MAIN
  ) +

  scale_fill_manual(
    values = color_map,
    labels = legend_labels,
    breaks = legend_order,
    name   = "Orders grouped by Class",
    guide  = guide_legend(
      ncol          = 1,
      reverse       = TRUE,
      keywidth      = unit(0.9, "lines"),
      keyheight     = unit(0.75, "lines"),
      title.hjust   = 0,
      label.hjust   = 0
    )
  ) +

  scale_x_continuous(
    expand = expansion(mult = c(0, 0.10)),
    labels = label_number(accuracy = 1)
  ) +

  labs(
    title    = "Genome Distribution by Class",
    subtitle = "Orders stacked within each Class; bars sorted by total genome count",
    x        = "Number of Genomes",
    y        = NULL,
    caption  = "One genome per unique genome_file entry."
  ) +

  theme_minimal(base_size = 11, base_family = FONT_MAIN) +
  theme(
    # ── Canvas ───────────────────────────────────────────────
    plot.background  = element_rect(fill = BG_COLOR, color = NA),
    panel.background = element_rect(fill = BG_COLOR, color = NA),
    plot.margin      = margin(20, 18, 14, 14),

    # ── Title / subtitle / caption ───────────────────────────
    plot.title = element_text(
      family   = FONT_TITLE,
      face     = "bold",
      size     = 15,
      color    = TEXT_DARK,
      margin   = margin(b = 3)
    ),
    plot.subtitle = element_text(
      size   = 9.5,
      color  = TEXT_MID,
      margin = margin(b = 14)
    ),
    plot.caption = element_text(
      size   = 8,
      color  = TEXT_LIGHT,
      hjust  = 0,
      margin = margin(t = 10)
    ),

    # ── Axes ─────────────────────────────────────────────────
    axis.title.x = element_text(
      size   = 9.5,
      color  = TEXT_MID,
      margin = margin(t = 8)
    ),
    axis.text.y  = element_text(
      size   = 9.5,
      color  = TEXT_DARK,
      face   = "italic",          # italics for taxonomic names
      hjust  = 1
    ),
    axis.text.x  = element_text(size = 9, color = TEXT_MID),
    axis.ticks   = element_blank(),

    # ── Grid ─────────────────────────────────────────────────
    panel.grid       = element_blank(),  # drawn manually via geom_vline above
    panel.border     = element_blank(),

    # ── Legend ───────────────────────────────────────────────
    legend.position      = "right",
    legend.justification = "top",
    legend.title         = element_text(
      size     = 8.5,
      face     = "bold",
      color    = TEXT_DARK,
      margin   = margin(b = 6)
    ),
    legend.text  = element_text(size = 7.8, color = TEXT_MID),
    legend.key   = element_rect(fill = NA, color = NA),
    legend.background = element_rect(fill = NA, color = NA),
    legend.spacing.y  = unit(0.18, "lines")
  )

# ==============================
# SAVE
# ==============================
ext        <- tolower(tools::file_ext(OUTPUT_FILE))
n_cls      <- length(unique(counts$class))
plot_h     <- max(6.5, n_cls * 0.72 + 2.5)
plot_w     <- 14

if (ext == "pdf") {
  ggsave(OUTPUT_FILE, plot = p, width = plot_w, height = plot_h,
         device = cairo_pdf, dpi = 300)
} else if (ext == "png") {
  ggsave(OUTPUT_FILE, plot = p, width = plot_w, height = plot_h,
         device = ragg::agg_png, dpi = 400, bg = BG_COLOR)
} else {
  ggsave(OUTPUT_FILE, plot = p, width = plot_w, height = plot_h,
         device = "svg", bg = BG_COLOR)
}

cat("[✓] Saved →", OUTPUT_FILE, "\n")