# Plot theme and color-blind friendly palette for the tutorial figures
# Usage: ggplot(...) + clean_plot(); extra arguments are passed to theme()

###############################################################################!
# Color Palette ----------------------------------------------------------------
###############################################################################!

# Okabe-Ito color-blind friendly palette
eam_palette <- c(
  "#0072B2",
  "#E69F00",
  "#009E73",
  "#D55E00",
  "#56B4E9",
  "#CC79A7",
  "#F0E442",
  "#000000"
)

# Scales that apply the palette
scale_color_eam <- function(...) {
  scale_color_manual(values = eam_palette, ...)
}

scale_fill_eam <- function(...) {
  scale_fill_manual(values = eam_palette, ...)
}

###############################################################################!
# Plot Theme -------------------------------------------------------------------
###############################################################################!

# base_size 11 suits figures saved 6.5 in wide; ... is passed on to theme()
clean_plot <- function(base_size = 11, base_family = "", ...) {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      # no grid or panel border
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border     = element_blank(),

      # axis lines and ticks
      axis.line  = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.4),
      axis.text  = element_text(color = "black", size = rel(0.82)),

      # legend below the plot; override with theme() after clean_plot()
      legend.position   = "bottom",
      legend.key        = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),

      # facet strip labels
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text       = element_text(face = "bold"),

      # plot title and subtitle
      plot.title    = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(color = "grey30", hjust = 0),
      plot.margin   = margin(10, 10, 10, 10),

      ...
    )
}
