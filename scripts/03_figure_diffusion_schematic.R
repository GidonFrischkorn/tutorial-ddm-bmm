#' Diffusion Process Schematic
#'
#' Creates a conceptual diagram of the Diffusion Decision Model (DDM)
#' showing evidence accumulation between two response boundaries.

pacman::p_load(here, tidyverse, patchwork)
source(here("functions", "clean_plot.R"))

set.seed(42)

# --- DDM parameters for the illustration ---
drift   <- 2.0
bound   <- 1.5
ndt     <- 0.30
zr      <- 0.5
dt      <- 0.001
n_steps <- 600

# --- Simulate a single diffusion path ---
noise <- rnorm(n_steps, mean = drift * dt, sd = sqrt(dt))
evidence <- cumsum(noise) + zr * bound
time <- seq(0, n_steps * dt, by = dt)[1:n_steps] + ndt

# Truncate at boundary crossing
hit_upper <- which(evidence >= bound)[1]
hit_lower <- which(evidence <= 0)[1]
hit <- min(hit_upper, hit_lower, na.rm = TRUE)
evidence <- evidence[1:hit]
time <- time[1:hit]

# --- Second path (error, for illustration) ---
set.seed(7)
noise2 <- rnorm(n_steps, mean = 0.3 * dt, sd = sqrt(dt))
evidence2 <- cumsum(noise2) + zr * bound
time2 <- seq(0, n_steps * dt, by = dt)[1:n_steps] + ndt

hit_upper2 <- which(evidence2 >= bound)[1]
hit_lower2 <- which(evidence2 <= 0)[1]
hit2 <- min(hit_upper2, hit_lower2, na.rm = TRUE)
evidence2 <- evidence2[1:hit2]
time2 <- time2[1:hit2]

# --- Build the figure ---
df_path1 <- tibble(time = time, evidence = evidence, path = "Correct")
df_path2 <- tibble(time = time2, evidence = evidence2, path = "Error")
df_paths <- bind_rows(df_path1, df_path2)

# Annotations
x_max <- max(df_paths$time) + 0.05

p_schematic <- ggplot() +
  # Boundary lines
  geom_hline(yintercept = bound, linewidth = 1, color = "grey30") +
  geom_hline(yintercept = 0, linewidth = 1, color = "grey30") +

  # Non-decision time shading
  annotate("rect", xmin = 0, xmax = ndt, ymin = -0.15, ymax = bound + 0.15,
           fill = "grey90", alpha = 0.6) +

  # Starting point
  geom_point(aes(x = ndt, y = zr * bound),
             size = 3, color = "grey30") +

  # Evidence paths
  geom_line(data = df_path1, aes(x = time, y = evidence),
            color = eam_palette[1], linewidth = 0.7, alpha = 0.9) +
  geom_line(data = df_path2, aes(x = time, y = evidence),
            color = eam_palette[4], linewidth = 0.7, alpha = 0.9) +

  # Boundary crossing markers
  geom_point(data = df_path1 |> slice_tail(n = 1),
             aes(x = time, y = evidence),
             size = 3, color = eam_palette[1]) +
  geom_point(data = df_path2 |> slice_tail(n = 1),
             aes(x = time, y = evidence),
             size = 3, color = eam_palette[4]) +

  # Parameter annotations
  # Bound (alpha)
  annotate("segment", x = 0.05, xend = 0.05, y = 0, yend = bound,
           arrow = arrow(ends = "both", length = unit(0.08, "inches")),
           color = "grey40", linewidth = 0.4) +
  annotate("text", x = 0.02, y = bound / 2, label = "italic(alpha)",
           parse = TRUE, hjust = 1, size = 4, color = "grey30") +

  # Starting point (zr)
  annotate("segment", x = ndt - 0.01, xend = ndt - 0.01,
           y = 0, yend = zr * bound,
           arrow = arrow(ends = "both", length = unit(0.06, "inches")),
           color = "grey40", linewidth = 0.4, linetype = "dotted") +
  annotate("text", x = ndt - 0.03, y = zr * bound / 2, label = "italic(zeta)",
           parse = TRUE, hjust = 1, size = 3.5, color = "grey30") +

  # Drift arrow
  annotate("segment", x = ndt + 0.05, xend = ndt + 0.15,
           y = zr * bound, yend = zr * bound + drift * 0.1 * 0.6,
           arrow = arrow(length = unit(0.12, "inches"), type = "closed"),
           color = "grey30", linewidth = 0.8) +
  annotate("text", x = ndt + 0.12, y = zr * bound + 0.17,
           label = "italic(delta)", parse = TRUE, size = 4, color = "grey30") +

  # Ndt bracket
  annotate("segment", x = 0, xend = ndt, y = -0.1, yend = -0.1,
           arrow = arrow(ends = "both", length = unit(0.06, "inches")),
           color = "grey40", linewidth = 0.4) +
  annotate("text", x = ndt / 2, y = -0.18, label = "italic(tau)",
           parse = TRUE, size = 4, color = "grey30") +

  # Boundary labels
  annotate("text", x = x_max, y = bound, label = "Upper boundary",
           hjust = 0, vjust = -0.5, size = 3.2, color = "grey40") +
  annotate("text", x = x_max, y = 0, label = "Lower boundary",
           hjust = 0, vjust = 1.5, size = 3.2, color = "grey40") +

  # Axis labels
  labs(x = "Time (s)", y = "Evidence") +

  # Scales
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  scale_y_continuous(breaks = c(0, bound), labels = c("0", expression(alpha))) +

  # Theme
  clean_plot() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank()
  )

ggsave(here("figures", "03_diffusion_schematic.pdf"), p_schematic,
       width = 6.5, height = 3.5, device = "pdf")

cat("Diffusion schematic saved.\n")
