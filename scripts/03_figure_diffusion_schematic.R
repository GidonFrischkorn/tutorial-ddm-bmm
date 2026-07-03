#' Diffusion Process Schematic
#'
#' Creates a conceptual diagram of the Diffusion Decision Model (DDM):
#' several noisy evidence paths accumulate from the starting point toward two
#' response boundaries, with the resulting RT distributions attached to each
#' boundary (correct responses above, errors below). All four parameters are
#' annotated: drift (delta), boundary separation (alpha), non-decision time
#' (tau), and relative starting point (zeta_r).

pacman::p_load(here, bmm, tidyverse, patchwork)
source(here("functions", "clean_plot.R"))

set.seed(42)

# --- DDM parameters for the illustration ---
drift   <- 2.0
bound   <- 1.5
ndt     <- 0.30
zr      <- 0.5
dt      <- 0.001
n_steps <- 1200
x_lim   <- 1.25   # right edge of the plot (seconds)

# --- Simulate candidate diffusion paths and keep a small illustrative set ----
sim_paths <- map(1:300, function(i) {
  noise <- rnorm(n_steps, mean = drift * dt, sd = sqrt(dt))
  evidence <- cumsum(noise) + zr * bound
  hit <- suppressWarnings(
    min(which(evidence >= bound)[1], which(evidence <= 0)[1], na.rm = TRUE)
  )
  if (is.infinite(hit) || is.na(hit)) return(NULL)
  tibble(
    path_id  = i,
    time     = seq_len(hit) * dt + ndt,
    evidence = pmin(pmax(evidence[seq_len(hit)], 0), bound),
    outcome  = if (evidence[hit] >= bound) "upper" else "lower"
  )
}) |>
  compact() |>
  bind_rows()

# Keep paths that terminate within the plotting window and are neither
# too short nor too long for a readable display
path_info <- sim_paths |>
  summarise(rt = max(time), outcome = first(outcome), .by = path_id) |>
  filter(rt > 0.45, rt < 1.05)

ids_upper <- path_info |> filter(outcome == "upper") |> slice_head(n = 4) |> pull(path_id)
ids_lower <- path_info |> filter(outcome == "lower") |> slice_head(n = 2) |> pull(path_id)

paths_bg   <- sim_paths |> filter(path_id %in% c(ids_upper[-1], ids_lower[-1]))
paths_main <- sim_paths |> filter(path_id %in% c(ids_upper[1], ids_lower[1]))
ends_main  <- paths_main |> slice_tail(n = 1, by = path_id)

# --- Simulate many trials for the RT distributions at the boundaries ---------
sim_rt <- rddm(20000, drift = drift, bound = bound, ndt = ndt, zr = zr)
p_upper <- mean(sim_rt$response == 1)

dens_height <- 0.5  # display height of the correct-RT distribution

# Display heights are schematic, not to scale: with these parameters errors
# make up only ~5% of trials, so the error distribution is scaled relative
# to its own peak (at 45% of the correct distribution's height) to stay
# visible — as is conventional in DDM illustrations
dens_scaled <- function(rts, height) {
  d <- density(rts, from = ndt, to = x_lim, adjust = 1.2)
  tibble(time = d$x, h = d$y / max(d$y) * height)
}

dens_upper <- dens_scaled(sim_rt$rt[sim_rt$response == 1], dens_height)
dens_lower <- dens_scaled(sim_rt$rt[sim_rt$response == 0], dens_height * 0.45)

y_top <- bound + dens_height + 0.12
y_bot <- -max(dens_lower$h) - 0.28

# --- Build the figure ---------------------------------------------------------
p_schematic <- ggplot() +

  # RT distributions attached to the boundaries
  geom_ribbon(data = dens_upper,
              aes(x = time, ymin = bound, ymax = bound + h),
              fill = eam_palette[1], alpha = 0.30) +
  geom_line(data = dens_upper, aes(x = time, y = bound + h),
            color = eam_palette[1], linewidth = 0.5) +
  geom_ribbon(data = dens_lower,
              aes(x = time, ymin = -h, ymax = 0),
              fill = eam_palette[4], alpha = 0.30) +
  geom_line(data = dens_lower, aes(x = time, y = -h),
            color = eam_palette[4], linewidth = 0.5) +
  annotate("label", x = 0.82, y = bound + dens_height * 0.8,
           label = "Correct RT distribution", hjust = 0,
           size = 3, color = eam_palette[1],
           fill = "#FFFFFFD9", border.colour = NA,
           label.padding = unit(0.12, "lines")) +
  annotate("label", x = 0.82, y = -max(dens_lower$h) - 0.13,
           label = "Error RT distribution", hjust = 0,
           size = 3, color = eam_palette[4],
           fill = "#FFFFFFD9", border.colour = NA,
           label.padding = unit(0.12, "lines")) +

  # Boundary lines
  geom_hline(yintercept = bound, linewidth = 0.8, color = "grey30") +
  geom_hline(yintercept = 0, linewidth = 0.8, color = "grey30") +

  # Non-decision time shading and flat pre-accumulation segment
  annotate("rect", xmin = 0, xmax = ndt, ymin = 0, ymax = bound,
           fill = "grey92", alpha = 0.8) +
  annotate("segment", x = 0, xend = ndt, y = zr * bound, yend = zr * bound,
           color = "grey30", linewidth = 0.6) +
  annotate("text", x = 0.035, y = 0.06,
           label = "Stimulus onset", angle = 90, hjust = 0,
           size = 2.8, color = "grey45") +

  # Evidence paths (background paths faint, one per boundary highlighted)
  geom_line(data = paths_bg,
            aes(x = time, y = evidence, group = path_id, color = outcome),
            linewidth = 0.35, alpha = 0.35) +
  geom_line(data = paths_main,
            aes(x = time, y = evidence, group = path_id, color = outcome),
            linewidth = 0.7, alpha = 0.95) +
  geom_point(data = ends_main,
             aes(x = time, y = evidence, color = outcome), size = 2.6) +
  scale_color_manual(values = c(upper = eam_palette[1],
                                lower = eam_palette[4]), guide = "none") +

  # Starting point
  geom_point(aes(x = ndt, y = zr * bound), size = 2.6, color = "grey30") +

  # Drift: arrow along the mean accumulation slope
  annotate("segment", x = ndt + 0.03, xend = ndt + 0.27,
           y = zr * bound + 0.06, yend = zr * bound + 0.06 + drift * 0.24,
           arrow = arrow(length = unit(0.11, "inches"), type = "closed"),
           color = "grey20", linewidth = 0.8) +
  annotate("label", x = ndt + 0.33, y = zr * bound + 0.06 + drift * 0.24,
           label = "Drift~(italic(delta))", parse = TRUE,
           hjust = 0, vjust = 1.1, size = 3.4, color = "grey20",
           fill = "#FFFFFFD9", border.colour = NA,
           label.padding = unit(0.12, "lines")) +

  # Boundary separation (alpha), right side
  annotate("segment", x = x_lim - 0.04, xend = x_lim - 0.04, y = 0, yend = bound,
           arrow = arrow(ends = "both", length = unit(0.08, "inches")),
           color = "grey40", linewidth = 0.4) +
  annotate("text", x = x_lim - 0.06, y = bound / 2,
           label = "italic(alpha)", parse = TRUE,
           hjust = 1, size = 4, color = "grey30") +

  # Relative starting point (zeta_r), dotted
  annotate("segment", x = ndt - 0.015, xend = ndt - 0.015,
           y = 0, yend = zr * bound,
           arrow = arrow(ends = "both", length = unit(0.06, "inches")),
           color = "grey40", linewidth = 0.4, linetype = "dotted") +
  annotate("text", x = ndt - 0.05, y = zr * bound / 2,
           label = "italic(zeta)[r]", parse = TRUE,
           hjust = 1, size = 3.4, color = "grey30") +

  # Non-decision time (tau) bracket under the flat segment
  annotate("segment", x = 0, xend = ndt,
           y = zr * bound - 0.14, yend = zr * bound - 0.14,
           arrow = arrow(ends = "both", length = unit(0.06, "inches")),
           color = "grey40", linewidth = 0.4) +
  annotate("text", x = ndt / 2, y = zr * bound - 0.26,
           label = "italic(tau)", parse = TRUE, size = 3.6, color = "grey30") +

  # Boundary labels (inside the panel, right-aligned)
  annotate("text", x = x_lim - 0.1, y = bound - 0.05,
           label = "Upper boundary (correct response)",
           hjust = 1, vjust = 1, size = 3, color = "grey40") +
  annotate("text", x = x_lim - 0.1, y = 0.05,
           label = "Lower boundary (error response)",
           hjust = 1, vjust = 0, size = 3, color = "grey40") +

  # Axes
  labs(x = "Time (s)", y = "Evidence") +
  scale_x_continuous(limits = c(0, x_lim), expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(breaks = NULL) +
  coord_cartesian(ylim = c(y_bot, y_top)) +
  clean_plot() +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank(),
    panel.grid   = element_blank()
  )

ggsave(here("figures", "03_diffusion_schematic.pdf"), p_schematic,
       width = 6.5, height = 4.8, device = "pdf")

cat("Diffusion schematic saved.\n")
