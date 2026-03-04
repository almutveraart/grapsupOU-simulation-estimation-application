# ============================================================
# Identifiability of (alpha, c) for the graph supOU process
# Gamma(alpha, 1) marginal, a* = 1, Delta = 1
#
# Autocorrelation eigenvalue function:
#   rho(h; alpha, c) = (1 + (1 + c) * h)^(1 - alpha)
#
# Population loss (cf. eq. (X) in the paper):
#   L(alpha, c) = (1 / N*) * sum_{h=1}^{N*}
#                 [ rho(h; alpha_0, c_0) - (1 + (1+c)*h)^(1-alpha) ]^2
#
# Output: 4 individual panels of log10 L(alpha, c) as smooth
#         heatmaps with contour overlays (EPS + PNG), sharing a
#         common colour scale, plus a standalone horizontal
#         colour bar (identifiability_colorbar.eps / .png).
#         Panels correspond to Figure 1 in the paper.
#
# Minimum requirements: R >= 4.1.0, ggplot2 >= 3.3.0
# Tested with: ggplot2 3.5.1, dplyr 1.1.4
# ============================================================

library(ggplot2)
library(dplyr)

# ── Output directory (change if needed) ───────────────────────────────────────
out_dir <- "."
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ── Common theme ──────────────────────────────────────────────────────────────
my_theme <- theme(
  axis.title   = element_text(size = 30),
  axis.text    = element_text(size = 30),
  legend.text  = element_text(size = 20),
  legend.title = element_text(size = 25)
)

# ── Model: ACF and population MSE loss ────────────────────────────────────────
rho_gamma <- function(h, alpha, c) {
  base <- 1 + (1 + c) * h
  ifelse(base > 0, base^(1 - alpha), 0)
}

# L(alpha, c): mean squared deviation from true ACF at lags 1, ..., Nstar
loss_fn <- function(alpha, c, alpha0, c0, Nstar = 30) {
  hs <- seq_len(Nstar)
  mean((rho_gamma(hs, alpha, c) - rho_gamma(hs, alpha0, c0))^2)
}

# ── Parameter grid ────────────────────────────────────────────────────────────
N_grid     <- 120
alpha_grid <- seq(1.05, 4.5,  length.out = N_grid)
c_grid     <- seq(-0.99, 0.99, length.out = N_grid)

# ── Scenarios (alpha_0, c_0) ──────────────────────────────────────────────────
scenarios <- list(
  list(alpha0 = 1.5, c0 = -0.8),   # panel 1: long memory
  list(alpha0 = 2.5, c0 = -0.8),   # panel 2: short memory
  list(alpha0 = 2.5, c0 =  0.0),   # panel 3: no network effect
  list(alpha0 = 4.0, c0 = -0.5)    # panel 4: well-identified
)

# ── Loss surface data ─────────────────────────────────────────────────────────
make_loss_df <- function(alpha0, c0) {
  expand.grid(alpha = alpha_grid, c = c_grid) |>
    mutate(
      loss     = mapply(loss_fn, alpha, c,
                        MoreArgs = list(alpha0 = alpha0, c0 = c0)),
      log_loss = log10(loss + 1e-14)   # floor avoids -Inf at exact minimum
    )
}

# ── Compute all panels and derive shared colour limits ────────────────────────
all_dfs <- lapply(scenarios, function(s) make_loss_df(s$alpha0, s$c0))

global_limits <- range(sapply(all_dfs, function(df) range(df$log_loss)))

# ── Build and save each panel (no legend) ────────────────────────────────────
file_names <- c(
  "identifiability_panel1",
  "identifiability_panel2",
  "identifiability_panel3",
  "identifiability_panel4"
)

for (k in seq_along(scenarios)) {
  s  <- scenarios[[k]]
  df <- all_dfs[[k]]
  
  p <- ggplot(df, aes(x = c, y = alpha)) +
    geom_raster(aes(fill = log_loss), interpolate = TRUE) +
    geom_contour(aes(z = log_loss), colour = "white",
                 alpha = 0.25, bins = 12) +
    annotate("point", x = s$c0, y = s$alpha0,
             color = "red", shape = 3, size = 5, stroke = 2) +
    scale_fill_viridis_c(
      option = "viridis",
      limits = global_limits      # shared scale across all panels
    ) +
    labs(
      x = expression(italic(c)),
      y = expression(italic(alpha))
    ) +
    coord_cartesian(xlim = c(-0.99, 0.99), ylim = c(1.05, 4.5)) +
    theme_minimal(base_size = 10) +
    my_theme +
    theme(
      plot.title      = element_blank(),
      plot.margin     = margin(6, 10, 6, 6),
      legend.position = "none"
    )
  
  ggsave(file.path(out_dir, paste0(file_names[k], ".eps")), p,
         width = 8, height = 8, device = cairo_ps)
  ggsave(file.path(out_dir, paste0(file_names[k], ".png")), p,
         width = 8, height = 8, dpi = 150)
  
  message("Panel ", k, " saved: ", file_names[k])
}

# ── Standalone horizontal colour bar ─────────────────────────────────────────
# Build a thin data frame spanning the full shared range and plot as a
# 1-row raster, keeping only the colour bar via theme_void()
colorbar_df <- data.frame(
  x = seq(global_limits[1], global_limits[2], length.out = 500),
  y = 1
)

p_cbar <- ggplot(colorbar_df, aes(x = x, y = y, fill = x)) +
  geom_raster(alpha = 0) +
  scale_fill_viridis_c(
    option = "viridis",
    limits = global_limits,
    name   = expression(log[10]*L(alpha, c)),
    guide  = guide_colorbar(
      barwidth       = unit(18, "cm"),
      barheight      = unit(0.6, "cm"),
      title.position = "top",
      title.hjust    = 0.5,
      ticks.colour   = "white",
      direction      = "horizontal"
    )
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title    = element_text(size = 25),
    legend.text     = element_text(size = 20),
    legend.margin   = margin(t = 4)
  )

ggsave(file.path(out_dir, "identifiability_colorbar.eps"), p_cbar,
       width = 10, height = 1.8, device = cairo_ps)
ggsave(file.path(out_dir, "identifiability_colorbar.pdf"), p_cbar,
       width = 10, height = 1.8, device = cairo_pdf)
ggsave(file.path(out_dir, "identifiability_colorbar.png"), p_cbar,
       width = 10, height = 1.8, dpi = 150)
message("Colour bar saved: identifiability_colorbar")