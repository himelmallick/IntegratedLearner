# =============================================================================
# Reusable hex-sticker generator for R packages (built on {hexSticker})
#
# Usage for ANY package:
#   source("make_hexsticker.R")
#   make_pkg_sticker("mypkg", subtitle = "WHAT IT DOES", filename = "logo.png")
#
# Then, inside the package project:
#   usethis::use_logo("logo.png")      # copies to man/figures/logo.png
#
# Requirements: install.packages(c("hexSticker", "ggplot2", "showtext"))
# =============================================================================

library(hexSticker)
library(ggplot2)
library(showtext)

# -----------------------------------------------------------------------------
# Font helper: use a Google font if online, otherwise fall back silently
# -----------------------------------------------------------------------------
load_sticker_font <- function(family = "Inter", weight = 800) {
  ok <- tryCatch({
    sysfonts::font_add_google(family, family = "sticker_font",
                              regular.wt = weight)
    TRUE
  }, error = function(e) FALSE)
  showtext::showtext_auto()
  if (ok) "sticker_font" else "sans"
}

# -----------------------------------------------------------------------------
# Default motif: several coloured streams (methods) converging into one node,
# followed by a single combined peak. Good for ensemble / integration / meta-
# analysis packages, but swap in any ggplot or an image path via `subplot`.
# -----------------------------------------------------------------------------
motif_converge <- function(n = 5,
                           colours = c("#2EC4B6", "#4EA8DE", "#9B8CFF",
                                       "#FF7A7A", "#F4B942"),
                           line_col = "white",
                           node_col = "#F4B942",
                           peak = TRUE) {
  colours <- rep_len(colours, n)
  node <- c(x = 0.62, y = 0)
  y0 <- seq(1, -1, length.out = n)
  t <- seq(0, 1, length.out = 100)

  # cubic Bezier from (0, y0) to node, flat at the start, converging at the end
  bez <- function(p0, p1, p2, p3) {
    (1 - t)^3 * p0 + 3 * (1 - t)^2 * t * p1 + 3 * (1 - t) * t^2 * p2 + t^3 * p3
  }
  streams <- do.call(rbind, lapply(seq_len(n), function(i) {
    data.frame(
      id = i,
      x  = bez(0, 0.30, node["x"] - 0.28, node["x"]),
      y  = bez(y0[i], y0[i], 0.15 * y0[i], node["y"])
    )
  }))
  starts <- data.frame(x = 0, y = y0, id = seq_len(n))

  p <- ggplot() +
    geom_path(data = streams, aes(x, y, group = id, colour = factor(id)),
              linewidth = 1.6, lineend = "round") +
    geom_point(data = starts, aes(x, y, colour = factor(id)), size = 2.6) +
    scale_colour_manual(values = colours, guide = "none")

  if (peak) {
    # Cauchy density: the heavy-tailed distribution behind the CCT combination
    xs <- seq(node["x"], 1.55, length.out = 200)
    mu <- 1.05; gamma <- 0.11; h <- 1.05
    cd <- data.frame(x = xs, y = h * gamma^2 / ((xs - mu)^2 + gamma^2))
    p <- p +
      geom_ribbon(data = cd, aes(x, ymin = 0, ymax = y),
                  fill = line_col, alpha = 0.18) +
      geom_line(data = cd, aes(x, y), colour = line_col,
                linewidth = 1.1, lineend = "round")
  }

  p +
    geom_point(aes(node["x"], node["y"]), shape = 21, size = 6.2,
               stroke = 0.9, colour = node_col, fill = NA) +
    geom_point(aes(node["x"], node["y"]), size = 3.8, colour = line_col) +
    coord_cartesian(xlim = c(-0.08, 1.6), ylim = c(-1.12, 1.15), expand = FALSE) +
    theme_void() +
    hexSticker::theme_transparent()
}

# -----------------------------------------------------------------------------
# Main wrapper: one call per package. Every argument has a sensible default,
# so typically you only change `package`, `subtitle` and maybe colours.
# -----------------------------------------------------------------------------
make_pkg_sticker <- function(package,
                             subtitle   = NULL,
                             subplot    = motif_converge(),
                             filename   = "logo.png",
                             fill       = "#0E2A47",   # hexagon fill
                             border     = "#F4B942",   # hexagon border
                             text_col   = "white",
                             sub_col    = border,
                             font       = "Inter",
                             p_size     = NULL,        # auto from name length
                             p_y        = 1.47,
                             line_gap   = 1,           # spacing for two-line names
                             s_x = 1.0, s_y = 0.86, s_width = 1.22, s_height = 0.72,
                             dpi        = 300) {
  family <- load_sticker_font(font)

  # long names need a smaller font so they stay inside the hexagon
  # a name with "\n" is drawn as stacked lines (e.g. "Integrated\nLearner")
  name_lines <- strsplit(package, "\n", fixed = TRUE)[[1]]
  longest <- max(nchar(name_lines))
  if (is.null(p_size)) p_size <- max(8, 19 - 1.1 * max(0, longest - 9))

  s <- hexSticker::sticker(
    subplot   = subplot,
    s_x = s_x, s_y = s_y, s_width = s_width, s_height = s_height,
    package   = if (length(name_lines) > 1) "" else package,
    p_family  = family, p_size = p_size, p_color = text_col, p_y = p_y,
    p_fontface = "bold",
    h_fill    = fill, h_color = border, h_size = 1.6,
    url       = "",
    white_around_sticker = FALSE,
    filename  = tempfile(fileext = ".png"), dpi = dpi
  )

  if (length(name_lines) > 1) {
    step <- 0.0105 * p_size * line_gap          # line spacing in hexagon units
    ys <- p_y + step * (rev(seq_along(name_lines)) - 1) - step * (length(name_lines) - 1) / 2
    s <- s + ggplot2::annotate("text", x = 1, y = ys, label = name_lines,
                               family = family, fontface = "bold",
                               colour = text_col, size = p_size)
  }

  # centred subtitle near the bottom edge (sticker()'s url is right-aligned)
  if (!is.null(subtitle)) {
    s <- s + ggplot2::annotate("text", x = 1, y = 0.36, label = subtitle,
                               family = family, fontface = "bold",
                               colour = sub_col, size = 5)
  }
  hexSticker::save_sticker(filename, s, dpi = dpi)
  message("Saved ", normalizePath(filename))
  invisible(s)
}

# =============================================================================
# Example: DAssemble (regenerates all five logo options)
#   Rscript hexsticker/make_hexsticker.R      (from the package root)
# =============================================================================
if (sys.nframe() == 0L) {
  here <- if (dir.exists("hexsticker")) "hexsticker" else "."
  source(file.path(here, "motifs_da.R"))
  out <- function(f) file.path(here, f)

  make_pkg_sticker("DAssemble", subplot = motif_bagging("volcano"),
                   s_width = 1.38, s_height = 0.86, s_y = 0.86,
                   filename = out("logo_bagging_volcano.png"))
  make_pkg_sticker("DAssemble", subplot = motif_bagging("bulb"),
                   s_width = 1.38, s_height = 0.86, s_y = 0.86,
                   filename = out("logo_bagging_bulb.png"))
  make_pkg_sticker("DAssemble", subplot = motif_funnel(),
                   s_width = 1.22, s_height = 0.92, s_y = 0.82,
                   filename = out("logo_funnel.png"))
  make_pkg_sticker("DAssemble", subplot = motif_hub(),
                   s_width = 1.40, s_height = 0.92, s_y = 0.84,
                   filename = out("logo_hub.png"))
  make_pkg_sticker("DAssemble", subplot = motif_deck(),
                   s_width = 1.40, s_height = 0.66, s_y = 0.88,
                   filename = out("logo_deck.png"))
}
