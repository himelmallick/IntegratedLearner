# =============================================================================
# Sticker motifs for multi-omics integration packages (IntegratedLearner).
# Visual language follows the IntegratedLearner flowchart: each omics view is a
# coloured data layer; the integrated prediction is gold.
# =============================================================================
library(ggplot2)

VIEW_COLS <- c("#6FA8FF", "#FF9F43", "#2EC4B6", "#B07CFF")   # blue, orange, teal, purple
GOLD <- "#F4B942"

# a matrix of tiles with its lower-left corner at (x0, y0)
.grid <- function(x0, y0, nr = 3, nc = 3, size = 0.42, gap = 0.08, col, alpha = 1) {
  g <- expand.grid(i = seq_len(nc) - 1, j = seq_len(nr) - 1)
  data.frame(xmin = x0 + g$i * (size + gap), ymin = y0 + g$j * (size + gap),
             xmax = x0 + g$i * (size + gap) + size, ymax = y0 + g$j * (size + gap) + size,
             col = col, alpha = alpha)
}
.tiles <- function(d) geom_rect(data = d, aes(xmin = xmin, xmax = xmax, ymin = ymin,
                                              ymax = ymax), fill = d$col, alpha = d$alpha)
.wire <- function(x, y, xend, yend, curv = 0, head = FALSE) {
  force(x); force(y); force(xend); force(yend)
  d <- data.frame(x = x, y = y, xend = xend, yend = yend)
  a <- if (head) grid::arrow(length = grid::unit(0.03, "npc"), type = "closed")
  if (curv == 0) geom_segment(data = d, aes(x = x, y = y, xend = xend, yend = yend),
                              colour = "white", linewidth = 0.7, arrow = a)
  else geom_curve(data = d, aes(x = x, y = y, xend = xend, yend = yend), colour = "white",
                  linewidth = 0.7, curvature = curv, arrow = a)
}
.done <- function(p, xlim, ylim) p + coord_fixed(xlim = xlim, ylim = ylim, expand = FALSE) +
  theme_void() + hexSticker::theme_transparent()

# -----------------------------------------------------------------------------
# Stacked omics layers -> integrated prediction (IntegratedLearner uses "converge")
#   style = "gold_top" : a gold "integrated" sheet sits on top of the omics sheets
#           "beam"     : a gold beam runs through every layer to one gold point
#           "converge" : each layer sends a line to one gold node on the right
#           "pyramid"  : lines rise from every layer to a gold node above the stack
# -----------------------------------------------------------------------------
.sheet_df <- function(z, col, n = 4, size = 0.5, gap = 0.1, skew = 0.55, flat = 0.42, x0 = 0) {
  g <- .grid(0, 0, nr = n, nc = n, size = size, gap = gap, col = col)
  do.call(rbind, lapply(seq_len(nrow(g)), function(i) {
    r <- g[i, ]
    xs <- c(r$xmin, r$xmax, r$xmax, r$xmin); ys <- c(r$ymin, r$ymin, r$ymax, r$ymax)
    data.frame(x = x0 + xs + skew * ys, y = flat * ys + z, id = paste(z, i))
  }))
}
# centre of a sheet at height z (same geometry as .sheet_df)
.sheet_centre <- function(z, n = 4, size = 0.5, gap = 0.1, skew = 0.55, flat = 0.42, x0 = 0) {
  w <- n * size + (n - 1) * gap
  c(x = unname(x0 + w / 2 + skew * w / 2), y = unname(z + flat * w / 2))
}

motif_stack <- function(style = c("gold_top", "beam", "converge", "pyramid"),
                        cols = VIEW_COLS, gap_z = 0.95) {
  style <- match.arg(style)
  zs <- (0:3) * gap_z
  p <- ggplot()
  draw_sheets <- function(p, zs, cols, alpha = 0.92) {
    for (k in seq_along(zs)) {
      s <- .sheet_df(zs[k], cols[k])
      p <- p + geom_polygon(data = s, aes(x, y, group = id), fill = cols[k], alpha = alpha)
    }
    p
  }

  if (style == "gold_top") {
    p <- draw_sheets(p, zs, cols)
    top <- .sheet_df(max(zs) + gap_z * 1.25, GOLD)
    p <- p + geom_polygon(data = top, aes(x, y, group = id), fill = GOLD)
    return(.done(p, c(-0.15, 4.2), c(-0.15, max(zs) + gap_z * 1.25 + 1.2)))
  }

  if (style == "beam") {
    ctr <- .sheet_centre(0)
    top_y <- unname(.sheet_centre(max(zs))["y"]) + 0.95
    # beam behind the sheets, then sheets, then the visible beam segments between them
    p <- p + geom_segment(aes(x = ctr["x"], xend = ctr["x"], y = ctr["y"] - 0.2, yend = top_y),
                          colour = GOLD, linewidth = 1.4)
    p <- draw_sheets(p, zs, cols)
    segs <- data.frame(y = sapply(zs, function(z) .sheet_centre(z)["y"]))
    p <- p +
      geom_segment(data = data.frame(y0 = segs$y[-nrow(segs)] + 0.05, y1 = segs$y[-1] - 0.05),
                   aes(x = ctr["x"], xend = ctr["x"], y = y0, yend = y1),
                   colour = GOLD, linewidth = 1.4) +
      geom_segment(aes(x = ctr["x"], xend = ctr["x"], y = max(segs$y), yend = top_y),
                   colour = GOLD, linewidth = 1.4) +
      geom_point(aes(ctr["x"], top_y), colour = GOLD, size = 8, alpha = 0.25) +
      geom_point(aes(ctr["x"], top_y), colour = GOLD, size = 4.6)
    return(.done(p, c(-0.15, 4.2), c(-0.2, top_y + 0.65)))
  }

  if (style == "converge") {
    mid <- unname(.sheet_centre(mean(zs))["y"])
    node <- c(x = 5.4, y = mid)
    starts <- data.frame(x = sapply(zs, function(z) .sheet_centre(z)["x"]) + 1.35,
                         y = sapply(zs, function(z) .sheet_centre(z)["y"]), k = seq_along(zs))
    p <- draw_sheets(p, zs, cols)
    for (k in seq_along(zs)) {
      p <- p + .wire(starts$x[k], starts$y[k], node[["x"]] - 0.2, node[["y"]],
                     curv = c(0.25, 0.12, -0.12, -0.25)[k])
    }
    p <- p + geom_point(data = starts, aes(x, y), colour = cols[starts$k], size = 2.2) +
      geom_point(data = data.frame(x = node[["x"]], y = node[["y"]]), aes(x, y),
                 colour = GOLD, size = 9, alpha = 0.25) +
      geom_point(data = data.frame(x = node[["x"]], y = node[["y"]]), aes(x, y),
                 colour = GOLD, size = 5.4)
    return(.done(p, c(-0.15, 6.45), c(-0.15, max(zs) + 1.25)))
  }

  # pyramid
  node <- c(x = unname(.sheet_centre(0)["x"]), y = max(zs) + 2.1)
  for (k in seq_along(zs)) {
    sc <- .sheet_centre(zs[k])
    for (sx in c(-1.2, 1.2)) p <- p + .wire(sc[["x"]] + sx, sc[["y"]], node[["x"]], node[["y"]])
  }
  nd <- data.frame(x = node[["x"]], y = node[["y"]])
  p <- draw_sheets(p, zs, cols) +
    geom_point(data = nd, aes(x, y), colour = GOLD, size = 9, alpha = 0.25) +
    geom_point(data = nd, aes(x, y), colour = GOLD, size = 5.4)
  .done(p, c(-0.15, 4.2), c(-0.15, node[["y"]] + 0.35))
}
