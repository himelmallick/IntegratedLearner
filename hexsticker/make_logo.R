# =============================================================================
# IntegratedLearner hex logo
#   Run from the package root (RStudio project):  source("hexsticker/make_logo.R")
# Requirements: install.packages(c("hexSticker", "showtext", "ggplot2", "usethis", "pkgdown"))
# =============================================================================
source("hexsticker/make_hexsticker.R")
source("hexsticker/motifs_il.R")

make_pkg_sticker(
  package  = "Integrated\nLearner",          # two lines, so the name stays legible
  subplot  = motif_stack("converge"),        # omics layers -> one integrated prediction
  p_size   = 17, p_y = 1.43, line_gap = 1.12,
  s_width  = 1.22, s_height = 0.80, s_y = 0.76,
  filename = "hexsticker/logo.png"
)

# install it as the package / website logo and refresh the browser-tab icons
usethis::use_logo("hexsticker/logo.png", retina = TRUE)
pkgdown::build_favicons(overwrite = TRUE)
