# install pacman if not already done
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  tidyverse, tidyterra, scales, lattice, latticeExtra, sp, terra, geodata, rasterVis,
  RColorBrewer, mapview, leaflet, leaflet.extras2, leafsync, viridis, plotly, here
)
