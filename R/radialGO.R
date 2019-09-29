# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'



require(visNetwork, quietly = TRUE)
require(GO.db, quietly = TRUE)

# topGO enrichment analysis returns a list of GO IDs with a score (KS or Fischer)
# need to grab nodes, organize with edges (relationships) and plot with visNetwork
# TODO maybe use topGO's functions which already get the relationships
# see: https://rdrr.io/bioc/topGO/src/R/topGOviz.R


hello <- function() {
  print("Hello, world!")
}
