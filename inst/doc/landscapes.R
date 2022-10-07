## ----col-trees-by-score-------------------------------------------------------
# Load required libraries
library("TreeTools", quietly = TRUE)
library("TreeDist")

# Generate a set of trees
trees <- as.phylo(as.TreeNumber(BalancedTree(16)) + 0:100 - 15, 16)

# Create a 2D mapping
distances <- ClusteringInfoDist(trees)
mapping <- cmdscale(distances, 2)

# Score trees according to their balance
scores <- TotalCopheneticIndex(trees)

# Normalize scores
scoreMax <- TCIContext(trees[[1]])[["maximum"]]
scoreMin <- TCIContext(trees[[1]])[["minimum"]]
scores <- scores - scoreMin
scores <- scores / (scoreMax - scoreMin)

# Generate colour palette
col <- colorRamp(c("orange", "blue"))(scores)
rgbCol <- rgb(col, maxColorValue = 255)

# Plot trees, coloured by their score
plot(mapping,
     asp = 1, # Preserve aspect ratio - do not distort distances
     ann = FALSE, axes = FALSE, # Don't label axes: dimensions are meaningless
     col = rgbCol, # Colour trees by score
     pch = 16
     )

# Add a legend
SpectrumLegend(palette = rgb(colorRamp(c("orange", "blue"))(0:100 / 100) / 255),
               legend = c(scoreMin, scoreMax))

## ----contoured-plot-----------------------------------------------------------
# Use an inverse distance weighting to interpolate between measured points
Predict <- function (x, y) {
  Distance <- function (a, b) {
    apply(a, 2, function (pt) sqrt(colSums((pt - b) ^ 2)))
  }
  predXY <- rbind(x, y)
  dists <- Distance(t(mapping), predXY)
  invDist <- 1 / dists
  weightings <- invDist / rowSums(invDist)

  # Return:
  colSums(scores * t(weightings))
}

# Generate grid for contour plot
resolution <- 32
xLim <- range(mapping[, 1]) * 1.1
yLim <- range(mapping[, 2]) * 1.11
x <- seq(xLim[1], xLim[2], length.out = resolution)
y <- seq(yLim[1], yLim[2], length.out = resolution)
z <- outer(x, y, Predict) # Predicted values for each grid square

# Plot
filled.contour(
  x, y, z,
  asp = 1, # Preserve aspect ratio - do not distort distances
  ann = FALSE, axes = FALSE, # Don't label axes: dimensions are meaningless
  plot.axes = {points(mapping, xpd = NA)} # Use filled.contour coordinates
)

## ----plot-3d, message = FALSE-------------------------------------------------
if (requireNamespace("plotly", quietly = TRUE)) {
  library("plotly", quietly = TRUE)
  fig <- plot_ly(x = x, y = y, z = z)
  fig <- fig %>% add_surface()
  fig
} else {
  print("Run `install.packages('plotly')` to view this output")
}

