## ----initialize---------------------------------------------------------------
library("TreeTools", quietly = TRUE)
balAL <- BalancedTree(letters[1:12])
balAJ <- DropTip(balAL, letters[11:12])
balCL <- DropTip(balAL, letters[1:2])
balGL <- DropTip(balAL, letters[1:6])

pecAL <- PectinateTree(letters[1:12])
pecAF <- DropTip(pecAL, letters[7:12])
pecAJ <- DropTip(pecAL, letters[11:12])
pecCL <- DropTip(pecAL, letters[1:2])

treeList <- list(balAL = balAL, balAJ = balAJ, balCL = balCL, balGL = balGL,
                 pecAL = pecAL, pecAJ = pecAJ, pecCL = pecCL, pecAF = pecAF)

# Define a function to plot two trees
Plot2 <- function(t1, t2, ..., main2 = "") {
  oPar <- par(mfrow = c(1, 2), mar = rep(1, 4), cex = 0.9)
  on.exit(par(oPar)) # Restore original parameters
  plot(t1, ...)
  plot(t2, main = main2)
}

## ----drop-some-plot, fig.width = 8--------------------------------------------
Plot2(balAL, balCL,
      main = "balAL", main2 = "balCL",
      font = c(rep(4, 2), rep(3, 10)) # Emphasize missing tips
      )

## ----drop-some----------------------------------------------------------------
library("TreeDist")
commonTips <- intersect(TipLabels(balAL), TipLabels(balCL))

# How much information is in tree balCL?
ClusteringEntropy(balCL)
# All the information in balCL is also present in balAL
MutualClusteringInfo(balCL, balAL)
# balAL also contains information about leaves A and B,
# so contains more information in total
ClusteringEntropy(balAL)

## ----full-list----------------------------------------------------------------
# Information in common (bits)
raw <- MutualClusteringInfo(treeList)
heatmap(raw, symm = TRUE)
# Normalized
normalized <- MutualClusteringInfo(treeList, normalize = TRUE)
heatmap(normalized, symm = TRUE)

