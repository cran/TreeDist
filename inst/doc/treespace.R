## ----secret-header, echo=FALSE------------------------------------------------
set.seed(0)
protoclust <- if (requireNamespace("protoclust", quietly = TRUE)) {
  protoclust::protoclust
} else {
  hclust
}

## ----generate-trees-----------------------------------------------------------
library("TreeTools", quietly = TRUE)
treeNumbers <- c(1:220)
trees <- as.phylo(treeNumbers, 8)
spectrum <- hcl.colors(220, "plasma")
treeCols <- spectrum[treeNumbers]

## ----calculate-distances, message=FALSE---------------------------------------
library("TreeDist")
distances <- ClusteringInfoDistance(trees)

## ----map----------------------------------------------------------------------
mapping <- cmdscale(distances, k = 12)

## ----plot-mapping-2d, fig.asp = 1, fig.width = 3, fig.align="center"----------
par(mar = rep(0, 4))
plot(mapping,
     asp = 1, # Preserve aspect ratio - do not distort distances
     ann = FALSE, axes = FALSE, # Don't label axes: dimensions are meaningless
     col = treeCols, pch = 16
     )

## ----clustering, fig.align="center"-------------------------------------------
possibleClusters <- 2:10

pamClusters <- lapply(possibleClusters, function(k) cluster::pam(distances, k = k))
pamSils <- vapply(pamClusters, function(pamCluster) {
  mean(cluster::silhouette(pamCluster)[, 3])
}, double(1))

bestPam <- which.max(pamSils)
pamSil <- pamSils[bestPam]
pamCluster <- pamClusters[[bestPam]]$cluster

hTree <- protoclust(distances)
hClusters <- lapply(possibleClusters, function(k) cutree(hTree, k = k))
hSils <- vapply(hClusters, function(hCluster) {
  mean(cluster::silhouette(hCluster, distances)[, 3])
}, double(1))


bestH <- which.max(hSils)
hSil <- hSils[bestH]
hCluster <- hClusters[[bestH]]

plot(pamSils ~ possibleClusters,
     xlab = "Number of clusters", ylab = "Silhouette coefficient",
     ylim = range(c(pamSils, hSils)))
points(hSils ~ possibleClusters, pch = 2)
legend("topright", c("PAM", "Hierarchical"), pch = 1:2)

## ----chosen-cluster-----------------------------------------------------------
cluster <- hClusters[[2 - 1]]

## ----h-tree, fig.align="center"-----------------------------------------------
class(hTree) <- "hclust"
par(mar = c(0, 0, 0, 0))
plot(hTree, labels = FALSE, main = "")
points(seq_along(trees), rep(1, length(trees)), pch = 16,
       col = spectrum[hTree$order])

## ----consensus, fig.align="center"--------------------------------------------
par(mfrow = c(1, 2), mar = rep(0.2, 4))
col1 <- spectrum[mean(treeNumbers[cluster == 1])]
col2 <- spectrum[mean(treeNumbers[cluster == 2])]
plot(consensus(trees[cluster == 1], p = 0.5),
     edge.color = col1, edge.width = 2, tip.color = col1)
plot(consensus(trees[cluster == 2], p = 0.5),
     edge.color = col2, edge.width = 2, tip.color = col2)

## ----how-many-dims, fig.align="center"----------------------------------------
txc <- vapply(seq_len(ncol(mapping)), function(k) {
  newDist <- dist(mapping[, seq_len(k)])
  MappingQuality(distances, newDist, 10)["TxC"]
}, 0)
plot(txc, xlab = "Dimension")
abline(h = 0.9, lty = 2)

## ----calculate-MST------------------------------------------------------------
mstEnds <- MSTEdges(distances)

## ----plot-mapping-5d, fig.asp = 1, fig.align="center"-------------------------
nDim <- which.max(txc > 0.9)
plotSeq <- matrix(0, nDim, nDim)
plotSeq[upper.tri(plotSeq)] <- seq_len(nDim * (nDim - 1) / 2)
plotSeq <- t(plotSeq[-nDim, -1])
plotSeq[nDim * 1:3] <- (nDim * (nDim - 1) / 2) + 1:3
layout(plotSeq)
par(mar = rep(0.1, 4))

for (i in 2:nDim) for (j in seq_len(i - 1)) {
  # Set up blank plot
  plot(mapping[, j], mapping[, i], ann = FALSE, axes = FALSE, frame.plot = TRUE,
       type = "n", asp = 1, xlim = range(mapping), ylim = range(mapping))
  
  # Plot MST
  MSTSegments(mapping[, c(j, i)], mstEnds,
              col = StrainCol(distances, mapping[, c(j, i)]))
  
  # Add points
  points(mapping[, j], mapping[, i], pch = 16, col = treeCols)

  # Mark clusters
  for (clI in unique(cluster)) {
    inCluster <- cluster == clI
    clusterX <- mapping[inCluster, j]
    clusterY <- mapping[inCluster, i]
    hull <- chull(clusterX, clusterY)
    polygon(clusterX[hull], clusterY[hull], lty = 1, lwd = 2,
            border = "#54de25bb")
    text(mean(clusterX), mean(clusterY), clI, col = "#54de25bb", font = 2)
  }
}
# Annotate dimensions
plot(0, 0, type = "n", ann = FALSE, axes = FALSE)
text(0, 0, "Dimension 2")
plot(0, 0, type = "n", ann = FALSE, axes = FALSE)
text(0, 0, "Dimension 3")
plot(0, 0, type = "n", ann = FALSE, axes = FALSE)
text(0, 0, "Dimension 4")

## ----pid, fig.asp = 1, fig.width = 4, fig.align = "center", echo = FALSE, message = FALSE----
library("TreeDist")
pid_distances <- PhylogeneticInfoDistance(trees)
pid_mapping <- cmdscale(pid_distances, k = 6)
pid_cluster <- cutree(protoclust(pid_distances), k = 2)

par(mar = rep(0, 4))
plot(pid_mapping, ann = FALSE, axes = FALSE, asp = 1,
     col = treeCols, pch = 16)
MSTEdges(pid_distances, TRUE, pid_mapping[, 1], pid_mapping[, 2],
         col = "#bbbbbb", lty = 1)
for (clI in 1:2) {
  inCluster <- pid_cluster == clI
  clusterX <- pid_mapping[inCluster, 1]
  clusterY <- pid_mapping[inCluster, 2]
  hull <- chull(clusterX, clusterY)
  polygon(clusterX[hull], clusterY[hull], lty = 1, lwd = 2,
          border = "#54de25bb")
}

## ----hypervolume, message = FALSE---------------------------------------------
hypervolumeInstalled <- requireNamespace("hypervolume", quietly = TRUE)
if (hypervolumeInstalled) {
  library("hypervolume")
  hv1 <- hypervolume_gaussian(pid_mapping[pid_cluster == 1, 1:3],
                              verbose = FALSE)
  hv2 <- hypervolume_gaussian(pid_mapping[pid_cluster == 2, 1:3],
                              verbose = FALSE)
  hv_dist <- hypervolume_distance(hv1, hv2)
  capture.output(
    hyperset <- hypervolume_set(hv1, hv2, verbose = FALSE,
                              check.memory = FALSE)
  ) -> XX_VerboseNotRespected
  hv_overlap <- hypervolume_overlap_statistics(hyperset)
  hv_dist
  hv_overlap
} else {
  print("Install the 'hypervolume' package to run this example")
}


## ----umatrix, fig.asp = 1, fig.align = "center"-------------------------------
umatrixInstalled <- requireNamespace("Umatrix", quietly = TRUE)
if (umatrixInstalled) {
  map <- Umatrix::esomTrain(as.matrix(distances), Key = seq_along(trees),
                            Epochs = 5, # Increase for better results
                            Lines = 42,
                            Columns = 42,
                            Toroid = FALSE)
  
  Umatrix::plotMatrix(Matrix = map$Umatrix,
                      Toroid = FALSE, FixedRatio = TRUE,
                      TransparentContours = FALSE, Clean = TRUE) +
  ggplot2::geom_point(data = data.frame(x = map$BestMatches[, 3],
                                        y = map$BestMatches[, 2]),
                      shape = 19, color = treeCols, size = 2)

} else {
  message("Install the 'Umatrix' package to run this example")
}

