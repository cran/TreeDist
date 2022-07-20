## ----generate-trees-----------------------------------------------------------
# Load trees
library("TreeTools", quietly = TRUE)
batch1 <- as.phylo(1:60, 8) # Generate 60 similar trees 
batch2 <- as.phylo(seq(200, 800, length.out = 30), 8) # A separate batch of 30 trees
styles <- c(1, 2) # Select plotting colours / symbols
treeStyle <- rep(styles, c(length(batch1), length(batch2)))

# Calculate distances
library("TreeDist")
distances <- ClusteringInfoDistance(c(batch1, batch2))

# Construct over-simple 2D PCoA mapping
mapping <- cmdscale(distances, k = 2)

## ----plot-mapping, fig.align = "center"---------------------------------------
# Plot mapping
par(mar = rep(0, 4))
plot(mapping,
     asp = 1, # Preserve aspect ratio - do not distort distances
     ann = FALSE, axes = FALSE, # Don't label axes: dimensions are meaningless
     col = treeStyle, # Colour
     pch = treeStyle # Plotting symbol
     )
legend("left", c("Batch 1", "Batch 2"), col = styles, pch = styles)

## ----compare-dist-from-median, fig.align = "center"---------------------------
# Calculate median trees
median1 <- median(batch1)
median2 <- median(batch2)

# Compute distance from each tree to the median of its batch
dist1 <- ClusteringInfoDist(batch1, median1)
dist2 <- ClusteringInfoDist(batch2, median2)

# Set resolution of histogram
nBreaks <- 10
breaks <- seq(0, max(dist1, dist2), length.out = nBreaks)

# Plot first distance set
hist(dist1, col = "#00000022", breaks = breaks,
     main = "Distance from median of batch",
     xlab = "Clustering information distance",
     ylim = c(0, 25) # Omit this line to infer Y axis limit from first batch.
     )

# Add second distance set
hist(dist2, col = "#ff000022", breaks = breaks, add = TRUE)

# Add legend
legend("topleft", c("Batch 1", "Batch 2"),
       fill = c("#00000022",  "#ff000022"))

## ----rogue-detection----------------------------------------------------------
# Create tree set with a rogue taxon
batch3 <- AddTipEverywhere(as.phylo(7, 7), "t8")

# Set up plotting area
par(mfrow = c(2, 2), mar = rep(0.4, 4))

# Plot naive strict consensus
plot(consensus(batch1, p = 1))
plot(consensus(batch3, p = 1))

if (requireNamespace("Rogue", quietly = TRUE)) {
  cons1 <- ConsensusWithout(batch1, p = 1,
                            Rogue::QuickRogue(batch1, p = 1)[-1, "taxon"])
  cons3 <- ConsensusWithout(batch3, p = 1,
                            Rogue::QuickRogue(batch3, p = 1)[-1, "taxon"])
  
  # The information content of each tree gives a measure of its resolution,
  # accounting for omitted rogue leaves
  SplitwiseInfo(cons1) # 8.5 bits
  SplitwiseInfo(cons3) # 15.1 bits: higher resolution indicates that these
  # trees are more similar, notwithstanding rogue taxa.
  
  # Plot the trees
  plot(cons1)
  plot(cons3)
} else {
  message("The package 'Rogue' is required to run this example.")
}

