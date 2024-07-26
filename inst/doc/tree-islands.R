## ----col-trees-by-score-------------------------------------------------------
# Load required libraries
library("TreeTools", quietly = TRUE)
library("TreeDist")

# Generate a set of trees
trees <- as.phylo(as.TreeNumber(BalancedTree(16)) + 0:100 - 15, 16)


