## ---- fig.width=6, out.width="90%", fig.align="center"------------------------
library("TreeDist")
tree1 <- ape::read.tree(text = '(A, ((B, ((C, D), (E, F))), (G, (H, (I, J, K)))));')
tree2 <- ape::read.tree(text = '(A, (B, (C, D, E, (J, K)), (F, (G, H, I))));')
VisualizeMatching(NyeSimilarity, tree1, tree2,
                  Plot = TreeDistPlot, matchZeros = FALSE)

## -----------------------------------------------------------------------------
NyeSimilarity(tree1, tree2, normalize = FALSE) / 8
NyeSimilarity(tree1, tree2, normalize = 8)

## -----------------------------------------------------------------------------
NyeSimilarity(tree1, tree2,
                  normalize = min(TreeTools::NSplits(list(tree1, tree2))))

## -----------------------------------------------------------------------------
NyeSimilarity(tree1, tree2, normalize = min)

## -----------------------------------------------------------------------------
NyeSimilarity(list(tree1, tree2), list(tree1, tree2), normalize = pmin)

## -----------------------------------------------------------------------------
NyeSimilarity(tree1, tree2, normalize = TRUE)

## -----------------------------------------------------------------------------
TreeTools::NSplits(tree1)

## -----------------------------------------------------------------------------
NyeSimilarity(tree1, tree2, normalize = TreeTools::NSplits(tree1))

## -----------------------------------------------------------------------------
library("Quartet", exclude = "RobinsonFoulds")
expectedQD <- 2 / 3
normalizedQD <- QuartetDivergence(QuartetStatus(tree1, tree2),
                                  similarity = FALSE) / expectedQD

## ---- fig.width=7, fig.height=4, message=FALSE--------------------------------
if (requireNamespace("TreeDistData", quietly = TRUE)) {
  library("TreeDistData", exclude = "PairwiseDistances")
  data("randomTreeDistances", package = "TreeDistData")
  
  methods <- c("pid", "cid", "nye", "qd")
  methodCol <- c(pid = "#e15659", cid = "#58a14e", nye = "#edc949",
                 qd = "#af7aa1")
  
  oldPar <- par(cex = 0.7, mar = c(5, 5, 0.01, 0.01))
  nLeaves <- as.integer(dimnames(randomTreeDistances)[[3]])
  plot(nLeaves, type = "n", randomTreeDistances["pid", "mean", ],
       ylim = c(0.54, 1),
       xlab = "Number of leaves",
       ylab = "Normalized distance between random tree pairs")
  
  for (method in methods) {
    dat <- randomTreeDistances[method, , ]
    lines(nLeaves, dat["50%", ], pch = 1, col = methodCol[method])
    polygon(c(nLeaves, rev(nLeaves)), c(dat["25%", ], rev(dat["75%", ])),
            border = NA, col = paste0(methodCol[method], "55"))
    
  }
  
  text(202, randomTreeDistances[methods, "50%", "200"] + 0.02, 
       c("Different phylogenetic information", 
         "Clustering information distance",
         expression(paste(plain("Nye "), italic("et al."))),
         "Quartet divergence"
         ), col = methodCol[methods], pos = 2)
  par(oldPar)
}

## ---- eval = FALSE------------------------------------------------------------
#  expectedCID <- randomTreeDistances["cid", "mean", "9"]
#  ClusteringInfoDistance(tree1, tree2, normalize = TRUE) / expectedCID

## ----fig.align="center", fig.height=1.8, fig.width=6, out.width="80%"---------
testTrees <- list(
  trueTree = ape::read.tree(text = '(a, (b, (c, (d, (e, (f, (g, h)))))));'),
  lackRes = ape::read.tree(text = '(a, (b, c, (d, e, (f, g, h))));'),
  smallErr = ape::read.tree(text = '(a, (c, (b, (d, (f, (e, (g, h)))))));'),
  bigErr = ape::read.tree(text = '(a, (c, (((b, d), (f, h)), (e, g))));')
)
VisualizeMatching(MutualClusteringInfo, testTrees$trueTree, testTrees$lackRes)
points(4, 7.5, pch = 2, cex = 3, col = "#E69F00")
VisualizeMatching(MutualClusteringInfo, testTrees$trueTree, testTrees$smallErr)
points(4, 7.5, pch = 3, cex = 3, col = "#56B4E9")
VisualizeMatching(MutualClusteringInfo, testTrees$trueTree, testTrees$bigErr)
points(4, 7.5, pch = 4, cex = 3, col = "#009E73")

## ---- fig.width=4, fig.align="center", fig.asp=1, out.width="50%"-------------
if (requireNamespace("Ternary", quietly = TRUE)) {
  library("Ternary")
  oldPar <- par(mar = rep(0.1, 4))
  TernaryPlot(alab = "Absent information", blab = "Shared information",
              clab = "Misinformation",
              lab.cex = 0.8, lab.offset = 0.18,
              point = "left", clockwise = FALSE,
              grid.col = "#dedede", grid.minor.lines = 0,
              axis.labels = 0:10 / 10, axis.col = "#aaaaaa")
  
  HorizontalGrid()
  correct <- MutualClusteringInfo(testTrees$trueTree, testTrees)
  resolved <- ClusteringEntropy(testTrees)
  unresolved <- resolved["trueTree"] - resolved
  incorrect <- resolved - correct
  TernaryPoints(cbind(unresolved, correct, incorrect), 
                pch = 1:4, cex = 2, col = Ternary::cbPalette8[1:4])
  par(oldPar)
}

## -----------------------------------------------------------------------------
set.seed(0)
trueTree <- TreeTools::RandomTree(20, root = TRUE)

## -----------------------------------------------------------------------------
treeSearchInstalled <- requireNamespace("TreeSearch", quietly = TRUE)
if (treeSearchInstalled) {
  library("TreeSearch", quietly = TRUE, warn.conflict = FALSE) # for TBR, NNI
  oneAway <- structure(lapply(seq_len(200), function(x) {
    tbrTree <- TBR(trueTree)
    ape::consensus(list(tbrTree,
                        NNI(tbrTree),
                        NNI(tbrTree),
                        NNI(tbrTree)))
  }), class = "multiPhylo")
} else {
  message("Install \"TreeSearch\" to run this example")
}

## -----------------------------------------------------------------------------
if (treeSearchInstalled) {
  threeAway <- structure(lapply(seq_len(200), function(x) {
    tbrTree <- TBR(TBR(TBR(trueTree)))
    ape::consensus(list(tbrTree, 
                        NNI(NNI(tbrTree)),
                        NNI(NNI(tbrTree)),
                        NNI(NNI(tbrTree))))
  }), class = "multiPhylo")
}

## -----------------------------------------------------------------------------
if (treeSearchInstalled) {
  correct1 <- MutualClusteringInfo(trueTree, oneAway)
  correct3 <- MutualClusteringInfo(trueTree, threeAway)
}

## -----------------------------------------------------------------------------
if (treeSearchInstalled) {
  infoInTree1 <- ClusteringEntropy(oneAway)
  infoInTree3 <- ClusteringEntropy(threeAway)
}

## -----------------------------------------------------------------------------
if (treeSearchInstalled) {
  unresolved1 <- ClusteringEntropy(trueTree) - infoInTree1
  unresolved3 <- ClusteringEntropy(trueTree) - infoInTree3
}

## -----------------------------------------------------------------------------
if (treeSearchInstalled) {
  incorrect1 <- infoInTree1 - correct1
  incorrect3 <- infoInTree3 - correct3
}

## ---- collapse=TRUE-----------------------------------------------------------
col1 <- hcl(200, alpha = 0.9)
col3 <- hcl(40, alpha = 0.9)
spec1 <- matrix(col2rgb(col1, alpha = TRUE), nrow = 4, ncol = 181)
spec3 <- matrix(col2rgb(col3, alpha = TRUE), nrow = 4, ncol = 181)
spec1[4, ] <- spec3[4, ] <- 0:180
ColToHex <- function(x) rgb(x[1], x[2], x[3], x[4], maxColorValue = 255)
spec1 <- apply(spec1, 2, ColToHex)
spec3 <- apply(spec3, 2, ColToHex)

## ---- fig.width=7, fig.align="center", fig.asp=5/7, out.width="70%"-----------
if (treeSearchInstalled && requireNamespace("Ternary", quietly = TRUE)) {
  layout(matrix(c(1, 2), ncol = 2), widths = c(5, 2))
  oldPar <- par(mar = rep(0, 4))
  TernaryPlot(alab = "Information absent in degraded tree", 
              blab = "\n\nCorrect information in degraded tree", 
              clab = "Misinformation in degraded tree",
              point = "left", clockwise = FALSE, grid.minor.lines = 0,
              axis.labels = 0:10 / 10)
  
  HorizontalGrid()
  
  coords1 <- cbind(unresolved1, correct1, incorrect1)
  coords3 <- cbind(unresolved3, correct3, incorrect3)
  
  ColourTernary(TernaryDensity(coords1, resolution = 20), spectrum = spec1)
  ColourTernary(TernaryDensity(coords3, resolution = 20), spectrum = spec3)
  
  TernaryDensityContour(coords3, col = col3, nlevels = 4)
  TernaryDensityContour(coords1, col = col1, nlevels = 4)
  
  if (requireNamespace("kdensity", quietly = TRUE)) {
    library("kdensity")
    HorizontalKDE <- function(dat, col, add = FALSE) {
      lty <- 1
      lwd <- 2
      kde <- kdensity(dat)
      kdeRange <- kdensity:::get_range(kde)
      if (add) {
        lines(kde(kdeRange), kdeRange, col = col, lty = lty, lwd = lwd)
      } else {
        plot(kde(kdeRange), kdeRange, col = col, lty = lty, lwd = lwd, 
             ylim = c(0, 1), main = "", axes = FALSE, type = "l")
      }
      # abline(h = 0:10 / 10) # Useful for confirming alignment
    }
  
    par(mar = c(1.8, 0, 1.8, 0)) # align plot limits with ternary plot
    HorizontalKDE(correct1 / infoInTree1, col1, add = FALSE)
    HorizontalKDE(correct3 / infoInTree3, col3, add = TRUE)
    mtext("\u2192 Normalized tree quality \u2192", 2)
  }
  par(oldPar)
} else {
  message("Install \"TreeSearch\" and \"Ternary\" to generate this plot")
}

