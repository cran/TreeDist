## ----init, message=FALSE, warning=FALSE, echo = FALSE-------------------------
library('ape')
library('TreeDist')
tree1 <- read.tree(text='((A, B), ((C, (D, E)), (F, (G, (H, I)))));')
tree2 <- read.tree(text='((A, B), ((C, D, (E, I)), (F, (G, H))));')
AtoJ <-  read.tree(text='(((((A, B), C), D), E), (F, (G, (H, (I, J)))));')
swapAJ <-  read.tree(text='(((((J, B), C), D), E), (F, (G, (H, (I, A)))));')

## ---- fig.width=6, fig.align='center'-----------------------------------------
VisualizeMatching(SharedPhylogeneticInfo, tree1, tree2, 
                  Plot = TreeDistPlot, matchZeros = FALSE)
SharedPhylogeneticInfo(tree1, tree2)

## ---- fig.width=6, fig.align='center'-----------------------------------------
VisualizeMatching(SharedPhylogeneticInfo, AtoJ, swapAJ,
                  Plot = TreeDistPlot, matchZeros = FALSE, prune = c(5, 18))

## ---- fig.width=6, fig.align='center'-----------------------------------------
VisualizeMatching(MutualClusteringInfo, AtoJ, swapAJ,
                  Plot = TreeDistPlot, matchZeros = FALSE, prune = c(5, 18))
MutualClusteringInfo(AtoJ, swapAJ)

## ---- fig.width=6, fig.align='center'-----------------------------------------
VisualizeMatching(MutualClusteringInfo, tree1, tree2, 
                  Plot = TreeDistPlot, matchZeros = FALSE)

## ---- fig.width=6, fig.align='center'-----------------------------------------
VisualizeMatching(NyeSimilarity, tree1, tree2, 
                  Plot = TreeDistPlot, matchZeros = FALSE)
NyeSimilarity(tree1, tree2, normalize = FALSE)

## ---- fig.width=6, fig.align='center'-----------------------------------------
JaccardRobinsonFoulds(tree1, tree2, k = 1)
VisualizeMatching(JaccardRobinsonFoulds, tree1, tree2,
                  Plot = TreeDistPlot, matchZeros = FALSE)

JRF2 <- function(...) JaccardRobinsonFoulds(k = 2, ...)
JRF2(tree1, tree2)
VisualizeMatching(JRF2, tree1, tree2,
                  Plot = TreeDistPlot, matchZeros = FALSE)

## ---- fig.width=6, fig.height = 4,fig.align = 'center', echo = FALSE----------
oldPar <- par(mar = c(4, 4, 0, 0))
x <- seq(1, 20, by = 0.1)
cbPalette8 <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2", "#D55E00", "#CC79A7") # dput(Ternary::cbPalette8)
plot(0, xlim = range(x), ylim = c(3.5, 9.5),
     xlab = expression(italic('k')),
     ylab = 'Jaccard\u2013Robinson\u2013Foulds distance',
     pch = 3, type = 'n')

abline(h = RobinsonFoulds(tree1, tree2),
       lty = 'dashed', col = cbPalette8[6])
abline(h = NyeSimilarity(tree1, tree2, similarity = FALSE),
       lty = 'dotted', col = cbPalette8[8])

lines(x = x, y = vapply(x, JaccardRobinsonFoulds, tree1 = tree1, tree2 = tree2,
                        allowConflict = TRUE, 0), col = cbPalette8[5])
lines(x = x, y = vapply(x, JaccardRobinsonFoulds, tree1 = tree1, tree2 = tree2,
                        0), col = cbPalette8[4])

legend('right', c('Robinson\u2013Foulds', 'JRF, no conflict', 'JRF, conflict ok', 'Nye et al.'),
       lty = c('dashed', 'solid', 'solid', 'dotted'), 
       col = cbPalette8[c(6, 4, 5, 8)], bty = 'n')
par(oldPar)

## ---- fig.width=6, fig.align='center'-----------------------------------------
MatchingSplitDistance(tree1, tree2)
VisualizeMatching(MatchingSplitDistance, tree1, tree2,
                  Plot = TreeDistPlot, matchZeros = FALSE)

## ---- fig.width=6, fig.align='center'-----------------------------------------
MatchingSplitDistance(read.tree(text='((a, b, c, d, e, f), (g, h, i, j));'),
                      read.tree(text='((a, b, c, d, e, i, j), (g, h, f));'))

## -----------------------------------------------------------------------------
TreeTools::SplitInformation(5, 2)

## ---- fig.width=6, fig.align='center'-----------------------------------------
MatchingSplitInfoDistance(tree1, tree2)
VisualizeMatching(MatchingSplitInfoDistance, tree1, tree2,
                  Plot = TreeDistPlot, matchZeros = FALSE)

