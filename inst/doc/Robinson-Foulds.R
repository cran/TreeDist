## ----init, echo=FALSE, warning=FALSE, message=FALSE---------------------------
library('ape')
library('TreeDist')
standardMargin <- c(0.4, 0.4, 0.8, 0.4)
cbPalette8 <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2", "#D55E00", "#CC79A7")
TwoTreePlot <- function () par(mfrow = c(1, 2), mar = standardMargin, cex = 0.9)

Plot <- function (tree, tree2 = NULL, highlight = character(0),
                  prune = if (is.null(tree2)) integer(0) 
                          else list(integer(0), integer(0))
                  ) {
  ABCDEFGHIJK <- LETTERS[1:11]
  if (is.null(tree2)) {
    if(!all(tree$tip.label %in% 1:11)) {
      tree$tip.label <- match(tree$tip.label, ABCDEFGHIJK)
    }
    
    highlightTip <- match(highlight, ABCDEFGHIJK)
    TreeDistPlot(tree, bold = highlightTip, prune = prune,
                 graft = which(tree$edge[, 2] == 
                                 match(highlightTip, tree$tip.label)))
  } else {
    par(mfrow = c(1, 2), mar = rep(0.4, 4))
    Plot(tree, highlight = highlight, prune = prune[[1]])
    Plot(tree2, highlight = highlight, prune = prune[[2]])
  }
}

## ---- echo=FALSE, fig.width = 6, fig.align = 'center'-------------------------
startPar <- TwoTreePlot()
palette <- colorspace::qualitative_hcl(5, c=42, l=88)
TreeDistPlot(balancedTree <- 
               read.tree(text="(((1, 2), (3, 4)), ((5, 6), (7, 8)));"))
edgelabels(1:5, c(1, 2, 5, 9, 12), adj = c(1.5, 0.5), bg = palette[1]) # Splits
edgelabels(1:6, edge = c(1, 2, 5, 8, 9, 12), bg = palette[2],
           adj = c(-0.5, 0.5)) # Internal edges
nodelabels(1:7, bg = palette[3], adj=c(1, 0.5)) # Nodes
nodelabels(1:6, node = 10:15, bg = palette[4], adj=c(-1, 0.5)) # Clades
legend('bottomleft', col=palette[1:2], pch=15, bty='n',
       c('Splits (5)', 'Internal edges (6)'))
legend('bottomright', col=palette[3:4], pch=15, bty='n',
       c('Nodes (7)', 'Clades (6)'))

TreeDistPlot(read.tree(text="((1, 2, 3, 4), (5, 6, 7, 8));"))
edgelabels(1, 1, adj = c(1.5, 0.5), bg = palette[1])
edgelabels(1:2, edge = c(1, 6), adj = c(-0.5, 0.5), bg = palette[2])
nodelabels(1:3, bg = palette[3], adj = c(1, 0.5))
nodelabels(1:2, node=10:11, bg = palette[4], adj = c(-1, 0.5))


## ----two-trees, echo=FALSE, fig.height = 3, fig.width=6, fig.align = 'center'----
TwoTreePlot()
TreeDistPlot(balancedTree, edge.color = cbPalette8[3])
ape::edgelabels(1:5, c(1, 2, 5, 9, 12), bg = palette[1], adj=c(0.5, -0.25))
ape::edgelabels(c('4:4', rep('2:6', 4)), c(1, 2, 5, 9, 12),
                adj=c(0.5, 1.25), frame = 'n')
ape::edgelabels(c(5.53, rep(3.46, 4)), c(1, 2, 5, 9, 12),
                adj=c(0.5, 2.75), frame = 'n')
legend('topleft', 'Balanced tree', bty='n', inset=c(-0.05, -0.03))

TreeDistPlot(caterpillarTree <- 
               ape::read.tree(text="(1, (2, (3, (4, (5, (6, (7, 8)))))));"),
             edge.color = cbPalette8[2])
ape::edgelabels(1:5, bg = palette[1], c(4, 6, 8, 10, 12), adj = c(0.5, -0.25))
ape::edgelabels(c('2|6', '3|5', '4|4', '3|5', '2|6'), edge = c(4, 6, 8, 10, 12),
                frame = 'n', adj = c(0.5, 1.25))
ape::edgelabels(c(3.46, 5.04, 5.53, 5.04, 3.46), edge = c(4, 6, 8, 10, 12),
                frame = 'n', adj = c(0.5, 2.75))
legend('topleft', 'Asymmetric tree', bty = 'n', inset = c(-0.05, -0.03))

## ----ic-of-splits, display='asis', echo=FALSE, warnings=FALSE-----------------
splitSmall <- 2:4
splitLarge <- 8L - splitSmall

rootedTrees <- c(
  '1' = 1,
  '2' = 1 * 1,
  '3' = 3 * 1,
  '4' = 5 * 3 * 1,
  '5' = 7 * 5 * 3 * 1,
  '6' = 9 * 7 * 5 * 3 * 1
)

matchingTrees <- rootedTrees[splitSmall] * rootedTrees[splitLarge]
names(matchingTrees) <- paste0("Split size: ", splitSmall, ':', splitLarge)

matchingP <- matchingTrees / (11 * 9 * 7 * 5 * 3 * 1 * 1)

ic <- -log(matchingP) / log(2)

knitr::kable(cbind(
      'Matching trees' = paste(matchingTrees, "/ 10 395"),
      '_P_(Match in random tree)' = signif(matchingP, 3),
      'Phylogenetic information content' = paste(signif(ic, 3), 'bits')))

## ----all-8-tip-trees, echo=FALSE, cache=TRUE, fig.width=4, fig.height=4, fig.align='center'----
calculate <- FALSE
if (calculate) {
  all8 <- phangorn::allTrees(8, tip.label = 1:8, rooted = FALSE)
  inBalanced <- RobinsonFoulds(all8, balancedTree, similarity = TRUE)
  inCaterpillar <- RobinsonFoulds(all8, caterpillarTree, similarity = TRUE)
} else {
  inBalanced <- rep(0:5, c(7088, 2708, 512, 76, 10, 1))
  inCaterpillar <- rep(0:5, c(8162, 1808, 350, 64, 10, 1))
}

par(cex = 0.8)
sch <- hist(inCaterpillar + 0.7, breaks = 0:18 / 3 - (1/6),
            main = '8-leaf trees with N common splits', cex.main = 1,
            xlim = c(0, 6), axes = FALSE,
            xlab = 'Splits in common', ylab = 'Number of trees')
sbh <- hist(inBalanced + 0.4, breaks = 0:18 / 3 - (1/6), plot = FALSE)
plot(sch, col = paste0(cbPalette8[2], '44'), add = TRUE)
plot(sbh, col = paste0(cbPalette8[3], '44'), add = TRUE)
text(1/6, 100, paste0('Balanced: ', sum(inBalanced == 0)), 
     pos = 4, srt = 90, cex = 0.7)
text(0.5, 100, paste0('Asymmetric: ', sum(inCaterpillar == 0)), 
     pos = 4, srt = 90, cex = 0.7)


legend('topright', pch = 22,
       pt.cex = 2, col = 'black',
       pt.bg = paste0(cbPalette8[2:3], '44'), bty = 'n',
       c('Asymmetric', 'Balanced'))

axis(1, at = 0:5 + 0.5, labels = 0:5)
axis(2)

## ---- fig.height = 3, fig.width=6, fig.align = 'center'-----------------------
tree1 <- ape::read.tree(text='(1, (2, (3, (4, (5, (6, (7, 8)))))));')
tree2 <- ape::read.tree(text='(1, (2, (3, (4, (5, (7, (6, 8)))))));')
tree3 <- ape::read.tree(text='(1, (2, (3, (5, (4, (6, (7, 8)))))));')

VisualizeMatching(InfoRobinsonFoulds, tree1, tree2, 
                  Plot = TreeDistPlot, prune = 12)

## ---- fig.height = 3, fig.width=6, fig.align = 'center'-----------------------
VisualizeMatching(InfoRobinsonFoulds, tree1, tree3, 
                  Plot = TreeDistPlot, prune = 8)

## ---- fig.height = 3, fig.width=6, fig.align = 'center'-----------------------
tree1 <- ape::read.tree(text='(1, (2, (3, (4, (5, (6, (7, 8)))))));')
tree2 <- ape::read.tree(text='(8, (1, (2, (3, (4, (5, (6, 7)))))));')

VisualizeMatching(RobinsonFouldsMatching, tree1, tree2, Plot = TreeDistPlot)

## ---- fig.height = 3, fig.width=6, fig.align = 'center'-----------------------
tree1 <- ape::read.tree(text='((A, B), ((C, (D, E)), (F, (G, (H, I)))));')
tree2 <- ape::read.tree(text='((A, B), ((C, D, (E, I)), (F, (G, H))));')

Plot(tree1, tree2, highlight = 'I', prune = list(8, integer(0)))

## ---- fig.height = 3, fig.width=6, fig.align = 'center'-----------------------
TwoTreePlot()
VisualizeMatching(RobinsonFouldsMatching, tree1, tree2)

## ---- fig.height = 3, fig.width=6, fig.align = 'center'-----------------------
TwoTreePlot()
VisualizeMatching(RobinsonFouldsMatching,
                  drop.tip(tree1, 'I'),
                  drop.tip(tree2, 'I'))

## ---- fig.height = 3, fig.width=6, fig.align = 'center'-----------------------
TwoTreePlot()
VisualizeMatching(SharedPhylogeneticInfo, tree1, tree2)

## -----------------------------------------------------------------------------
summary(TreeTools::as.Splits(tree1, LETTERS[1:9]))
summary(TreeTools::as.Splits(tree2, LETTERS[1:9]))

## -----------------------------------------------------------------------------
attributes(SharedPhylogeneticInfo(tree1, tree2, reportMatching = TRUE))

## ----echo=FALSE---------------------------------------------------------------
par(startPar)

