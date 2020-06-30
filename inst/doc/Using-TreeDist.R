## ----set-up-------------------------------------------------------------------
tree1 <- ape::read.tree(text = '(A, ((B, (C, (D, E))), ((F, G), (H, I))));')
tree2 <- ape::read.tree(text = '(A, ((B, (C, (D, (H, I)))), ((F, G), E)));')

## ----load-package, message=FALSE----------------------------------------------
library('TreeDist')

## ----measure-distance---------------------------------------------------------
distance <- TreeDistance(tree1, tree2)

## ----multi-trees--------------------------------------------------------------
oneTree <- ape::rtree(11)
twoTrees <- structure(list(one = ape::rtree(11), two = ape::rtree(11)),
                      class = 'multiPhylo')
threeTrees <- list(a = ape::rtree(11), b = ape::rtree(11), c = ape::rtree(11))

TreeDistance(oneTree, twoTrees)
TreeDistance(twoTrees, threeTrees)

## ----visualise-matching, fig.align='center', fig.width=8, out.width='90%'-----
VisualizeMatching(ClusteringInfoDistance, tree1, tree2)

## ----write-vis-matching, fig.align='center', fig.width=8, out.width='90%'-----
ClusteringInfoDistance(tree1, tree2, reportMatching = TRUE)

## -----------------------------------------------------------------------------
splits <- as.character(TreeTools::as.Splits(tree2))
splits

## ----named-splits, fig.align='center', fig.width=7, fig.height = 5, out.width='80%'----
oldPar <- par(mar = rep(0, 4))
plot(tree2)
ape::nodelabels()
ape::nodelabels(splits, as.integer(names(splits)), 
                adj = c(1.1, -0.2), cex = 0.8, frame = 'none')

## ---- echo=FALSE--------------------------------------------------------------
par(oldPar)

