#' Kendall&ndash;Colijn distance
#' 
#' Calculate the Kendall&ndash;Colijn tree distance, a measure related to the 
#' path difference. 
#' 
#' The Kendall&ndash;Colijn distance works by measuring, for each pair of
#' leaves, the distance from the most recent common ancestor of those leaves
#' and the root node.
#' For a given tree, this produces a vector of values recording the
#' distance-from-the-root of each most recent common ancestor of each pair of
#' leaves.
#' 
#' Two trees are compared by taking the Euclidean distance between the
#' respective vectors.  This is calculated by taking the square root of the sum 
#' of the squares of the differences between the vectors.
#' 
#' An analogous distance can be created from any vector representation of a
#' tree.
#' The split size vector metric \insertCite{SmithSpace}{TreeDist} is an attempt
#' to mimic the Kendall&nbsp;Colijn metric in situations where the position of
#' the root should not be afforded special significance; and the path distance
#' \insertCite{Steel1993}{TreeDist} is a familiar alternative whose underlying
#' vector measures the distance of the last common ancestor of each pair
#' of leaves from the leaves themselves, i.e. the length of the path from one 
#' leaf to another.
#' 
#' 
#' None of these vector-based methods performs as well as other tree distances
#' in measuring similarities in the relationships implied by a pair of trees
#' \insertCite{SmithDist}{TreeDist}; in particular, the Kendall&nbsp;Colijn
#' metric is strongly influenced by tree balance, and may not be appropriate
#' for a suite of common applications \insertCite{SmithSpace}{TreeDist}.
#' 
#' @template tree12ListParams
#' @param Vector Function converting a tree to a numeric vector.
#' 
#' `KCVector`, the default, returns the number of edges between the common
#' ancestor of each pair of leaves and the root of the tree
#' \insertCite{@per @Kendall2016}{TreeDist}.
#' 
#' `PathVector` returns the number of edges between each pair of leaves
#' \insertCite{@per @Steel1993}{TreeDist}.
#' 
#' `SplitVector` returns the size of the smallest split that contains each
#' pair of leaves (per \insertCite{SmithSpace;nobrackets}{TreeDist}).
#' 
#' @templateVar returns `KendallColijn()` returns
#' @template distReturn
#' 
#' @examples 
#' KendallColijn(TreeTools::BalancedTree(8), TreeTools::PectinateTree(8))
#'
#' set.seed(0)
#' KendallColijn(TreeTools::BalancedTree(8), lapply(rep(8, 3), ape::rtree))
#' KendallColijn(lapply(rep(8, 4), ape::rtree))
#' 
#' KendallColijn(lapply(rep(8, 4), ape::rtree), Vector = SplitVector)
#' 
#' # Notice that changing tree shape close to the root results in much
#' # larger differences
#' tree1 <- ape::read.tree(text = "(a, (b, (c, (d, (e, (f, (g, h)))))));")
#' tree2 <- ape::read.tree(text = "(a, ((b, c), (d, (e, (f, (g, h))))));")
#' tree3 <- ape::read.tree(text = "(a, (b, (c, (d, (e, ((f, g), h))))));")
#' trees <- c(tree1, tree2, tree3)
#' KendallColijn(trees)
#' KendallColijn(trees, Vector = SplitVector)
#' @template MRS
#' 
#' @seealso [`treespace::treeDist`](https://CRAN.R-project.org/package=treespace/vignettes/introduction.html)
#' is a more sophisticated, if more cumbersome, implementation that supports 
#' lambda > 0, i.e. use of edge lengths in tree comparison.
#' 
#' @references \insertAllCited{}
#' 
#' @family tree distances
#' @importFrom utils combn
#' @encoding UTF-8
#' @export
KendallColijn <- function(tree1, tree2 = NULL, Vector = KCVector) {
  
  FunValue <- function(nTip) double(nTip * (nTip - 1L) / 2L)
  
  if (inherits(tree1, "phylo")) {
    if (inherits(tree2, "phylo")) {
      if (length(tree1[["tip.label"]]) != length(tree2[["tip.label"]]) || 
          length(setdiff(tree1[["tip.label"]], tree2[["tip.label"]])) > 0) {
        stop("Leaves must bear identical labels.")
      }
      .EuclideanDistance(Vector(tree1) - Vector(tree2))
    } else {
      if (is.null(tree2)) {
        0
      } else {
        apply(Vector(tree1) - vapply(tree2, Vector,
                                     FunValue(length(tree1[["tip.label"]]))),
              2L, .EuclideanDistance)
      }
    }
  } else {
    if (inherits(tree2, "phylo")) {
      apply(Vector(tree2) - vapply(tree1, Vector,
                                   FunValue(length(tree2[["tip.label"]]))),
            2L, .EuclideanDistance)
    } else if (is.null(tree2)) {
      
      treeVec <- vapply(tree1, Vector, FunValue(length(tree1[[1]][["tip.label"]])))
      nTree <- length(tree1)
      ret <- matrix(0, nTree, nTree)
      is <- combn(seq_len(nTree), 2)
      
      ret <- structure(class = "dist", Size = nTree,
                       Diag = FALSE, Upper = FALSE,
                       apply(is, 2, function(i)
                         .EuclideanDistance(treeVec[, i[[1]]] - treeVec[, i[[2]]])))
      # Return:
      ret
    }
      else {
      vector1 <- vapply(tree1, Vector, FunValue(length(tree1[[1]][["tip.label"]])))
      vector2 <- vapply(tree2, Vector, FunValue(length(tree2[[1]][["tip.label"]])))
      apply(vector2, 2, function(i) 
        apply(vector1, 2, function(j) 
          .EuclideanDistance(i - j)))
    }
  }
}

.EuclideanDistance <- function(x) sqrt(sum(x * x))

#' @describeIn KendallColijn Creates a vector that characterises a rooted tree,
#' as described in \insertCite{Kendall2016;textual}{TreeDist}.
#' @param tree A tree of class \code{\link[ape:read.tree]{phylo}}.
#' @importFrom TreeTools AllAncestors Preorder
#' @importFrom utils combn
#' @export
KCVector <- function(tree) {
  tree <- Preorder(tree)
  edge <- tree[["edge"]]
  parent <- edge[, 1L]
  child <- edge[, 2L]
  root <- parent[[1]]
  nTip <- root - 1L
  tipOrder <- order(tree[["tip.label"]])
  is <- combn(tipOrder, 2)
  
  ancestors <- AllAncestors(parent, child)
  
  mrca <- apply(is, 2, function(i) 
    max(intersect(ancestors[[i[[1]]]], ancestors[[i[[2]]]])))
  
  rootDist <- lengths(ancestors)
  structure(rootDist[mrca], Size = nTip, class = "dist")
}

#' @describeIn KendallColijn Creates a vector reporting the number of edges
#' between each pair of leaves, per the path metric of
#' \insertCite{Steel1993;textual}{TreeDist}.
#' @importFrom TreeTools AllAncestors Preorder
#' @importFrom utils combn
#' @export
PathVector <- function(tree) {
  if (!inherits(tree, "phylo")) {
    stop("`tree` must be of class `phylo`")
  }
  edge <- tree[["edge"]]
  structure(path_vector(edge), Size = NTip(tree), class = "dist")
}

#' @describeIn KendallColijn Creates a vector reporting the smallest split
#' containing each pair of leaves, per the metric proposed in
#' \insertCite{SmithSpace;textual}{TreeDist}.
#' @importFrom TreeTools as.Splits
#' @export
SplitVector <- function(tree) {
  tipLabel <- tree[["tip.label"]]
  nTip <- length(tipLabel)
  splits <- as.logical(as.Splits(tree, tipLabel[order(tipLabel)]))
  splits <- rbind(splits, !splits)
  inSplit <- rowSums(splits)
  is <- combn(seq_len(nTip), 2)
  
  smallestSplit <- apply(is, 2, function(i)
    min(inSplit[splits[, i[[1]]] & splits[, i[[2]]]], nTip - 1L))
  
  structure(smallestSplit, Size = nTip, class = "dist")
}

#' Approximate diameter of the Kendall&ndash;Collijn metric
#'
#' @return `KCDiameter()` returns the value of the Kendall & Colijn's (2016)
#' metric distance between two pectinate trees with _n_ leaves ordered in
#' the opposite direction, which I suggest (without any attempt at a proof) may
#' be a useful proxy for the diameter (i.e. maximum value) of the K&ndash;C
#' metric.
#'
#' @examples
#' KCDiameter(4)
#' KCDiameter(trees)
#' @importFrom TreeTools PectinateTree
#' @rdname KendallColijn
#' @export
KCDiameter <- function(tree) UseMethod("KCDiameter")

#' @importFrom TreeTools NTip
#' @export
KCDiameter.phylo <- function(tree) {
  KCDiameter.numeric(NTip(tree))
}

#' @export
KCDiameter.numeric <- function(tree) {
  nTip <- as.integer(tree)
  mat <- matrix(seq_len(nTip), nTip, nTip)
  Euclid <- function(x, y) sqrt(sum((x - y) ^ 2))
  
  # Return:
  Euclid(nTip - mat[lower.tri(mat)] + 1L, t(mat)[lower.tri(mat)])
}

#' @export
KCDiameter.list <- function(tree) {
  vapply(tree, KCDiameter, double(1))
}

#' @export
KCDiameter.multiPhylo <- function(tree) {
  vapply(tree, KCDiameter, double(1))
}
