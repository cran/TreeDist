#' Matching Split Distance
#' 
#' Calculate the 
#' [Matching Split Distance](https://ms609.github.io/TreeDist/articles/Generalized-RF.html#matching-split-distance)
#' \insertCite{Bogdanowicz2012,Lin2012}{TreeDist} for unrooted binary trees.
#' 
#' Trees need not contain identical leaves; scores are based on the leaves that
#' trees hold in common.  Check for unexpected differences in tip labelling
#' with `setdiff(TipLabels(tree1), TipLabels(tree2))`.
#' 
#' @inheritParams TreeDistance
#' 
#' @templateVar returns `MatchingSplitDistance()` returns
#' @template distReturn
#' 
#' @section Normalization:
#' 
#' A normalization value or function must be provided in order to return a
#' normalized value.  If you are aware of a generalised formula, please
#' let me know by
#' \href{https://github.com/ms609/TreeDist/issues/new}{creating a GitHub issue}
#' so that it can be implemented.
#' 
#' @examples 
#' MatchingSplitDistance(lapply(rep(8, 5), ape::rtree), normalize = 16)
#' 
#' MatchingSplitDistance(TreeTools::BalancedTree(6),
#'                       TreeTools::PectinateTree(6),
#'                       reportMatching = TRUE)
#' 
#' VisualizeMatching(MatchingSplitDistance, TreeTools::BalancedTree(6),
#'                   TreeTools::PectinateTree(6))
#' @template MRS
#'  
#' @references \insertAllCited{}
#' 
#' @family tree distances
#' 
#' @export
MatchingSplitDistance <- function(tree1, tree2 = NULL, normalize = FALSE,
                                   reportMatching = FALSE) {
  unnormalized <- CalculateTreeDistance(MatchingSplitDistanceSplits, tree1,
                                        tree2, reportMatching)
  
  # Return:
  NormalizeInfo(unnormalized, tree1, tree2, how = normalize,
                InfoInTree = function(X) stop(
                  "Please specify a function to generate a normalizing constant"
                  ),
                Combine = max)
}

#' @rdname MatchingSplitDistance
#' @inheritParams SharedPhylogeneticInfoSplits
#' @useDynLib TreeDist, .registration = TRUE
#' @export
MatchingSplitDistanceSplits <- function(splits1, splits2, 
                                         nTip = attr(splits1, "nTip"),
                                         normalize = TRUE, 
                                         reportMatching = FALSE) {
  GeneralizedRF(splits1, splits2, nTip, cpp_matching_split_distance,
                maximize = FALSE, reportMatching = reportMatching)
}
