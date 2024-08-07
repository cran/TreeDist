% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree_distance_nye.R
\encoding{UTF-8}
\name{JaccardRobinsonFoulds}
\alias{JaccardRobinsonFoulds}
\alias{JaccardSplitSimilarity}
\title{Jaccard–Robinson–Foulds metric}
\usage{
JaccardRobinsonFoulds(
  tree1,
  tree2 = NULL,
  k = 1L,
  allowConflict = TRUE,
  similarity = FALSE,
  normalize = FALSE,
  reportMatching = FALSE
)

JaccardSplitSimilarity(
  splits1,
  splits2,
  nTip = attr(splits1, "nTip"),
  k = 1L,
  allowConflict = TRUE,
  reportMatching = FALSE
)
}
\arguments{
\item{tree1, tree2}{Trees of class \code{phylo}, with leaves labelled identically,
or lists of such trees to undergo pairwise comparison.  Where implemented,
\code{tree2 = NULL} will compute distances between each pair of trees in the list
\code{tree1} using a fast algorithm based on
\insertCite{Day1985;textual}{TreeDist}.}

\item{k}{An arbitrary exponent to which to raise the Jaccard index.
Integer values greater than one are anticipated by Böcker \emph{et al}.
The Nye \emph{et al}. metric uses \code{k = 1}.
As k increases towards infinity, the metric converges to the Robinson–Foulds
metric.}

\item{allowConflict}{Logical specifying whether to allow conflicting splits
to be paired. If \code{FALSE}, such pairings will be allocated a similarity
score of zero.}

\item{similarity}{Logical specifying whether to report the result as a tree
similarity, rather than a difference.}

\item{normalize}{If a numeric value is provided, this will be used as a
maximum value against which to rescale results.
If \code{TRUE}, results will be rescaled against a maximum value calculated from
the specified tree sizes and topology, as specified in the "Normalization"
section below.
If \code{FALSE}, results will not be rescaled.}

\item{reportMatching}{Logical specifying whether to return the clade
matchings as an attribute of the score.}

\item{splits1, splits2}{Logical matrices where each row corresponds to a leaf,
either listed in the same order or bearing identical names (in any sequence),
and each column corresponds to a split, such that each leaf is identified as
a member of the ingroup (\code{TRUE}) or outgroup (\code{FALSE}) of the respective
split.}

\item{nTip}{(Optional) Integer specifying the number of leaves in each split.}
}
\value{
\code{JaccardRobinsonFoulds()} returns an array of numerics providing the
distances between each  pair of trees in \code{tree1} and \code{tree2},
or \code{splits1} and \code{splits2}.
}
\description{
Calculate the
\href{https://ms609.github.io/TreeDist/articles/Generalized-RF.html#jaccard-robinson-foulds-metric}{Jaccard–Robinson–Foulds metric}
\insertCite{Bocker2013}{TreeDist}, a
\href{https://ms609.github.io/TreeDist/articles/Robinson-Foulds.html#generalized-robinson-foulds-distances}{Generalized Robinson–Foulds metric}.
}
\details{
In short, the Jaccard–Robinson–Foulds
metric is a generalized Robinson-Foulds
metric: it finds the optimal matching that pairs each split in one tree with
a similar split in the second.
Matchings are scored according to the size of the largest split that is
consistent with both of them, normalized against the Jaccard index, and
raised to an arbitrary exponent.
A more detailed explanation is provided in the
\href{https://ms609.github.io/TreeDist/articles/Generalized-RF.html#jaccard-robinson-foulds-metric}{vignettes}.

By default, conflicting splits may be paired.

Note that the settings \verb{k = 1, allowConflict = TRUE, similarity = TRUE}
give the similarity metric of \insertCite{Nye2006;textual}{TreeDist};
a slightly faster implementation of this metric is available as
\code{\link[=NyeSimilarity]{NyeSimilarity()}}.

The examples section below details how to visualize matchings with
non-default parameter values.

Trees need not contain identical leaves; scores are based on the leaves that
trees hold in common.  Check for unexpected differences in tip labelling
with \code{setdiff(TipLabels(tree1), TipLabels(tree2))}.
}
\section{Normalization}{


If \code{normalize = TRUE}, then results will be rescaled from zero to one by
dividing by the maximum possible value for trees of the given topologies,
which is equal to the sum of the number of splits in each tree.
You may wish to normalize instead against the maximum number of splits
present in a pair of trees with \emph{n} leaves, by specifying
\code{normalize = n - 3}.
}

\examples{
set.seed(2)
tree1 <- ape::rtree(10)
tree2 <- ape::rtree(10)
JaccardRobinsonFoulds(tree1, tree2, k = 2, allowConflict = FALSE)
JaccardRobinsonFoulds(tree1, tree2, k = 2, allowConflict = TRUE)

JRF2 <- function(tree1, tree2, ...) 
  JaccardRobinsonFoulds(tree1, tree2, k = 2, allowConflict = FALSE, ...)
  
VisualizeMatching(JRF2, tree1, tree2, matchZeros = FALSE)
}
\references{
\insertAllCited{}
}
\seealso{
Other tree distances: 
\code{\link{KendallColijn}()},
\code{\link{MASTSize}()},
\code{\link{MatchingSplitDistance}()},
\code{\link{NNIDist}()},
\code{\link{NyeSimilarity}()},
\code{\link{PathDist}()},
\code{\link{Robinson-Foulds}},
\code{\link{SPRDist}()},
\code{\link{TreeDistance}()}
}
\author{
\href{https://orcid.org/0000-0001-5660-1727}{Martin R. Smith}
(\href{mailto:martin.smith@durham.ac.uk}{martin.smith@durham.ac.uk})
}
\concept{tree distances}
