% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree_information.R
\encoding{UTF-8}
\name{TreeInfo}
\alias{TreeInfo}
\alias{SplitwiseInfo}
\alias{ClusteringEntropy}
\alias{ClusteringInfo}
\alias{ClusteringEntropy.phylo}
\alias{ClusteringEntropy.list}
\alias{ClusteringEntropy.multiPhylo}
\alias{ClusteringEntropy.Splits}
\alias{ClusteringInfo.phylo}
\alias{ClusteringInfo.list}
\alias{ClusteringInfo.multiPhylo}
\alias{ClusteringInfo.Splits}
\alias{ConsensusInfo}
\title{Information content of splits within a tree}
\usage{
SplitwiseInfo(x, p = NULL, sum = TRUE)

ClusteringEntropy(x, p = NULL, sum = TRUE)

ClusteringInfo(x, p = NULL, sum = TRUE)

\method{ClusteringEntropy}{phylo}(x, p = NULL, sum = TRUE)

\method{ClusteringEntropy}{list}(x, p = NULL, sum = TRUE)

\method{ClusteringEntropy}{multiPhylo}(x, p = NULL, sum = TRUE)

\method{ClusteringEntropy}{Splits}(x, p = NULL, sum = TRUE)

\method{ClusteringInfo}{phylo}(x, p = NULL, sum = TRUE)

\method{ClusteringInfo}{list}(x, p = NULL, sum = TRUE)

\method{ClusteringInfo}{multiPhylo}(x, p = NULL, sum = TRUE)

\method{ClusteringInfo}{Splits}(x, p = NULL, sum = TRUE)

ConsensusInfo(trees, info = "phylogenetic", p = 0.5, check.tips = TRUE)
}
\arguments{
\item{x}{A tree of class \code{phylo}, a list of trees, or a \code{multiPhylo} object.}

\item{p}{Scalar from 0.5 to 1 specifying minimum proportion of trees that
must contain a split for it to appear within the consensus.}

\item{sum}{Logical: if \code{TRUE}, sum the information content of each split to
provide the total splitwise information content of the tree.}

\item{trees}{List of \code{phylo} objects, optionally with class \code{multiPhylo}.}

\item{info}{Abbreviation of "phylogenetic" or "clustering", specifying
the concept of information to employ.}

\item{check.tips}{Logical specifying whether to renumber leaves such that
leaf numbering is consistent in all trees.}
}
\value{
\code{SplitwiseInfo()}, \code{ClusteringInfo()} and \code{ClusteringEntropy()}
return the splitwise information content of the tree -- or of each split
in turn, if \code{sum = FALSE} -- in bits.

\code{ConsensusInfo()} returns the splitwise information content of the
majority rule consensus of \code{trees}.
}
\description{
Sum the entropy (\code{ClusteringEntropy()}), clustering information content
(\code{ClusteringInfo()}), or phylogenetic information content (\code{SplitwiseInfo()})
across each split within a phylogenetic tree,
or the consensus of a set of phylogenetic trees (\code{ConsensusInfo()}).
This value will be greater than the total information
content of the tree where a tree contains multiple splits, as
these splits are not independent and thus contain mutual information that is
counted more than once
}
\section{Clustering information}{


Clustering entropy addresses the question "how much information is contained
in the splits within a tree". Its approach is complementary to the
phylogenetic information content, used in \code{\link[=SplitwiseInfo]{SplitwiseInfo()}}.
In essence, it asks, given a split that subdivides the leaves of a tree into
two partitions, how easy it is to predict which partition a randomly drawn
leaf belongs to \insertCite{@Meila2007; @Vinh2010}{TreeDist}.

Formally, the entropy of a split \emph{S} that divides \emph{n} leaves into two
partitions of sizes \emph{a} and \emph{b} is given by
\emph{H(S)} = - \emph{a/n} log \emph{a/n} - \emph{b/n} log \emph{b/n}.

Base 2 logarithms are conventionally used, such that entropy is measured in
bits.
Entropy denotes the number of bits that are necessary to encode the outcome
of a random variable: here, the random variable is "what partition does a
randomly selected leaf belong to".

An even split has an entropy of 1 bit: there is no better way of encoding
an outcome than using one bit to specify which of the two partitions the
randomly selected leaf belongs to.

An uneven split has a lower entropy: membership of the larger partition is
common, and thus less surprising; it can be signified using fewer bits in an
optimal compression system.

If this sounds confusing, let's consider creating a code to transmit the
cluster label of two randomly selected leaves.  One straightforward
option would be to use
\itemize{
\item \code{00} = "Both leaves belong to partition A"
\item \code{11} = "Both leaves belong to partition B"
\item \code{01} = 'First leaf in A, second in B`
\item \code{10} = 'First leaf in B, second in A`
}

This code uses two bits to transmit the partition labels of two leaves.
If partitions A and B are equiprobable, this is the optimal code; our
entropy -- the average information content required per leaf -- is 1 bit.

Alternatively, we could use the (suboptimal) code
\itemize{
\item \code{0} = "Both leaves belong to partition A"
\item \code{111} = "Both leaves belong to partition B"
\item \code{101} = 'First leaf in A, second in B`
\item \code{110} = 'First leaf in B, second in A`
}

If A is much larger than B, then most pairs of leaves will require just
a single bit (code \code{0}). The additional bits when 1+ leaf belongs to B
may be required sufficiently rarely that the average message
requires fewer than two bits for two leaves, so the entropy is less than
1 bit.  (The optimal coding strategy will depend on the exact sizes
of A and B.)

As entropy measures the bits required to transmit the cluster label of each
leaf \insertCite{@Vinh2010: p. 2840}{TreeDist}, the information content of
a split is its entropy multiplied by the number of leaves.
}

\section{Phylogenetic information}{


Phylogenetic information expresses the information content of a split
in terms of the probability that a uniformly selected tree will contain it
\insertCite{Thorley1998}{TreeDist}.
}

\section{Consensus information}{


The information content of splits in a consensus tree is calculated by
interpreting support values (i.e. the proportion of trees containing
each split in the consensus) as probabilities that the true tree contains
that split, following \insertCite{SmithCons;textual}{TreeDist}.
}

\examples{
library("TreeTools", quietly = TRUE)

SplitwiseInfo(PectinateTree(8))
tree <- read.tree(text = "(a, b, (c, (d, e, (f, g)0.8))0.9);")
SplitwiseInfo(tree)
SplitwiseInfo(tree, TRUE)

# Clustering entropy of an even split = 1 bit
ClusteringEntropy(TreeTools::as.Splits(c(rep(TRUE, 4), rep(FALSE, 4))))

# Clustering entropy of an uneven split < 1 bit
ClusteringEntropy(TreeTools::as.Splits(c(rep(TRUE, 2), rep(FALSE, 6))))

tree1 <- TreeTools::BalancedTree(8)
tree2 <- TreeTools::PectinateTree(8)

ClusteringInfo(tree1)
ClusteringEntropy(tree1)
ClusteringInfo(list(one = tree1, two = tree2))

ClusteringInfo(tree1) + ClusteringInfo(tree2)
ClusteringEntropy(tree1) + ClusteringEntropy(tree2)
ClusteringInfoDistance(tree1, tree2)
MutualClusteringInfo(tree1, tree2)

# Clustering entropy with uncertain splits
tree <- ape::read.tree(text = "(a, b, (c, (d, e, (f, g)0.8))0.9);")
ClusteringInfo(tree)
ClusteringInfo(tree, TRUE)

# Support-weighted information content of a consensus tree
set.seed(0)
trees <- list(RandomTree(8), RootTree(BalancedTree(8), 1), PectinateTree(8))
cons <- consensus(trees, p = 0.5)
p <- SplitFrequency(cons, trees) / length(trees)
plot(cons)
LabelSplits(cons, signif(SplitwiseInfo(cons, p, sum = FALSE), 4))
ConsensusInfo(trees)
LabelSplits(cons, signif(ClusteringInfo(cons, p, sum = FALSE), 4))
ConsensusInfo(trees, "clustering")
}
\references{
\insertAllCited{}
}
\seealso{
An introduction to the phylogenetic information content of a split is given
in \href{https://ms609.github.io/TreeTools/reference/SplitInformation.html}{\code{SplitInformation()}}
and in a \href{https://ms609.github.io/TreeDist/articles/information.html}{package vignette}.

Other information functions: 
\code{\link{SplitEntropy}()},
\code{\link{SplitSharedInformation}()}
}
\author{
\href{https://orcid.org/0000-0001-5660-1727}{Martin R. Smith}
(\href{mailto:martin.smith@durham.ac.uk}{martin.smith@durham.ac.uk})
}
\concept{information functions}
