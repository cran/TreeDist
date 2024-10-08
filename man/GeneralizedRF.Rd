% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree_distance.R
\encoding{UTF-8}
\name{GeneralizedRF}
\alias{GeneralizedRF}
\title{Generalized Robinson–Foulds distance}
\usage{
GeneralizedRF(
  splits1,
  splits2,
  nTip,
  PairScorer,
  maximize,
  reportMatching,
  ...
)
}
\arguments{
\item{splits1, splits2}{Logical matrices where each row corresponds to a leaf,
either listed in the same order or bearing identical names (in any sequence),
and each column corresponds to a split, such that each leaf is identified as
a member of the ingroup (\code{TRUE}) or outgroup (\code{FALSE}) of the respective
split.}

\item{nTip}{Integer specifying the number of leaves in each split.}

\item{PairScorer}{function taking four arguments, \code{splits1}, \code{splits2},
\code{nSplits1}, \code{nSplits2}, which should return the score of each pair of splits
in a two-dimensional matrix.  Additional parameters may be specified via
\dots.}

\item{reportMatching}{Logical specifying whether to return the clade
matchings as an attribute of the score.}

\item{\dots}{Additional parameters to \code{PairScorer}}
}
\value{
A numeric value specifying the score of the tree pairs under the
specified pair scorer. If \code{reportMatching = TRUE}, attribute also list:
\itemize{
\item \code{matching}: which split in \code{splits2} is optimally matched to each split in
\code{split1} (\code{NA} if not matched);
\item \code{matchedSplits}: Textual representation of each match
\item \code{matchedScores}: Scores for matched split.
\item \code{pairScores}: Calculated scores for each possible matching of each split.
}
}
\description{
An internal function to calculate Generalized Robinson–Foulds
distances from splits.
}
\details{
Note that no checks will be made to confirm that \code{splits1} and \code{splits2}
contain the same leaves in the same order.
This is the responsibility of the calling function.
}
\author{
\href{https://orcid.org/0000-0001-5660-1727}{Martin R. Smith}
(\href{mailto:martin.smith@durham.ac.uk}{martin.smith@durham.ac.uk})
}
\keyword{internal}
