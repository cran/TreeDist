% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{TreeDistPlot}
\alias{TreeDistPlot}
\title{Plot a simple tree}
\usage{
TreeDistPlot(
  tr,
  title = NULL,
  bold = NULL,
  leaveRoom = FALSE,
  prune = integer(0),
  graft = integer(0),
  edge.color = "black",
  edge.width = NULL,
  ...
)
}
\arguments{
\item{tr}{A tree of class \code{phylo}, with leaves labelled as integers.}

\item{title}{\code{main} title for the plot.}

\item{bold}{Integer specifying which leaves to print in bold.}

\item{leaveRoom}{Logical specifying whether to leave space to print
tree distances beneath the plotted tree.}

\item{prune, graft}{Integer vectors specifying which edges to highlight as
pruned and grafted.}

\item{edge.color}{Additional parameter to \code{plot.phylo}; will be overridden
by \code{prune} and \code{graft} as requested.}

\item{edge.width, \dots}{Additional parameters to \code{plot.phylo}.}
}
\description{
Convenience plotting function used in vignettes and documentation.
}
\author{
\href{https://orcid.org/0000-0001-5660-1727}{Martin R. Smith}
(\href{mailto:martin.smith@durham.ac.uk}{martin.smith@durham.ac.uk})
}
\keyword{internal}
