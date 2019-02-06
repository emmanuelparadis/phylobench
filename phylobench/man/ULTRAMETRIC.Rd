\name{ULTRAMETRIC}
\alias{ULTRAMETRIC}
\title{Ultrametric Trees}
\description{
  For each value in \code{n}, one tree is simulated with
  \code{\link[ape]{rtree}} and one tree with
  \code{\link[ape]{rcoal}}. Both trees are tested with
  \code{\link[ape]{is.ultrametric}}. This is replicated \code{N}
  times. If all tests are as expected, ``OK'' is returned; otherwise,
  a message with the number of unexpected results is returned.
}
\usage{
ULTRAMETRIC(N = 100, n = c(5, 10, 20, 50, 100))
}
\arguments{
  \item{N}{the number of replications.}
  \item{n}{the tree sizes.}
}
\author{Emmanuel Paradis}
\keyword{utilities}
