\name{VCVBM}
\alias{VCVBM}
\title{Phylogenetic Variance-Covariance Matrix Under Brownian Motion Model}
\description{
  This benchmark assesses the variance-covariance (VCV) matrix for a
  trait evolving under Brownian motion on an ultrametric tree with five
  tips and four branching times equal to four, three, two, and one unit
  of time (from the root to the most recent node). In this case, the
  values in the VCV matrix can be calculated by hand. These values are
  compared with the values returned by the function \code{\link{vcv}}.
}
\usage{
VCVBM()
}
\author{Emmanuel Paradis}
\references{
  Felsenstein, J. (1985) Phylogenies and the comparative method.
  \emph{American Naturalist}, \bold{125}, 1--15.

  Martins, E. P. and Hansen, T. F. (1997) Phylogenies and the comparative
  method: a general approach to incorporating phylogenetic information
  into the analysis of interspecific data. \emph{American Naturalist},
  \bold{149}, 646--667.
}
\examples{
## the tree used in the benchmark:
tr <- compute.brtime(stree(5, "l"), 4:1)
plot(tr)
}
\keyword{utilities}
