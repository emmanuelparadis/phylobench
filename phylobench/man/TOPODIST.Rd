\name{TOPODIST}
\alias{TOPODIST}
\title{Topological Distances}
\description{
  This benchmark assesses the topological distances between pairs of
  unrooted trees. Three trees are considered which represent the three
  possible unrooted topologies with three tips. The distances are
  calculated with the function \code{\link[ape]{dist.topo}}. If all
  distances are equal to two, ``OK'' is returned.
}
\usage{
TOPODIST()
}
\author{Emmanuel Paradis}
\references{
  Saitou, N. and Nei, M. (1987) The neighbor-joining method: a new
  method for reconstructing phylogenetic trees. \emph{Molecular Biology
    and Evolution}, \bold{4}, 406--425.
}
\seealso{\code{\link{SPLITS}}}
\keyword{utilities}
