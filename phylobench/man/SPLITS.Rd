\name{SPLITS}
\alias{SPLITS}
\title{Splits from Unrooted Trees}
\description{
  This benchmark assesses the extraction of splits (or bipartitions)
  from unrooted trees. The test trees are the three possible unrooted
  topologies with four tips. The splits are extracted with the functions
  \code{\link[ape]{prop.part}} and \code{\link[ape]{bitsplits}} both
  from \pkg{ape}. It is expected that three splits are observed in
  relative frequencies 1/3.
}
\usage{
SPLITS()
}
\author{Emmanuel Paradis}
\seealso{\code{\link{TOPODIST}}}
\keyword{utilities}
