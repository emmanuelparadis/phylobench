\name{MANTEL}
\alias{MANTEL}
\title{Type I Error Rate of the Mantel Test}
\description{
  This benchmark assesses the type I error rate of the Mantel test
  computed by the function \code{\link[ape]{mantel.test}}. The test is
  performed on two random square matrices of size \code{n} and is
  replicated \code{N} times. If more than 5\% of the tests are
  significant, this cycle (simulation + test) is repeated; otherwise
  ``OK'' is returned. If after 200 cycles, there are more than 5\% of
  significant tests, an appropriate message is returned.

  If at least one missing value (NA or NaN) is observed, an appropriate
  message is returned.
}
\usage{
MANTEL(N = 100, n = 10)
}
\arguments{
  \item{N}{the number of replications.}
  \item{n}{the number of rows in the matrices.}
}
\author{Emmanuel Paradis}
\references{
  Mantel, N. (1967) The detection of disease clustering and a
  generalized regression approach. \emph{Cancer Research}, \bold{27},
  209--220.
}
\keyword{utilities}
