\name{YULE}
\alias{YULE}
\title{Random Yule Trees}
\description{
  This benchmark simulated random phylogenies under the Yule model
  (i.e., without extinction) with the function
  \code{\link[ape]{rlineage}}. The process is repeated \code{N} times
  and the number of species is extracted and the frequencies of these
  values are compared with the expected values calculated with
  \code{\link[ape]{dbd}}.

  This benchmark is quite sensitive to the parameter values. For
  instance, \code{YULE(N = 10)} fails most of the time.

  Because the number of species can grow to a very large number if
  \code{lambda} and/or \code{Tmax} are large, this may result in an
  error if the simulated tree has more than 100,000 branches.
}
\usage{
YULE(N = 1000, lambda = 0.05, Tmax = 50, threshold = c(0.8, 1.2))
}
\arguments{
  \item{N}{the number of simulated trees.}
  \item{lambda}{the value the of speciation rate.}
  \item{Tmax}{the time span of the simulation.}
  \item{threshold}{the lower and upper bounds when comparing the
    observed and predicted numbers of species.}
}
\author{Emmanuel Paradis}
\references{
  Kendall, D. G. (1948) On the generalized ``birth-and-death''
  process. \emph{Annals of Mathematical Statistics}, \bold{19}, 1--15.

  Yule, G. U. (1924) A mathematical theory of evolution, based on the
  conclusions of Dr. J. C. Willis, F.R.S.. \emph{Philosophical
    Transactions of the Royal Society of London. Series B}, \bold{213},
  21--87.
}
\keyword{utilities}

