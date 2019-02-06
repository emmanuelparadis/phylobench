\name{PIC}
\alias{PIC}
\title{Phylogenetically Independent Contrasts Benchmark}
\description{
  This benchmark compares the phylogenetically independent contrasts
  calculated by \pkg{ape} with those output from Phylip (Felsenstein
  2004). The data are a phylogeny of five species of primates and two
  variables: body mass and brain mass, both log-transformed.

  The comparisons are done to the nearest 1e-5 since the results
  reported in Phylip have five digits.
}
\usage{
PIC()
}
\author{Emmanuel Paradis}
\references{
  Felsenstein, J. (2004) Phylip (Phylogeny Inference Package) version
  3.68. Department of Genetics, University of Washington, Seattle, USA.
  \url{http://evolution.genetics.washington.edu/phylip/phylip.html}.
}
\keyword{utilities}
