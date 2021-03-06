\name{REORDERPHYLO}
\alias{REORDERPHYLO}
\title{Test Reordering of Edge Matrix}
\description{
  This benchmark tests whether reodering the edges of a \code{"phylo"}
  tree works correctly.

  This benchmark is quite critical as reordering the edges of a
  \code{"phylo"} tree is an important operation used in many functions.
}
\usage{
REORDERPHYLO(Nmin = 3L, Nmax = 1000L, ProbRooted = 0.5,
             ProbMultichotomy = 0.5, nrep = 1e4L)
}
\arguments{
  \item{Nmin, Nmax}{the smallest and largest values allowed for the
    number of tips, by default between 3 and 1000.}
  \item{ProbRooted}{the probability that the tree is simulated rooted;
    0.5 by default.}
  \item{ProbMultichotomy}{the probability that multichotomies are
    introduced into the tree, 0.5 by default.}
  \item{nrep}{the number of trees simulated (10,000 by default).}
}
\details{
  The edges (or branches) in a tree of class \code{"phylo"} can be
  ordered in cladewise, pruningwise, or postorder order (see the
  definition of the class on \pkg{ape}'s web site for mode details).

  The idea of this benchmark is to simulate a tree with the function
  \code{\link[ape]{rtree}} which outputs trees in cladewise order. The
  simulated tree can be rooted or not and include multichotomies, both
  being controlled by the above options. The simulated tree is then
  reordered into pruningwise order, and back into cladewise order: the
  edge matrix of the final tree is expected to be identical to the one
  of the original tree. The same operation is performed with the
  postorder order. Both operations are repeated \code{nrep} on a newly
  simulated tree.

  If at least one test does not run as expected, a data frame is
  returned with information on all simulated trees (see below).
}
\value{
  \code{"OK"} if no problem is detected; otherwise, a data frame is
  returned with \code{nrep} rows and the following columns:

  \itemize{
    \item{Ntip}{the number of tips in the simulated tree.}
    \item{Nnode}{the number of nodes in the simulated tree}
    \item{Rooted}{a logical value indicating whether the tree was rooted
      or not.}
    \item{Test1}{a logical value indicating whether the first test
      described above run as expected.}
    \item{Test2}{idem for the second test.}
  }
}
\author{Emmanuel Paradis}
\keyword{utilities}
