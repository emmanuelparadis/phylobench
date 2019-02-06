\name{runTests}
\alias{runTests}
\alias{codeTests}
\alias{listTests}
\alias{fileTests}
\title{Phylogenetic Benchmarking}
\description{
  \code{runTests} runs a series of phylogenetic benchmark tests. The
  other functions are utilities.
}
\usage{
runTests(verbose = TRUE)
listTests()
fileTests(which)
codeTests(which)
}
\arguments{
  \item{verbose}{a logical value specifying whether to print the
    progress of the tests. Set it to \code{FALSE} for a completely
    silent testing.}
  \item{which}{a number giving the number of the test.}
}
\author{Emmanuel Paradis}
\value{
  \code{runTests} returns a named list with the results of the tests.

  \code{listTests} returns a data frame with the titles and functions of
  the tests.

  \code{fileTests} and \code{codeTests} simply prints their results.
}
\examples{
## This is a test to check that the functions match with
## those listed in the test list inside the package:
fun.db <- unlist(phylobench:::.list_of_tests)
fun.pkg <- ls(env = asNamespace("phylobench"))
del <- match(c("runTests", "listTests", "fileTests", "codeTests", "eps"),
             fun.pkg)
fun.pkg <- fun.pkg[-del]
test1 <- length(fun.db) == length(fun.pkg)
test2 <- all(sort(fun.db) == sort(fun.pkg))
if (!(test1 && test2)) {
    cat("Function(s) in the package not in the list:",
        fun.pkg[is.na(match(fun.pkg, fun.db))], "\n")
    cat("Function(s) in the list not in the package:",
        fun.db[is.na(match(fun.db, fun.pkg))], "\n")
    stop("Check the functions and the list of functions.")
}
}
\keyword{utilities}
