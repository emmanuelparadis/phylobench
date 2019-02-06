\name{RCOAL}
\alias{RCOAL}
\title{Random Coalescent Trees}
\description{
This benchmark simulates random coalescent trees with the function
\code{\link[ape]{rcoal}}, calculates the time to the MRCA (root age)
with \code{\link[ape]{branching.times}} for different sample sizes
(\emph{n} = 5, 10, 20, 50, 75, 100). The expected mean and variance are
calculated with:

\deqn{2\sum_{i=2}^n\frac{1}{i(i-1)}}{2 \sum 1/(i(i-1)), i=2, \dots, n}

\deqn{4\sum_{i=2}^n\frac{1}{[i(i-1)]^2}}{4 \sum 1/(i(i-1))^2, i=2, \dots, n}

For each value of \emph{n}, 100 trees are simulated and the mean root
age is centered and scaled with the above formulas, so the transformed
values are expected to follow a standard normal distribution (assuming
the central limit theorem applies). If more than 1\% of these
transformed values are smaller or larger than expected as indicated by
the quantile of the normal distribution (see
\code{\link[stats]{qnorm}}), the whole process is repeated, otherwise
the test returns ``OK''. The maximum number of repetitions is 200
(equivalent to 120,000 simulated trees). If after these 200 repetitions,
more than 1\% out of the 1200 means are out of range, an appropriate
message is returned.

If at least one missing value (NA or NaN) is observed, an appropriate
message is returned.
}
\usage{
RCOAL()
}
\author{Emmanuel Paradis}
\references{
Hudson, R. R. (1991) Gene genealogies and the coalescent process. \emph{Oxford Surveys in Evolutionary Biology}, \bold{7}, 1--44.

Kingman, J. F. C. (1982) On the genealogy of large populations. \emph{Journal of Applied Probability}, \bold{19A}, 27--43.

Wakeley, J. (2009) Coalescent theory: an introduction. Roberts \& Company Publishers, Greenwood Village CO.
}
\keyword{utilities}

