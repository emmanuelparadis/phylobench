## phylobench: Phylogenetic Benchmarking

phylobench implements tests of phylogenetic analyses and compares the output with expected outputs.

Currently, the implemented benchmarks are:

- Branching times calculation
- Base frequencies from DNA sequences
- Phylogenetically independent contrasts
- Variance-covariance under Brownian motion
- Neiborgh-joining
- Random coalescent trees
- Random Yule trees
- Type I error rate of the Mantel test
- Ultrametric trees
- Topological distances
- Splits from unrooted trees

New benchmarks are added easily by:

- Writing a function, say `FUN`, which runs the benchmark and returns `"OK"` if the result is as expected. This function must have no argument or its arguments defined by default so that it can be called with `FUN()`.
- Modify the list `.list_of_tests` (in the file phylobench/R/phylobench.R) by adding the new benchmark like this:
```r
.list_of_tests <- list(...., "Title of the benchmark" = "FUN")
```
- Optionally, create a file phylobench/man/FUN.Rd describing the new benchmark.
