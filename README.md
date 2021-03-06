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
- Reordering of the edge matrix

New benchmarks are added easily by:

- Writing a function, say `FUN`, which runs the benchmark and returns `"OK"` if the results are as expected. This function must have no argument or its arguments defined by default so that it can be called with `FUN()`.
- Modifying the list `.list_of_tests` (in the file phylobench/R/phylobench.R) by adding the new benchmark like this:
```r
.list_of_tests <- list(...., "Title of the benchmark" = "FUN")
```
- Optionally, a file phylobench/man/FUN.Rd describing the new benchmark may be created.

If the benchmark requires files, these must be placed in phylobench/inst/extdata/.

All the benchmarks are run when building the vignette, so they can be visualized once the package is installed with:

```r
vignette("PhylogeneticBenchmarks")
```
