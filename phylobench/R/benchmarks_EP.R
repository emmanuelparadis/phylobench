## benchmarks_EP.R (2019-02-06)

##   Phylogenetic Benchmarks

## Copyright 2019 Emmanuel Paradis

## This file is part of the R-package `phylobench'.
## See the file ../COPYING for licensing issues.

## read a tree, calculate its branching times, and compare them to values in a file
BTIMES <- function()
{
    nwk1 <- system.file("extdata/input/Newick/tree1_Newick.tre",
                        package = "phylobench")
    tr1 <- read.tree(nwk1)
    bt1 <- branching.times(tr1)
    bt1.0 <- scan(system.file("extdata/output/bt1.txt", package = "phylobench"),
                  sep = "\n", quiet = TRUE)
    if (all(abs(bt1 - bt1.0) < eps)) "OK"
    else "problem in branching times calculcation"
}

## base frequencies
BF <- function()
{
    fas1 <- system.file("extdata/input/FASTA/seq1_DNA.fas", package = "phylobench")
    dna1 <- read.dna(fas1, format = "f")
    BF1 <- base.freq(dna1, TRUE, TRUE)
    out1 <- system.file("extdata/output/BF1.txt", package = "phylobench")
    BF1.0 <- read.table(out1, header = TRUE)
    if (all(BF1[c("a", "c", "g", "t", "n")] == BF1.0)) "OK" else "problem when clculating base frequencies"
}

## Phylogenetically independent contrasts
PIC <- function()
{
    treefile <- system.file("extdata/input/Newick/tree_primates.tre", package = "phylobench")
    datfile <- system.file("extdata/input/Table/data_primates.txt", package = "phylobench")
    tree.primates <- read.tree(treefile)
    DATA <- read.table(datfile, header = TRUE)
    pic.body <- pic(DATA$body, tree.primates)
    pic.brain <- pic(DATA$brain, tree.primates)
    outfile <- system.file("extdata/output/PIC_primates.txt", package = "phylobench")
    PIC.0 <- read.table(outfile, header = TRUE)
    ## only 6 digits in PHYLIP's output
    test1 <- all(abs(sort(pic.body) - sort(PIC.0$body)) < 1e-5)
    test2 <- all(abs(sort(pic.brain) - sort(PIC.0$brain)) < 1e-5)
    if (test1 && test2) return("OK")
    return("disagreement between the values of PICs")
}

## Variance-covariance under Brownian motion
VCVBM <- function()
{
    tr <- compute.brtime(stree(5, "l"), 4:1)
    vcvape <- vcv(tr)
    expected.vcv <- diag(4, 5, 5)
    expected.vcv[lower.tri(expected.vcv)] <- offdiag <- rep(0:3, 4:1)
    expected.vcv <- t(expected.vcv)
    expected.vcv[lower.tri(expected.vcv)] <- offdiag
    if (all(expected.vcv == vcvape)) return("OK")
}

## Neiborgh-joining
NJ_SaitouNei <- function()
{
    matfile <- system.file("extdata/input/Table/M_SaitouNei.txt", package = "phylobench")
    M <- as.matrix(read.table(matfile))
    tr.nj <- nj(M)
    outfile <- system.file("extdata/output/tree_NJ_SaitouNei.tre", package = "phylobench")
    tr.ref <- read.tree(outfile)
    test <- all.equal(tr.nj, tr.ref)
    if (test) return("OK")
    return("disagreement between the reconstructed and reference NJ trees")
}

## Random coalescent trees
RCOAL <- function()
{
    BOUND <- qnorm(0.995)
    N <- 100
    tree.sizes <- c(5, 10, 20, 50, 75, 100)
    res <- numeric()
    for (i in 1:200) {
        for (n in tree.sizes) {
            k <- 2:n
            expected.mean <- 2 * sum(1/(k * (k - 1)))
            expected.var <- 4 * sum(1/(k * (k - 1))^2)
            x <- replicate(N, branching.times(rcoal(n))[1])
            res <- c(res, (mean(x) - expected.mean) * sqrt(N/expected.var))
        }
        if (anyNA(res))
            return(paste("some missing values returned after",
                         length(res), "simulations"))
        tab <- tabulate((abs(res) > BOUND) + 1L, 2L)[2]
        if (tab < length(res)/100) return("OK")
    }
    paste("found", tab, "replications", "out of range out of",
          length(res), "(1% expected)")
}

## Random Yule trees
YULE <- function(N = 1000, lambda = 0.05, Tmax = 50, threshold = c(0.8, 1.2))
{
    x <- replicate(floor(N/2), balance(rlineage(lambda, 0, Tmax))[1, ])
    dim(x) <- NULL
    mx <- max(x)
    O <- tabulate(x, mx)
    P <- length(x) * dyule(1:mx, lambda, Tmax)
    r <- cor(P, O)
    if (r < threshold[1] || r > threshold[2])
        return(paste("observed and predicted numbers of species seem too different: cor =", round(r, 3)))
    "OK"
}

## Type I error rate of the Mantel test
MANTEL <- function(N = 100, n = 10)
{
    rmat <- function(n) {
        x <- runif(n * (n - 1) / 2)
        m <- matrix(0, n, n)
        m[lower.tri(m)] <- x
        m <- t(m)
        m[lower.tri(m)] <- x
        m
    }
    res <- numeric()
    for (i in 1:200) {
        res <- c(res, replicate(N, {
            ma <- rmat(n)
            mb <- rmat(n)
            mantel.test(ma, mb)$p
            }))
        if (anyNA(res))
            return(paste("some missing values returned after",
                         length(res), "replications"))
        nsig <- sum(res < 0.05)
        if (nsig / length(res) <= 0.05) return("OK")
    }
    paste("number of significant tests seem too high after",
          length(res), "replications")
}

## Ultrametric trees
ULTRAMETRIC <- function(N = 100, n = c(5, 10, 20, 50, 100))
{
    res <- logical()
    for (k in n) {
        res <- c(res, !replicate(N, is.ultrametric(rcoal(k))))
        res <- c(res, replicate(N, is.ultrametric(rtree(k))))
    }
    if (any(res))
        return(paste(sum(res), "test(s) incorrect out of", length(res)))
    "OK"
}

