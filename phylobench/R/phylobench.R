## phylobench.R (2021-04-20)

##   Phylogenetic Benchmarks

## Copyright 2019-2021 Emmanuel Paradis

## This file is part of the R-package `phylobench'.
## See the file ../COPYING for licensing issues.

.list_of_tests <- list("Branching times calculation" = "BTIMES",
                       "Base frequencies from DNA sequences" = "BF",
                       "Phylogenetically independent contrasts" = "PIC",
                       "Variance-covariance under Brownian motion" = "VCVBM",
                       "Neiborgh-joining" = "NJ_SaitouNei",
                       "Random coalescent trees" = "RCOAL",
                       "Random Yule trees" = "YULE",
                       "Type I error rate of the Mantel test" = "MANTEL",
                       "Ultrametric trees" = "ULTRAMETRIC",
                       "Topological distances" = "TOPODIST",
                       "Splits from unrooted trees" = "SPLITS",
                       "Test reordering of edge matrix" = "REORDERPHYLO")

eps <- .Machine$double.eps

runTests <- function(verbose = TRUE)
{
    FUN <- .list_of_tests
    tl <- names(FUN)
    ntests <- length(tl)
    res <- vector("list", ntests)
    if (verbose)
        cat("    Phylogenetic benchmarking: starting", ntests, "tests...\n\n")
    for (i in 1:ntests) {
        if (verbose) cat("Running test no.", i, ": ", tl[i], "...", sep = "")
        out <- try(eval(parse(text = paste0(FUN[[i]], "()"))))
        if (identical(out, "OK")) cat(" OK.\n")
        else cat(" problem!\n")
        res[[i]] <- out
    }
    if (verbose) {
        allok <- sapply(res, identical, y = "OK")
        Nnotok <- sum(!allok)
        if (!Nnotok) cat("\nAll tests were OK.\n")
        else cat("\n", Nnotok, "tests out of", ntests, "were not OK. See details in the returned list.\n")
    }
    names(res) <- tl
    res
}

listTests <- function()
{
    FUN <- .list_of_tests
    DF <- data.frame(Function = unlist(FUN))
    row.names(DF) <- paste(1:nrow(DF), names(FUN), sep = ": ")
    DF
}

fileTests <- function(which)
{
    if (missing(which))
        stop("give the number of the test to show its code (see listTests())")
    FUN <- .list_of_tests[which]
    tl <- names(FUN)
    code <- deparse(get(FUN[[1]]))
    files <- grep("system\\.file", code, value = TRUE)
    prefix <- paste0(system.file(package = "phylobench"), "/")
    files <- gsub(".*system\\.file\\(\"", prefix, files)
    files <- gsub("\".*$", "", files)
    infiles <- grep("/extdata/input/", files, value = TRUE)
    outfiles <- grep("/extdata/output/", files, value = TRUE)
    cat("Test no.", which, ": ", tl, ":\n\n", sep = "")
    cat("Input (data) files:", infiles, sep = "\n")
    cat("\nOutput (result) files:", outfiles, sep = "\n")
}

codeTests <- function(which)
{
    if (missing(which))
        stop("give the number of the test to show its code (see listTests())")
    FUN <- .list_of_tests[which]
    tl <- names(FUN)
    cat("Test no.", which, ": ", tl, "\n\n", FUN[[1]], " <- ", sep = "")
    get(FUN[[1]])
}
