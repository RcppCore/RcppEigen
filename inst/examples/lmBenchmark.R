## lmBenchmark.R: Benchmark different implementations of linear model solutions
##
## Copyright (C)  2011 Douglas Bates, Dirk Eddelbuettel and Romain Francois
##
## This file is part of RcppEigen.

require("stats", character=TRUE, quietly=TRUE)
require("rbenchmark", character=TRUE, quietly=TRUE)
require("RcppEigen", character=TRUE, quietly=TRUE)

## define different versions of lm
exprs <- list()

## These versions use rank-revealing decompositions and thus can
## handle rank-deficient cases.

                                        # default version used in lm()
exprs$lm.fit <- expression(stats::lm.fit(mm, y))
                                        # versions from RcppEigen
## column-pivoted QR decomposition - similar to lm.fit
exprs$PivQR <- expression(.Call("fastLm", mm, y, 0L, PACKAGE="RcppEigen"))
## LDLt Cholesky decomposition with rank detection
exprs$LDLt <- expression(.Call("fastLm", mm, y, 3L, PACKAGE="RcppEigen"))
## SVD
exprs$SVD <- expression(.Call("fastLm", mm, y, 4L, PACKAGE="RcppEigen"))
## eigenvalues and eigenvectors of X'X
exprs$SymmEig <- expression(.Call("fastLm", mm, y, 5L, PACKAGE="RcppEigen"))

## Non-rank-revealing decompositions.  These work fine except when
## they don't.

## Unpivoted  QR decomposition
exprs$QR <- expression(.Call("fastLm", mm, y, 1L, PACKAGE="RcppEigen"))
## LLt Cholesky decomposition
exprs$LLt <- expression(.Call("fastLm", mm, y, 2L, PACKAGE="RcppEigen"))

if (require("RcppArmadillo", character=TRUE, quietly=TRUE)) {
    exprs$arma <- expression(.Call("fastLm", mm, y, PACKAGE="RcppArmadillo"))
}

if (require("RcppGSL", character=TRUE, quietly=TRUE)) {
    exprs$GSL <- expression(.Call("fastLm", mm, y, PACKAGE="RcppGSL"))
}

do_bench <- function(n=100000L, p=40L, nrep=20L) {
    mm <- cbind(1, matrix(rnorm(n * (p - 1L)), nc=p-1L))
    y <- rnorm(n)
    cat(paste("lm benchmark for n = ", n, " and p = ", p, ": nrep = ",
              nrep, "\n", sep=''))
    do.call(benchmark, c(exprs,
                         list(order="relative",
                              columns = c("test", "relative",
                              "elapsed", "user.self", "sys.self"),
                              replications = nrep)))
}

do_bench()
