## init.R: Startup
##
## Copyright (C)  2025  Dirk Eddelbuettel
##
## This file is part of RcppEigen.
##
## RcppEigen is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## RcppEigen is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with RcppEigen.  If not, see <http://www.gnu.org/licenses/>.

.pkgenv <- new.env(parent=emptyenv())

.onLoad <- function(libname, pkgname) {
    ## simple fallback: 'Ncpus' (if set) or else all cpus seen by OpenMP
    ncores <- getOption("Ncpus", EigenNbThreads())
    ## consider OMP_THREAD_LIMIT (cf Writing R Extensions), gets NA if envvar unset
    ompcores <- as.integer(Sys.getenv("OMP_THREAD_LIMIT"))
    ## keep the smaller value, omitting NA
    ncores <- min(na.omit(c(ncores, ompcores)))
    .pkgenv[["nb_threads"]] <- ncores		# #nocov
    RcppEigen_throttle_cores(ncores)
}

.onAttach <- function(libname, pkgname) {
    if (interactive()) {
        packageStartupMessage("RcppEigen ", packageVersion("RcppEigen"),
                              " using ", .pkgenv[["nb_threads"]], " cores. See ",
                              "'help(\"RcppEigen-package\")' for details.")
    }
}

##' Throttle (or Reset) (Rcpp)Eigen Core Usage
##'
##' Helper functions to throttle use of cores by RcppEigen-internal code.
##' On package load, the initial value is saved and used to reset the value.
##' @param n Integer value of desired cores, default is the value set at package
##' startup reflecting the smallest value among the total number of available
##' cores (or one if compiled without OpenMP support), the value of option
##' \code{Ncpus} and the value of environment variable \code{OMP_THREAD_LIMIT}.
##' @return Only \code{EigenNbThreads()} returns a value, the current value of
##' the number of cores used. The other functions are invoked for their side
##' effect of affecting the count of cores used.
##' @seealso \code{\link{RcppEigen-package}}
RcppEigen_throttle_cores <- function(n) {
    if (missing(n)) n <- .pkgenv[["nb_threads"]]
    EigenSetNbThreads(n)
}

##' @rdname RcppEigen_throttle_cores
RcppEigen_reset_cores <- function() {
    EigenSetNbThreads(.pkgenv[["nb_threads"]])
}
