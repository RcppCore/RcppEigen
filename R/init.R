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
    .pkgenv[["nb_threads"]] <- EigenNbThreads()		# #nocov
}

##' Throttle (or Reset) (Rcpp)Eigen to Two Cores
##'
##' Helper functions to throttle use of cores by RcppEigen-internal code.
##' On package load, the initial value is saved and used to reset the value.
##' @param n Integer value of desired cores, default is two
RcppEigen_throttle_cores <- function(n = 2) {
    EigenSetNbThreads(n)
}

##' @rdname RcppEigen_throttle_cores
eigen_reset_cores <- function() {
    n <- .pkgenv[["nb_threads"]]
    EigenSetNbThreads(n)
}
