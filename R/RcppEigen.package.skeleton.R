## RcppEigen.package.skeleton.R: makes a skeleton for a package that wants to use RcppEigen
##
## Copyright (C) 2011 - 2019  Douglas Bates, Dirk Eddelbuettel and Romain Francois
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

RcppEigen.package.skeleton <- function(name= "anRpackage", list = character(),
                                       environment = .GlobalEnv,
                                       path = ".", force = FALSE,
                                       code_files = character(),
                                       example_code = TRUE) {

    env <- parent.frame(1)

    if (!length(list)) {
        fake <- TRUE
        assign("Rcpp.fake.fun", function() {}, envir = env)
    } else {
        fake <- FALSE
    }

    haveKitten <- requireNamespace("pkgKitten", quietly=TRUE)
    skelFunUsed <- ifelse(haveKitten, pkgKitten::kitten, package.skeleton)
    skelFunName <- ifelse(haveKitten, "kitten", "package.skeleton")
    message("\nCalling ", skelFunName, " to create basic package.")

    ## first let the traditional version do its business
    call <- match.call()
    call[[1]] <- skelFunUsed
    if ("example_code" %in% names(call)) {
        call[["example_code"]] <- NULL    # remove the example_code argument
    }
    if (! haveKitten) {                 # in the package.skeleton() case
        if (fake) {
            call[["list"]] <- "Rcpp.fake.fun"
        }
    } else {
        if (force) {
            call[["force"]] <- NULL
        }
    }

    tryCatch(eval(call, envir = env),
             error = function(e) {
                 cat(paste(e, "\n")) # print error
                 stop(paste("error while calling `", skelFunName, "`", sep=""))
             })

    message("\nAdding RcppEigen settings")

    ## now pick things up
    root <- file.path(path, name)

    ## Add Rcpp to the DESCRIPTION
    DESCRIPTION <- file.path(root, "DESCRIPTION")
    if (file.exists(DESCRIPTION)) {
        x <- cbind(read.dcf(DESCRIPTION),
                   "Imports" = sprintf("Rcpp (>= %s)",
                   packageDescription("Rcpp")[["Version"]]),
                   "LinkingTo" = "Rcpp, RcppEigen")
        write.dcf(x, file = DESCRIPTION)
        message(" >> added Imports: Rcpp")
        message(" >> added LinkingTo: Rcpp, RcppEigen")
    }

    ## add a useDynLib to NAMESPACE,
    NAMESPACE <- file.path(root, "NAMESPACE")
    lines <- readLines(NAMESPACE)
    if (! grepl("useDynLib", lines)) {
        lines <- c(sprintf("useDynLib(%s)", name),
                   "import(RcppEigen)",
                   "importFrom(Rcpp, evalCpp)",        ## ensures Rcpp instantiation
                   lines)
        writeLines(lines, con = NAMESPACE)
        message(" >> added useDynLib and importFrom directives to NAMESPACE")
    }

    ## lay things out in the src directory
    src <- file.path(root, "src")
    if (!file.exists(src)) {
        dir.create(src)
    }
    man <- file.path(root, "man")
    if (!file.exists(man)) {
        dir.create(man)
    }
    skeleton <- system.file("skeleton", package = "RcppEigen")
    Makevars <- file.path(src, "Makevars")
    if (!file.exists(Makevars)) {
        file.copy(file.path(skeleton, "Makevars"), Makevars)
        message(" >> added Makevars file with RcppEigen settings")
    }

    Makevars.win <- file.path(src, "Makevars.win")
    if (!file.exists(Makevars.win)) {
        file.copy(file.path(skeleton, "Makevars.win"), Makevars.win)
        message(" >> added Makevars.win file with RcppEigen settings")
    }

    if (example_code) {
        file.copy(file.path(skeleton, "rcppeigen_hello_world.cpp"), src)
        message(" >> added example src file using Eigen classes")
        file.copy(file.path(skeleton, "rcppeigen_hello_world.Rd"), man)
        message(" >> added example Rd file for using Eigen classes")

	Rcpp::compileAttributes(root)
        message(" >> invoked Rcpp::compileAttributes to create wrappers")
    }

    if (fake) {
        rm("Rcpp.fake.fun", envir = env)
        unlink(file.path(root, "R"  , "Rcpp.fake.fun.R"))
        unlink(file.path(root, "man", "Rcpp.fake.fun.Rd"))
    }

    invisible(NULL)
}
