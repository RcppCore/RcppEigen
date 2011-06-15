## fastLm.R: Rcpp/Eigen implementation of lm()
##
## Copyright (C)  2011 Douglas Bates, Dirk Eddelbuettel and Romain Francois
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

fastLmPure <- function(X, y) {

    stopifnot(is.matrix(X), is.numeric(y), NROW(y)==nrow(X))

    .Call("fastLm", X, y, package="RcppEigen")
}

fastLm <- function(X, ...) UseMethod("fastLm")

fastLm.default <- function(X, y, ...) {

    X <- as.matrix(X)
    y <- as.numeric(y)

    res <- fastLmPure(X, y)
    names(res$coefficients) <- colnames(X)
    res$call <- match.call()

    class(res) <- "fastLm"
    res
}

print.fastLm <- function(x, ...) {
    cat("\nCall:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    print(x$coefficients, digits=5)
}

summary.fastLm <- function(object, ...) {
    se <- object$stderr
    tval <- coef(object)/se

    object$coefficients <- cbind(Estimate = object$coefficients,
                                 StdErr   = se,
                                 t.value  = tval,
                                 p.value  = 2*pt(-abs(tval), df=object$df))

    ## cf src/stats/R/lm.R and case with no weights and an intercept
    f <- object$fitted.values
    r <- object$residuals
    mss <- sum((f - mean(f))^2)
    rss <- sum(r^2)

    object$r.squared <- mss/(mss + rss)
    df.int <- 1 		# case of intercept
    n <- length(f)
    rdf <- object$df
    object$adj.r.squared <- 1 - (1 - object$r.squared) * ((n - df.int)/rdf)
    class(object) <- "summary.fastLm"
    object
}

print.summary.fastLm <- function(x, ...) {
    cat("\nCall:\n")
    print(x$call)
    cat("\n")

    printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE)
    digits <- max(3, getOption("digits") - 3)
    cat("\nResidual standard error: ", formatC(x$s, digits=digits), " on ",
        formatC(x$df), " degrees of freedom\n", sep="")
    cat("Multiple R-squared: ", formatC(x$r.squared, digits=digits),
        ",\tAdjusted R-squared: ",formatC(x$adj.r.squared, digits=digits),
        "\n", sep="")
    invisible(x)
}

fastLm.formula <- function(formula, data=list(), ...) {
    mf <- model.frame(formula=formula, data=data)
    X <- model.matrix(attr(mf, "terms"), data=mf)
    y <- model.response(mf)

    res <- fastLm.default(X, y, ...)
    res$call <- match.call()
    # I think this is redundant.  The formula is available as res$call$formula
    res$formula <- formula
    res
}

predict.fastLm <- function(object, newdata=NULL, ...) {
    if (is.null(newdata)) {
        y <- fitted(object)
    } else {
        if (!is.null(object$formula)) {
            x <- model.matrix(object$formula, newdata)
        } else {
            x <- newdata
        }
        y <- as.vector(x %*% coef(object))
    }
    y
}
