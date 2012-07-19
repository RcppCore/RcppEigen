## fastLm.R: Rcpp/Eigen implementation of lm()
##
## Copyright (C)  2011 - 2012  Douglas Bates, Dirk Eddelbuettel and Romain Francois
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

fastLmPure <- function(X, y, method = 0L) {

    stopifnot(is.matrix(X), is.numeric(y), NROW(y)==nrow(X))

    .Call("fastLm", X, y, as.integer(method[1]), PACKAGE="RcppEigen")
}

fastLm <- function(X, ...) UseMethod("fastLm")

fastLm.default <- function(X, y, method = 0L, ...) {

    X <- as.matrix(X)
    y <- as.numeric(y)

    res <- fastLmPure(X, y, as.integer(method[1]))
    res$call <- match.call()
    res$intercept <- any(apply(X, 2, function(x) all(x == x[1])))

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
    coef <- object$coefficients
    se   <- object$se
    tval <- coef/se

    object$coefficients <- cbind(Estimate     = coef,
                                 "Std. Error" = se,
                                 "t value"    = tval,
                                 "Pr(>|t|)"   = 2*pt(-abs(tval), df=object$df))

    ## cf src/stats/R/lm.R and case with no weights and an intercept
    f <- object$fitted.values
    r <- object$residuals
    #mss <- sum((f - mean(f))^2)
    mss <- if (object$intercept) sum((f - mean(f))^2) else sum(f^2)
    rss <- sum(r^2)

    object$r.squared <- mss/(mss + rss)
    df.int <- if (object$intercept) 1L else 0L
    n <- length(f)
    rdf <- object$df
    object$adj.r.squared <- 1 - (1 - object$r.squared) * ((n - df.int)/rdf)
    class(object) <- "summary.fastLm"
    object
}

print.summary.fastLm <- function(x, ...) {
    cat("\nCall:\n")
    print(x$call)
    cat("\nResiduals:\n")
    digits <- max(3, getOption("digits") - 3)
    print(summary(x$residuals, digits=digits)[-4])
    cat("\n")

    printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE, ...)
    cat("\nResidual standard error: ", formatC(x$s, digits=digits), " on ",
        formatC(x$df), " degrees of freedom\n", sep="")
    cat("Multiple R-squared: ", formatC(x$r.squared, digits=digits),
        ",\tAdjusted R-squared: ",formatC(x$adj.r.squared, digits=digits),
        "\n", sep="")
    invisible(x)
}

fastLm.formula <- function(formula, data=list(), method = 0L, ...) {
    mf <- model.frame(formula=formula, data=data)
    X <- model.matrix(attr(mf, "terms"), data=mf)
    y <- model.response(mf)

    res <- fastLm.default(X, y, method=method, ...)
    res$call <- match.call()
    ## I think this is redundant.  The formula is available as res$call$formula
    res$formula <- formula
    res$intercept <- attr(attr(mf, "terms"), "intercept")
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
