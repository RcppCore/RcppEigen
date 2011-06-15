#!/usr/bin/r -t
#
# Copyright (C) 2011	Dirk Eddelbuettel, Romain Francois and Douglas Bates
#
# This file is part of RcppEigen
#
# RcppEigen is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# RcppEigen is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with RcppEigen.  If not, see <http://www.gnu.org/licenses/>.

.setUp <- function(){
    suppressMessages(require(datasets))
    suppressMessages(require(RcppEigen))
}

test.fastLm <- function() {
    data(trees)
    flm <- .Call("fastLm",
                 cbind(1, log(trees$Girth)),
                 log(trees$Volume),
                 package="RcppEigen")
    fit <- lm(log(Volume) ~ log(Girth), data=trees)

    checkEquals(as.numeric(flm$coefficients), as.numeric(coef(fit)),
                msg="fastLm.coef")
    checkEquals(as.numeric(flm$stderr), as.numeric(coef(summary(fit))[,2]),
                msg="fastLm.stderr")
}

test.fastLm.Bench <- function() {
    data(trees)
    flm <- .Call("fastLmBench",
                 cbind(1, log(trees$Girth)),
                 log(trees$Volume),
                 package="RcppEigen")
    fit <- lm(log(Volume) ~ log(Girth), data=trees)

    checkEquals(as.numeric(flm$coefficients), as.numeric(coef(fit)),
                msg="fastLm.coef")
    checkEquals(as.numeric(flm$se), as.numeric(coef(summary(fit))[,2]),
                msg="fastLm.stderr")
}

test.fastLm.Chol1 <- function() {
    data(trees)
    flm <- .Call("fastLmChol1",
                 cbind(1, log(trees$Girth)),
                 log(trees$Volume),
                 package="RcppEigen")
    fit <- lm(log(Volume) ~ log(Girth), data=trees)

    checkEquals(as.numeric(flm$coefficients), as.numeric(coef(fit)),
                msg="fastLm.coef")
    checkEquals(as.numeric(flm$se), as.numeric(coef(summary(fit))[,2]),
                msg="fastLm.stderr")
}

test.fastLm.Chol2 <- function() {
    data(trees)
    flm <- .Call("fastLmChol2",
                 cbind(1, log(trees$Girth)),
                 log(trees$Volume),
                 package="RcppEigen")
    fit <- lm(log(Volume) ~ log(Girth), data=trees)

    checkEquals(as.numeric(flm$coefficients), as.numeric(coef(fit)),
                msg="fastLm.coef")
    checkEquals(as.numeric(flm$se), as.numeric(coef(summary(fit))[,2]),
                msg="fastLm.stderr")
}

test.fastLm.formula <- function() {
    data(trees)
    flm <- fastLm(log(Volume) ~ log(Girth), data=trees)
    fit <- lm(log(Volume) ~ log(Girth), data=trees)

    checkEquals(flm$coefficients, coef(fit), msg="fastLm.formula.coef")
    checkEquals(as.numeric(flm$stderr), as.numeric(coef(summary(fit))[,2]),
                msg="fastLm.formula.stderr")
}

