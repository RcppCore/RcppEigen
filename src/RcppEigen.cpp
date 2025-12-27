// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// RcppEigen.cpp: Rcpp/Eigen glue
//
// Copyright (C) 2011 - 2023  Douglas Bates, Dirk Eddelbuettel and Romain Francois
//
// This file is part of RcppEigen.
//
// RcppEigen is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// RcppEigen is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RcppEigen.  If not, see <http://www.gnu.org/licenses/>.

#include <RcppEigen.h>

// [[Rcpp::export]]
Rcpp::IntegerVector eigen_version(bool single) {
    using Rcpp::_;
    using Rcpp::IntegerVector;

    if (single) {
        return Rcpp::wrap( 10000 * EIGEN_WORLD_VERSION +
                           100 * EIGEN_MAJOR_VERSION +
                           EIGEN_MINOR_VERSION ) ;
    }

    return IntegerVector::create(_["major"] = EIGEN_WORLD_VERSION,
                                 _["minor"] = EIGEN_MAJOR_VERSION,
                                 _["patch"] = EIGEN_MINOR_VERSION);
}

// [[Rcpp::export]]
Rcpp::List eigen_version_typed() {
    // create a vector of major, minor, patch
    auto v = Rcpp::IntegerVector::create(EIGEN_WORLD_VERSION, EIGEN_MAJOR_VERSION, EIGEN_MINOR_VERSION);
    // and place it in a list (as e.g. packageVersion() in R returns)
    auto l = Rcpp::List::create(v);
    // and class it as 'package_version' accessing print() etc methods
    l.attr("class") = Rcpp::CharacterVector::create("package_version", "numeric_version");
    return l;
}

// [[Rcpp::export]]
bool Eigen_SSE() {
    return Rcpp::wrap(Eigen::SimdInstructionSetsInUse());
}

// [[Rcpp::export]]
int EigenNbThreads() {
    return Eigen::nbThreads();
}

// [[Rcpp::export]]
void EigenSetNbThreads(int n) {
    Eigen::setNbThreads(n);
}
