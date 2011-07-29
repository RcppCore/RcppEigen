#!/usr/bin/r -t
#
# Copyright (C)      2011 Douglas Bates, Dirk Eddelbuettel and Romain Francois
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
    suppressMessages(require(inline))
}

test.wrapSparse.R <- function(){

    fx <- cxxfunction( , '

    Eigen::SparseMatrix<double>  mm(9,3);
    mm.reserve(9);
    for (int j = 0; j < 3; ++j) {
        mm.startVec(j);
        for (int i = 3 * j; i < 3 * (j + 1); ++i)
            mm.insertBack(i, j) = 1.;
    }
    mm.finalize();
    return wrap(mm);
' , plugin = "RcppEigen" )

    res <- fx()
    rr <- Matrix::t(as(gl(3,3), "sparseMatrix"))
    colnames(rr) <- NULL
    checkEquals( res, rr, msg = "Sparsematrix wrap")
}

test.solveCholmod.R <- function() {
    suppressMessages(require("Matrix", character.only=TRUE))
    data("KNex", package = "Matrix")

    fx <- cxxfunction( signature(input_ = "list"), '
    using Eigen::VectorXd;
    using Eigen::MatrixXd;
    using Eigen::Lower;
    using Eigen::Map;
    using Eigen::MappedSparseMatrix;
    using Eigen::SparseMatrix;
    using Eigen::CholmodDecomposition;
    using Eigen::CholmodAuto;
    using Eigen::Success;

    List input(input_);
    const MappedSparseMatrix<double> m1 = input[0];
    const Map<VectorXd>         v1 = input[1];
    SparseMatrix<double>        m2(m1.cols(), m1.cols());
    m2.selfadjointView<Lower>().rankUpdate(m1.adjoint());

    CholmodDecomposition<SparseMatrix<double> > ff(m2);
    VectorXd                   res = ff.solve(m1.adjoint() * v1);
    
    return List::create(_["res"]   = res,
                        _["rows"]  = ff.rows(),
                        _["cols"]  = ff.cols());
',
                      plugin = "RcppEigen")

    rr <- fx(KNex)
    checkEquals(rr[[1]], as.vector(solve(crossprod(KNex[[1]]),
                                         crossprod(KNex[[1]], KNex[[2]]))),
                                       "Cholmod solution")
}

test.solveCholmodRect.R <- function() {
    suppressMessages(require("Matrix", character.only=TRUE))
    data("KNex", package = "Matrix")

    fx <- cxxfunction( signature(input_ = "list"), '
    using Eigen::VectorXd;
    using Eigen::MatrixXd;
    using Eigen::Lower;
    using Eigen::Map;
    using Eigen::MappedSparseMatrix;
    using Eigen::SparseMatrix;
    using Eigen::CholmodDecomposition;
    using Eigen::CholmodAuto;
    using Eigen::Success;

    List input(input_);
    const MappedSparseMatrix<double> m1 = input[0];
    const Map<VectorXd>         v1 = input[1];
    SparseMatrix<double>        m2(m1.adjoint());

    CholmodDecomposition<SparseMatrix<double> > ff(m2);
    VectorXd                   res = ff.solve(m2 * v1);
    
    return List::create(_["res"]   = res,
                        _["rows"]  = ff.rows(),
                        _["cols"]  = ff.cols());
',
                      plugin = "RcppEigen")

    rr <- fx(KNex)
    checkEquals(rr[[1]], as.vector(solve(crossprod(KNex[[1]]),
                                         crossprod(KNex[[1]], KNex[[2]]))),
                                       "Cholmod solution")
}
