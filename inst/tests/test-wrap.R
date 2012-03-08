#
# Copyright (C) 2012	Douglas Bates, Dirk Eddelbuettel and Romain Francois
#
# This file is part of RcppEigen.
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
# along with Rcpp.  If not, see <http://www.gnu.org/licenses/>.

library("inline")
library("RcppEigen")
incl <- '
typedef Eigen::ArrayXd                   Ar1;
typedef Eigen::Map<Ar1>                 MAr1;
typedef Eigen::ArrayXXd                  Ar2;
typedef Eigen::Map<Ar2>                 MAr2;
typedef Eigen::MatrixXd                  Mat;
typedef Eigen::Map<Mat>                 MMat;
typedef Eigen::VectorXd                  Vec;
typedef Eigen::Map<Vec>                 MVec;
typedef Eigen::ArrayXi                  iAr1;
typedef Eigen::Map<iAr1>               MiAr1;
typedef Eigen::ArrayXXi                 iAr2;
typedef Eigen::Map<iAr2>               MiAr2;
typedef Eigen::MatrixXi                 iMat;
typedef Eigen::Map<iMat>               MiMat;
typedef Eigen::VectorXi                 iVec;
typedef Eigen::Map<iVec>               MiVec;
typedef Eigen::ArrayXf                  fAr1;
typedef Eigen::Map<fAr1>               MfAr1;
typedef Eigen::ArrayXXf                 fAr2;
typedef Eigen::Map<fAr2>               MfAr2;
typedef Eigen::MatrixXf                 fMat;
typedef Eigen::Map<fMat>               MfMat;
typedef Eigen::VectorXf                 fVec;
typedef Eigen::Map<fVec>               MfVec;
typedef Eigen::ArrayXcd                cdAr1;
typedef Eigen::Map<cdAr1>             McdAr1;
typedef Eigen::ArrayXXcd               cdAr2;
typedef Eigen::Map<cdAr2>             McdAr2;
typedef Eigen::MatrixXcd               cdMat;
typedef Eigen::Map<cdMat>             McdMat;
typedef Eigen::VectorXcd               cdVec;
typedef Eigen::Map<cdVec>             McdVec;
'

definitions <-
    list("wrap_vectors" = list(signature(), '
List vecs = List::create(Named("Vec<complex>", cdVec::Zero(5)),
			 Named("Vec<double>",    Vec::Zero(5)),
			 Named("Vec<float>",    fVec::Zero(5)),
			 Named("Vec<int>",      iVec::Zero(5))
    );

    // A VectorX<T> behaves as a matrix with one column but is converted to
    // a vector object in R, not a matrix of one column.  The distinction is
    // that VectorX<T> objects are defined at compile time to have one column,
    // whereas a MatrixX<T> has a dynamic number of columns that is set to 1
    // during execution of the code.  A MatrixX<T> object can be resized to have
    // a different number of columns.  A VectorX<T> object cannot.

List cols = List::create(Named("Col<complex>", cdMat::Zero(5, 1)),
			 Named("Col<double>",    Mat::Zero(5, 1)),
			 Named("Col<float>",    fMat::Zero(5, 1)),
			 Named("Col<int>",      iMat::Zero(5, 1))
    );

    List rows = List::create(
        _["Row<complex>"] = Eigen::RowVectorXcd::Zero(5),
        _["Row<double>"]  = Eigen::RowVectorXd::Zero(5),
        _["Row<float>"]   = Eigen::RowVectorXf::Zero(5),
        _["Row<int>"]     = Eigen::RowVectorXi::Zero(5)
    );

    List matrices = List::create(
        _["Mat<complex>"] = Eigen::MatrixXcd::Identity(3, 3),
        _["Mat<double>"]  = Eigen::MatrixXd::Identity(3, 3),
        _["Mat<float>"]   = Eigen::MatrixXf::Identity(3, 3),
        _["Mat<int>"]     = Eigen::MatrixXi::Identity(3, 3)
    );

    // ArrayXX<t> objects have the same structure as matrices but allow
    // componentwise arithmetic.  A * B is matrix multiplication for
    // matrices and componentwise multiplication for arrays.
    List arrays2 = List::create(
        _["Arr2<complex>"] = Eigen::ArrayXXcd::Zero(3, 3),
        _["Arr2<double>"]  = Eigen::ArrayXXd::Zero(3, 3),
        _["Arr2<float>"]   = Eigen::ArrayXXf::Zero(3, 3),
        _["Arr2<int>"]     = Eigen::ArrayXXi::Zero(3, 3)
    );

    // ArrayX<t> objects have the same structure as VectorX<T> objects
    // but allow componentwise arithmetic, including functions like exp, log,
    // sqrt, ...
    List arrays1 = List::create(
        _["Arr1<complex>"] = Eigen::ArrayXcd::Zero(5),
        _["Arr1<double>"]  = Eigen::ArrayXd::Zero(5),
        _["Arr1<float>"]   = Eigen::ArrayXf::Zero(5),
        _["Arr1<int>"]     = Eigen::ArrayXi::Zero(5)
    );

    List operations = List::create(
        _["Op_seq"]  = Eigen::ArrayXd::LinSpaced(6, 1, 10),  // arguments are length.out, start, end
        _["Op_log"]  = Eigen::ArrayXd::LinSpaced(6, 1, 10).log(),
        _["Op_exp"]  = Eigen::ArrayXd::LinSpaced(6, 1, 10).exp(),
        _["Op_sqrt"] = Eigen::ArrayXd::LinSpaced(6, 1, 10).sqrt(),
        _["Op_cos"]  = Eigen::ArrayXd::LinSpaced(6, 1, 10).cos()
    );

    List output = List::create(
    	_["vectors : VectorX<T>"]   = vecs,
    	_["matrices : MatrixX<T>"]  = matrices,
    	_["rows : RowVectorX<T>"]   = rows,
    	_["columns : MatrixX<T>"]   = cols,
        _["arrays2d : ArrayXX<T>"]  = arrays2,
        _["arrays1d : ArrayX<T>"]   = arrays1,
        _["operations : ArrayXd"]   = operations
        );
    return output;
'),
         "as_Vec" = list(signature(input_ = "list"), '
    List input(input_) ;
    Eigen::VectorXi                                m1 = input[0] ; /* implicit as */
    Eigen::VectorXd                                m2 = input[1] ; /* implicit as */
    Eigen::Matrix<unsigned int, Eigen::Dynamic, 1> m3 = input[0] ; /* implicit as */
    Eigen::VectorXf                                m4 = input[1] ; /* implicit as */

    List res = List::create(m1.sum(), m2.sum(), m3.sum(), m4.sum());

    return res ;

    ')
         )

cxxargs <- ifelse(Rcpp:::capabilities()[["initializer lists"]],"-std=c++0x","")

.setUp <- function() {
    tests <- ".rcppeigen.wrap"
    if( ! exists( tests, globalenv() )) {
        fun <- RcppEigen:::compile_unit_tests(definitions,
                                              includes=incl,
                                              cxxargs = cxxargs)
        names(fun) <- names(definitions)
        assign( tests, fun, globalenv() )
    }
}

.setUp()

res <- .rcppeigen.wrap$wrap_vectors()

context("wrap")
test_that("array", {
    expect_that(res[[1]][[1]], equals(complex(5)))
    expect_that(res[[1]][[2]], equals(double(5)))
    expect_that(res[[1]][[3]], equals(double(5)))
    expect_that(res[[1]][[4]], equals(integer(5)))

    expect_that(res[[2]][[1]], equals((1+0i) * diag(nr=3L)))
    expect_that(res[[2]][[2]], equals(diag(nr=3L)))
    expect_that(res[[2]][[3]], equals(diag(nr=3L)))
    expect_that(res[[2]][[4]], equals(matrix(as.integer((diag(nr=3L))),nr=3L)))

    expect_that(res[[3]][[1]], equals(matrix(complex(5), nr=1L)))
    expect_that(res[[3]][[2]], equals(matrix(numeric(5), nr=1L)))
    expect_that(res[[3]][[3]], equals(matrix(numeric(5), nr=1L)))
    expect_that(res[[3]][[4]], equals(matrix(integer(5), nr=1L)))

    expect_that(res[[4]][[1]], equals(as.matrix(complex(5))))
    expect_that(res[[4]][[2]], equals(as.matrix(numeric(5))))
    expect_that(res[[4]][[3]], equals(as.matrix(numeric(5))))
    expect_that(res[[4]][[4]], equals(as.matrix(integer(5))))

    expect_that(res[[5]][[1]], equals(matrix(complex(9L), nc=3L)))
    expect_that(res[[5]][[2]], equals(matrix(numeric(9L), nc=3L)))
    expect_that(res[[5]][[3]], equals(matrix(numeric(9L), nc=3L)))
    expect_that(res[[5]][[4]], equals(matrix(integer(9L), nc=3L)))

    expect_that(res[[6]][[1]], equals(complex(5)))
    expect_that(res[[6]][[2]], equals(double(5)))
    expect_that(res[[6]][[3]], equals(double(5)))
    expect_that(res[[6]][[4]], equals(integer(5)))

    oneTen <- seq(1, 10, length.out=6L)

    expect_that(res[[7]][[1]], equals(oneTen))
    expect_that(res[[7]][[2]], equals(log(oneTen)))
    expect_that(res[[7]][[3]], equals(exp(oneTen)))
    expect_that(res[[7]][[4]], equals(sqrt(oneTen)))
    expect_that(res[[7]][[5]], equals(cos(oneTen)))
})

res <- .rcppeigen.wrap$as_Vec(list(1:10, as.numeric(1:10)))

context("as_Vec")
test_that("asVec",
    expect_that(unlist(.rcppeigen.wrap$as_Vec(list(1:10, as.numeric(1:10)))),
                equals(rep.int(55, 4L)))
          )
