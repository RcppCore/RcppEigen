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

test.wrap.R <- function(){

    fx <- cxxfunction( , '

    // The eval() method is necessary because the static methods Zero()
    // and Identity() return expression objects, so that more general
    // expressions involving them can be evaluated more effectively.  We
    // have not yet written a wrap method for expression objects.
    List vecs = List::create(
        _["Vec<complex>"] = Eigen::VectorXcd::Zero(5).eval(),
        _["Vec<double>"]  = Eigen::VectorXd::Zero(5).eval(),
        _["Vec<float>"]   = Eigen::VectorXf::Zero(5).eval(),
        _["Vec<int>"]     = Eigen::VectorXi::Zero(5).eval()
    );

    // A VectorX<T> behaves as a matrix with one column but is converted to
    // a vector object in R, not a matrix of one column.  The distinction is
    // that VectorX<T> objects are defined at compile time to have one column,
    // whereas a MatrixX<T> has a dynamic number of columns that is set to 1
    // during execution of the code.  A MatrixX<T> object can be resized to have
    // a different number of columns.  A VectorX<T> object cannot.
    List cols = List::create(
        _["Col<complex>"] = Eigen::MatrixXcd::Zero(5, 1).eval(),
        _["Col<double>"]  = Eigen::MatrixXd::Zero(5, 1).eval(),
        _["Col<float>"]   = Eigen::MatrixXf::Zero(5, 1).eval(),
        _["Col<int>"]     = Eigen::MatrixXi::Zero(5, 1).eval()
    );

    List rows = List::create(
        _["Row<complex>"] = Eigen::RowVectorXcd::Zero(5).eval(),
        _["Row<double>"]  = Eigen::RowVectorXd::Zero(5).eval(),
        _["Row<float>"]   = Eigen::RowVectorXf::Zero(5).eval(),
        _["Row<int>"]     = Eigen::RowVectorXi::Zero(5).eval()
    );

    List matrices = List::create(
        _["Mat<complex>"] = Eigen::MatrixXcd::Identity(3, 3).eval(),
        _["Mat<double>"]  = Eigen::MatrixXd::Identity(3, 3).eval(),
        _["Mat<float>"]   = Eigen::MatrixXf::Identity(3, 3).eval(),
        _["Mat<int>"]     = Eigen::MatrixXi::Identity(3, 3).eval()
    );

    // ArrayXX<t> objects have the same structure as matrices but allow
    // componentwise arithmetic.  A * B is matrix multiplication for
    // matrices and componentwise multiplication for arrays.
    List arrays2 = List::create(
        _["Arr2<complex>"] = Eigen::ArrayXXcd::Zero(3, 3).eval(),
        _["Arr2<double>"]  = Eigen::ArrayXXd::Zero(3, 3).eval(),
        _["Arr2<float>"]   = Eigen::ArrayXXf::Zero(3, 3).eval(),
        _["Arr2<int>"]     = Eigen::ArrayXXi::Zero(3, 3).eval()
    );

    // ArrayX<t> objects have the same structure as VectorX<T> objects
    // but allow componentwise arithmetic, including functions like exp, log,
    // sqrt, ...
    List arrays1 = List::create(
        _["Arr1<complex>"] = Eigen::ArrayXcd::Zero(5).eval(),
        _["Arr1<double>"]  = Eigen::ArrayXd::Zero(5).eval(),
        _["Arr1<float>"]   = Eigen::ArrayXf::Zero(5).eval(),
        _["Arr1<int>"]     = Eigen::ArrayXi::Zero(5).eval()
    );

    List operations = List::create(
        _["Op_seq"]  = Eigen::ArrayXd::LinSpaced(6, 1, 10).eval(),  // arguments are length.out, start, end
        _["Op_log"]  = Eigen::ArrayXd::LinSpaced(6, 1, 10).log().eval(),
        _["Op_exp"]  = Eigen::ArrayXd::LinSpaced(6, 1, 10).exp().eval(),
        _["Op_sqrt"] = Eigen::ArrayXd::LinSpaced(6, 1, 10).sqrt().eval(),
        _["Op_cos"]  = Eigen::ArrayXd::LinSpaced(6, 1, 10).cos().eval()
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

    return output ;
	' , plugin = "RcppEigen" )

    res <- fx()

    checkEquals( res[[1]][[1]], complex(5), msg = "VectorXcd::Zero(5)")
    checkEquals( res[[1]][[2]], double(5), msg = "VectorXd::Zero(5)")
    checkEquals( res[[1]][[3]], double(5), msg = "VectorXf::Zero(5)")
    checkEquals( res[[1]][[4]], integer(5), msg = "VectorXi::Zero(5)")
    
    checkEquals( res[[2]][[1]], (1+0i) * diag(nr=3L), msg = "MatrixXcd::Identity(3,3)")
    checkEquals( res[[2]][[2]], diag(nr=3L), msg = "MatrixXd::Identity(3,3)")
    checkEquals( res[[2]][[3]], diag(nr=3L), msg = "MatrixXf::Identity(3,3)")
    checkEquals( res[[2]][[4]], matrix(as.integer((diag(nr=3L))),nr=3L), msg = "MatrixXi::Identity(3,3)")

    checkEquals( res[[3]][[1]], matrix(complex(5), nr=1L), msg = "RowVectorXcd::Zero(5)" )
    checkEquals( res[[3]][[2]], matrix(numeric(5), nr=1L), msg = "RowVectorXd::Zero(5)" )
    checkEquals( res[[3]][[3]], matrix(numeric(5), nr=1L), msg = "RowVectorXf::Zero(5)" )
    checkEquals( res[[3]][[4]], matrix(integer(5), nr=1L), msg = "RowVectorXi::Zero(5)" )

    checkEquals( res[[4]][[1]], as.matrix(complex(5)), msg = "MatrixXcd::Zero(5, 1)")
    checkEquals( res[[4]][[2]], as.matrix(numeric(5)), msg = "MatrixXd::Zero(5, 1)")
    checkEquals( res[[4]][[3]], as.matrix(numeric(5)), msg = "MatrixXf::Zero(5, 1)")
    checkEquals( res[[4]][[4]], as.matrix(integer(5)), msg = "MatrixXi::Zero(5, 1)")

    checkEquals( res[[5]][[1]], matrix(complex(9L), nc=3L), msg = "ArrayXXcd::Zero(3,3)")
    checkEquals( res[[5]][[2]], matrix(numeric(9L), nc=3L), msg = "ArrayXXd::Zero(3,3)")
    checkEquals( res[[5]][[3]], matrix(numeric(9L), nc=3L), msg = "ArrayXXf::Zero(3,3)")
    checkEquals( res[[5]][[4]], matrix(integer(9L), nc=3L), msg = "ArrayXXi::Zero(3,3)")

    checkEquals( res[[6]][[1]], complex(5), msg = "ArrayXcd::Zero(5)")
    checkEquals( res[[6]][[2]], double(5), msg = "ArrayXd::Zero(5)")
    checkEquals( res[[6]][[3]], double(5), msg = "ArrayXf::Zero(5)")
    checkEquals( res[[6]][[4]], integer(5), msg = "ArrayXi::Zero(5)")

    oneTen <- seq(1, 10, length.out=6L)
    
    checkEquals( res[[7]][[1]], oneTen,       msg = "Op_seq")
    checkEquals( res[[7]][[2]], log(oneTen),  msg = "Op_log")
    checkEquals( res[[7]][[3]], exp(oneTen),  msg = "Op_exp")
    checkEquals( res[[7]][[4]], sqrt(oneTen), msg = "Op_sqrt")
    checkEquals( res[[7]][[5]], cos(oneTen),  msg = "Op_cos")
    
}

test.as.Col <- function(){
    fx <- cxxfunction( signature(input_ = "list" ) , '

    List input(input_) ;
    Eigen::VectorXi                                m1 = input[0] ; /* implicit as */
    Eigen::VectorXd                                m2 = input[1] ; /* implicit as */
    Eigen::Matrix<unsigned int, Eigen::Dynamic, 1> m3 = input[0] ; /* implicit as */
    Eigen::VectorXf                                m4 = input[1] ; /* implicit as */

    List res = List::create(m1.sum(), m2.sum(), m3.sum(), m4.sum());

    return res ;

    ', plugin = "RcppEigen" )

    res <- fx( list( 1:10, as.numeric(1:10) ) )
    checkEquals( unlist( res ), rep(55.0, 4 ), msg = "as<Col>" )
}

if (FALSE) {

test.as.Mat <- function(){

    fx <- cxxfunction( signature(input_ = "list" ) , '
    List input(input_) ;
    Eigen::MatrixXi                                             m1 = input[0] ; /* implicit as */
    Eigen::MatrixXd                                             m2 = input[1] ; /* implicit as */
    Eigen::Matrix<unsigned int, Eigen::Dynamic, Eigen::Dynamic> m3 = input[0] ; /* implicit as */
    Eigen::MatrixXf                                             m4 = input[1] ; /* implicit as */

    List res = List::create(m1.sum(), m2.sum(), m3.sum(), m4.sum());

    return res ;
    ', plugin = "RcppEigen" )

    integer_mat <- matrix( as.integer(diag(4)), ncol = 4, nrow = 4 )
    numeric_mat <- diag(5)
    res <- fx( list( integer_mat, numeric_mat ) )
    checkEquals( unlist( res), c(4L, 5L, 4L, 5L ), msg = "as<Mat>" )
}

test.wrap.Glue <- function(){

    fx <- cxxfunction( , '

    arma::mat m1 = arma::eye<arma::mat>( 3, 3 ) ;
    arma::mat m2 = arma::eye<arma::mat>( 3, 3 ) ;

    List res ;
    res["mat+mat"] = m1 + m2 ;
    return res ;

    ', plugin = "RcppArmadillo" )

	res <- fx()
    checkEquals( res[[1]], 2*diag(3), msg = "wrap(Glue)" )
}

test.wrap.Op <- function(){

    fx <- cxxfunction( , '

    arma::mat m1 = arma::eye<arma::mat>( 3, 3 ) ;

    List res ;
    res["- mat"] = - m1 ;
    return res ;

    ', plugin = "RcppArmadillo" )
    res <- fx()
    checkEquals( res[[1]], -1*diag(3), msg = "wrap(Op)" )
}

test.as.Col <- function(){
    fx <- cxxfunction( signature(input_ = "list" ) , '

    List input(input_) ;
    arma::icolvec m1 = input[0] ; /* implicit as */
    arma::colvec  m2 = input[1] ; /* implicit as */
    arma::ucolvec m3 = input[0] ; /* implicit as */
    arma::fcolvec m4 = input[1] ; /* implicit as */

    List res = List::create(
    	arma::accu( m1 ),
    	arma::accu( m2 ),
    	arma::accu( m3 ),
    	arma::accu( m4 ) ) ;

    return res ;

    ', plugin = "RcppArmadillo" )

    res <- fx( list( 1:10, as.numeric(1:10) ) )
    checkEquals( unlist( res ), rep(55.0, 4 ), msg = "as<Col>" )
}

test.as.Row <- function(){
    fx <- cxxfunction( signature(input_ = "list" ) , '

    List input(input_) ;
    arma::irowvec m1 = input[0] ; /* implicit as */
    arma::rowvec  m2 = input[1] ; /* implicit as */
    arma::urowvec m3 = input[0] ; /* implicit as */
    arma::frowvec m4 = input[1] ; /* implicit as */

    List res = List::create(
    	arma::accu( m1 ),
    	arma::accu( m2 ),
    	arma::accu( m3 ),
    	arma::accu( m4 ) ) ;
    return res ;

	', plugin = "RcppArmadillo" )

    res <- fx( list( 1:10, as.numeric(1:10) ) )
    checkEquals( unlist( res ), rep(55.0, 4 ), msg = "as<Row>" )
}

test.cxmat <- function(){

    fx <- cxxfunction( signature() , '

    arma::cx_mat m1  = arma::eye<arma::cx_mat> ( 3, 3 ) ;
    arma::cx_fmat m2 = arma::eye<arma::cx_fmat>( 3, 3 ) ;
    return List::create( _["double"] = m1, _["float"] = m2 ) ;

    ', plugin = "RcppArmadillo" )
    checkEquals( fx(),
		list( double = (1+0i)*diag(3), float = (1+0i)*diag(3) ),
		msg = "support for complex matrices" )

}

test.mtOp <- function(){

    fx <- cxxfunction( signature() , '

    std::complex<double> x( 1.0, 2.0 ) ;
    arma::mat m1  = arma::eye<arma::mat> ( 3, 3 ) ;

    return wrap( x * m1 ) ;

    ', plugin = "RcppArmadillo" )
    checkEquals( fx(),
		(1+2i)*diag(3),
		msg = "support for mtOp" )

}

test.mtGlue <- function(){

    fx <- cxxfunction( signature() , '

    arma::imat m2 = arma::eye<arma::imat> ( 3, 3 ) ;
    arma::mat m1  = arma::eye<arma::mat> ( 3, 3 ) ;

    return wrap( m1 + m2 ) ;

    ', plugin = "RcppArmadillo" )
    checkEquals( fx(),
		2.0 * diag(3) ,
		msg = "support for mtOp" )

}


test.sugar <- function(){

    fx <- cxxfunction( signature(x= "numeric") , '
    NumericVector xx(x) ;
    arma::mat m = forward( xx + xx ) ;
    return wrap( m ) ;

    ', plugin = "RcppArmadillo" )
    checkEquals( fx(1:10),
		matrix( 2*(1:10), nrow = 10 ) ,
		msg = "RcppArmadillo and sugar" )

}

test.sugar.cplx <- function(){

    fx <- cxxfunction( signature(x= "complex") , '
    ComplexVector xx(x) ;
    arma::cx_mat m = forward( exp( xx ) ) ;

    return wrap( m ) ;

    ', plugin = "RcppArmadillo" )
    x <- 1:10*(1+1i)
    checkEquals( fx(x),
		matrix( exp(x), nrow = 10 ) ,
		msg = "RcppArmadillo and sugar (complex)" )

}

test.armadillo.sugar.ctor <- function(){

    fx <- cxxfunction( signature(x= "numeric") , '
    NumericVector xx(x) ;
    arma::mat m = xx + xx ;
    arma::colvec co = xx ;
    arma::rowvec ro = xx ;
    return List::create(
    	_["mat"] = m + m,
    	_["rowvec"] = ro,
    	_["colvec"] = co
    );
    ', plugin = "RcppArmadillo" )
    checkEquals( fx(1:10),
		list(
                     mat = matrix( 4*(1:10), nrow = 10 ),
                     rowvec = matrix( 1:10, nrow = 1 ),
                     colvec = matrix( 1:10, ncol = 1 )
                     )
		,
		msg = "Mat( sugar expression )" )

}


test.armadillo.sugar.matrix.ctor <- function(){

    inc <- '
    double norm( double x, double y){
		return ::sqrt( x*x + y*y );
    }
    '
    fx <- cxxfunction( signature(x= "numeric") , '
    NumericVector xx(x) ;
    NumericVector yy = NumericVector::create( 1 ) ;
    arma::mat m = diag( xx ) ;
    arma::colvec co = outer( xx, yy, ::norm ) ;
    arma::rowvec ro = outer( yy, xx, ::norm ) ;
    return List::create(
    	_["mat"] = m + m ,
    	_["rowvec"] = ro,
    	_["colvec"] = co
    );
    ', plugin = "RcppArmadillo", includes = inc )
    res <- fx(1:10)
    norm <- function(x, y) sqrt( x*x + y*y )
    checkEquals( res,
		list(
                     mat = diag(2*(1:10)),
                     rowvec = outer( 1, 1:10, norm ),
                     colvec = outer( 1:10, 1, norm )
                     ),
		msg = "Mat( sugar expression )" )

}

test.armadillo.rtti.check <- function() {

    inc <- '
    void blah(arma::mat& X) {
         X.set_size(5,5);
    }
    '
    src <- '
    arma::vec V;
    blah(V); // if blah() worked, we have a problem
    '
    fun <- cxxfunction(signature(), body=src, inc=inc, plugin = "RcppArmadillo")

    checkException(fun(), msg="RTTI check on matrix constructor exception")

}
}
