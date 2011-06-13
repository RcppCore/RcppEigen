// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// fastLm.cpp: Rcpp/Eigen example of a simple lm() alternative
//
// Copyright (C)  2010 Dirk Eddelbuettel, Romain Francois and Douglas Bates
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
// along with RcppArmadillo.  If not, see <http://www.gnu.org/licenses/>.

#include <RcppEigen.h>

extern "C" SEXP fastLm(SEXP ys, SEXP Xs) {
    using namespace Eigen;
    using namespace Rcpp;
    try {
	NumericMatrix X(Xs);
	NumericVector y(ys);
	size_t n = X.nrow(), p = X.ncol();
	if ((size_t)y.size() != n)
	    throw std::invalid_argument("size mismatch");
	
	MatrixXd A = Map<MatrixXd>(X.begin(), n, p); // shares storage
	VectorXd b = Map<VectorXd>(y.begin(), n);
    
	ColPivHouseholderQR<MatrixXd> mqr(A);      // decompose the model matrix
	size_t        r = mqr.rank();
	ptrdiff_t    df = n - r;
	VectorXd   coef(p);
	VectorXd fitted(n);
	VectorXd     se(p);
	MatrixXd   Rinv;
	if (r < p) {		// Handle the rank-deficient case
	    MatrixXd Atrunc = (A * mqr.colsPermutation()).leftCols(r);
	    HouseholderQR<MatrixXd> QR(Atrunc);
	    VectorXd coefTrunc = QR.solve(b);
	    fitted = Atrunc * coefTrunc;
	    MatrixXd   R = QR.matrixQR().topLeftCorner(r, r);
	    Rinv = TriangularView<MatrixXd, Upper>(R).solve(MatrixXd::Identity(r, r));
	    coef.topRows(r) = coefTrunc;
	    std::fill(coef.data() + r, coef.data() + p,
		      std::numeric_limits<double>::quiet_NaN());
	    coef = mqr.colsPermutation() * coef;
	    std::fill(se.data() + r, se.data() + p,
		      std::numeric_limits<double>::quiet_NaN());
	} else {
	    coef = mqr.solve(b); 
	    fitted = A * coef;
	    MatrixXd  R = mqr.matrixQR().topRows(p);
	    Rinv = TriangularView<MatrixXd, Upper>(R).solve(MatrixXd::Identity(p, p));
	}
	VectorXd      res = b - fitted;
	double          s = std::sqrt(res.squaredNorm()/df);
	se.topRows(r) = Rinv.rowwise().norm() * s;
	se = mqr.colsPermutation() * se;

	return List::create(_["coefficients"]  = coef,
			    _["rank"]          = r,
			    _["df"]            = df,
			    _["perm"]          = mqr.colsPermutation().indices(),
			    _["stderr"]        = se,
			    _["s"]             = s,
			    _["residuals"]     = res,
			    _["fitted.values"] = fitted,
			    _["Rinv"]          = Rinv
	    );
    } catch( std::exception &ex ) {
	forward_exception_to_r( ex );
    } catch(...) { 
	::Rf_error( "c++ exception (unknown reason)" ); 
    }
    return R_NilValue; // -Wall
}





