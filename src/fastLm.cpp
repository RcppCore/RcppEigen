// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// fastLm.cpp: Rcpp/Eigen example of a simple lm() alternative
//
// Copyright (C)       2011 Douglas Bates, Dirk Eddelbuettel and Romain Francois
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
// in file.path(R.home("share"), "licenses").  If not, see
// <http://www.gnu.org/licenses/>.

#include "fastLm.h"

using namespace Rcpp;

lm::lm(const MMatrixXd &X, const MVectorXd &y) :
    m_n(X.rows()), m_p(X.cols()), m_coef(m_p), m_r(NA_INTEGER), m_df(m_n - m_p),
    m_perm(m_p), m_fitted(m_n), m_unsc(m_p, m_p) {}

ColPivQR::ColPivQR(const MMatrixXd &X, const MVectorXd &y) : lm(X, y) {
    PivQRType           PQR(X);	// decompose the model matrix
    PermutationType    Pmat = PQR.colsPermutation();
    m_perm                  = Pmat.indices();
    m_r                     = PQR.rank();
    m_df                    = m_n - m_r;
    MatrixXd              R = PQR.matrixQR().topLeftCorner(m_p, m_p);

    if (m_r < (int)m_p) {		// The rank-deficient case
	int           nsing = m_p - m_r;
	MatrixXd     Atrunc = (X * Pmat).leftCols(m_r);
	QRType           QR(Atrunc);
	VectorXd  coefTrunc = QR.solve(y);
	m_fitted            = Atrunc * coefTrunc;
	m_coef.topRows(m_r) = coefTrunc;
	std::fill(m_coef.data() + m_r, m_coef.data() + m_p, NA_REAL);
	m_coef              = Pmat * m_coef;
	MatrixXd     Rtrunc = R.topLeftCorner(m_r, m_r);
	MatrixXd  Rinvtrunc = Rtrunc.triangularView<Eigen::Upper>().solve(MatrixXd::Identity(m_r, m_r));
	m_unsc.topLeftCorner(m_r, m_r) = Rtrunc.setZero().selfadjointView<Eigen::Upper>().rankUpdate(Rinvtrunc);
	m_unsc.rightCols(nsing).fill(NA_REAL);
	m_unsc.bottomRows(nsing).fill(NA_REAL);
    } else {			// full rank
	MatrixXd       Rinv = R.triangularView<Eigen::Upper>().solve(MatrixXd::Identity(m_p, m_p));
	m_unsc              = R.setZero().selfadjointView<Eigen::Upper>().rankUpdate(Rinv);
	m_coef              = PQR.solve(y); 
	m_fitted            = X * m_coef;
    }
}

QR::QR(const MMatrixXd &X, const MVectorXd &y) : lm(X, y) {
    QRType     QR(X);
    m_coef        = QR.solve(y);
    m_fitted      = X * m_coef;
    MatrixXd    R = QR.matrixQR().topLeftCorner(m_p, m_p);
    MatrixXd Rinv = R.triangularView<Eigen::Upper>().solve(MatrixXd::Identity(m_p, m_p));
    m_unsc        = R.setZero().selfadjointView<Eigen::Upper>().rankUpdate(Rinv);
    for (Index i = 0; i < m_p; ++i) m_perm[i] = i;
}

LLT::LLT(const MMatrixXd &X, const MVectorXd &y) : lm(X, y) {
    LLTType  Ch(m_unsc.setZero().selfadjointView<Eigen::Lower>().rankUpdate(X.adjoint()));
    m_coef      = Ch.solve(X.adjoint() * y);
    m_fitted    = X * m_coef;
    m_unsc      = Ch.solve(MatrixXd::Identity(m_p, m_p));
    for (Index i = 0; i < m_p; ++i) m_perm[i] = i;
}

LDLT::LDLT(const MMatrixXd &X, const MVectorXd &y) : lm(X, y) {
    LDLTType Ch(m_unsc.setZero().selfadjointView<Eigen::Lower>().rankUpdate(X.adjoint()));
    m_coef      = Ch.solve(X.adjoint() * y);
    for (Index i = 0; i < m_p; ++i) m_perm[i] = i; // for now.  Maybe add a boolean flag to indicate if unsc is already reordered
//    m_perm      = PermutationType(Ch.transpositionsP()).indices();
    m_fitted    = X * m_coef;
    m_unsc      = Ch.solve(MatrixXd::Identity(m_p, m_p));
}

enum {ColPivQR_t = 0, QR_t, LDLT_t, LLT_t, SVD_t, SymmEigen_t};

static inline lm do_lm(const MMatrixXd &X, const MVectorXd &y, int type)
    throw (std::invalid_argument) {
    switch(type) {
    case ColPivQR_t:
	return ColPivQR(X, y);
    case QR_t:
	return QR(X, y);
    case LLT_t:
	return LLT(X, y);
    case LDLT_t:
	return LDLT(X, y);
    // case SVD_t:
    // 	return SVD(X, y);
    // case SymmEigen_t:
    // 	return SymmEigen(X, y);
    }
    throw std::invalid_argument("invalid type");
    return ColPivQR(X, y);	// -Wall
}

extern "C" SEXP fastLm(SEXP Xs, SEXP ys, SEXP type) {
    try {
	const NumericMatrix    X(Xs);
	const NumericVector    y(ys);
	Index    n = X.nrow(), p = X.ncol();
	if ((Index)y.size() != n)
	    throw std::invalid_argument("size mismatch");
	const MVectorXd       yy(y.begin(), n);
	lm                   ans = do_lm(MMatrixXd(X.begin(), n, p), yy, as<int>(type));
	
	NumericVector       coef(ans.coef().data(), ans.coef().data() + p);
				// install the names, if available
	List            dimnames = X.attr("dimnames");
	RObject         colnames = dimnames[1];
	if (!(colnames).isNULL())
	    coef.attr("names") = clone(CharacterVector(colnames));

	VectorXd           resid = yy - ans.fitted();
	double                s2 = resid.squaredNorm()/ans.df();
	PermutationType     Pmat = PermutationType(p);
	Pmat.indices()           = ans.perm();
	VectorXd              dd = Pmat * ans.unsc().diagonal();
	ArrayXd               se = (dd.array() * s2).sqrt();

	return List::create(_["coefficients"]  = coef,
			    _["se"]            = se,
			    _["rank"]          = ans.rank(),
			    _["df.residual"]   = ans.df(),
			    _["perm"]          = IntegerVector(ans.perm().data(), ans.perm().data() + p),
			    _["residuals"]     = resid,
			    _["s2"]            = s2,
			    _["fitted.values"] = ans.fitted(),
			    _["unsc"]          = ans.unsc());

    } catch( std::exception &ex ) {
	forward_exception_to_r( ex );
    } catch(...) { 
	::Rf_error( "c++ exception (unknown reason)" ); 
    }
    return R_NilValue; // -Wall
}

#if 0
extern "C" SEXP fastLmBench(SEXP Xs, SEXP ys) {
    using namespace Eigen;
    using namespace Rcpp;
    try {
	NumericMatrix X(Xs);
	NumericVector y(ys);
	Indices       n = X.nrow(), p = X.ncol();
	int          df = n - p;
	if ((size_t)y.size() != n)
	    throw std::invalid_argument("size mismatch");
	
	MatrixXd      A = Map<MatrixXd>(X.begin(), n, p); // shares storage
	VectorXd      b = Map<VectorXd>(y.begin(), n);

	HouseholderQR<MatrixXd> Aqr(A);
	VectorXd   coef = Aqr.solve(b);
	double        s = std::sqrt((b - A*coef).squaredNorm()/df);

	MatrixXd   Rinv((TriangularView<MatrixXd, Upper>(Aqr.matrixQR().topRows(p)))
			.solve(MatrixXd::Identity(p, p)));
	VectorXd     se = Rinv.rowwise().norm() * s;

	return List::create(_["coefficients"] = coef,
			    _["se"]           = se,
			    _["df.residual"]  = df);
    } catch( std::exception &ex ) {
	forward_exception_to_r( ex );
    } catch(...) { 
	::Rf_error( "c++ exception (unknown reason)" ); 
    }
    return R_NilValue; // -Wall
}

extern "C" SEXP fastLmChol1(SEXP Xs, SEXP ys) {
    using namespace Eigen;
    using namespace Rcpp;
    try {
	NumericMatrix X(Xs);
	NumericVector y(ys);
	size_t          n = X.nrow(), p = X.ncol();
	ptrdiff_t      df = n - p;
	if (df <= 0l)
	    throw std::invalid_argument("nrow(X) > ncol(X) not satisfied");
	if ((size_t)y.size() != n)
	    throw std::invalid_argument("size mismatch");
	
	MatrixXd        A = Map<MatrixXd>(X.begin(), n, p); // shares storage
	VectorXd        b = Map<VectorXd>(y.begin(), n);

	LLT<MatrixXd>  Ch(SelfAdjointView<MatrixXd, Lower>(MatrixXd::Zero(p, p)).rankUpdate(A.adjoint()));
	VectorXd     coef = Ch.solve(A.adjoint() * b);
	double          s = std::sqrt((b - A*coef).squaredNorm()/df);

	VectorXd       se = (Ch.matrixL().solve(MatrixXd::Identity(p, p))).colwise().norm() * s;;

	return List::create(_["coefficients"] = coef,
			    _["se"]           = se,
			    _["df.residual"]  = df);
    } catch( std::exception &ex ) {
	forward_exception_to_r( ex );
    } catch(...) { 
	::Rf_error( "c++ exception (unknown reason)" ); 
    }
    return R_NilValue; // -Wall
}

extern "C" SEXP fastLmChol2(SEXP Xs, SEXP ys) {
    using namespace Eigen;
    using namespace Rcpp;
    try {
	NumericMatrix   X(Xs);
	NumericVector   y(ys);
	size_t          n = X.nrow(), p = X.ncol();
	ptrdiff_t      df = n - p;
	if ((size_t)y.size() != n)
	    throw std::invalid_argument("size mismatch");
	
	MatrixXd        A = Map<MatrixXd>(X.begin(), n, p); // shares storage
	VectorXd        b = Map<VectorXd>(y.begin(), n);

	LDLT           Ch(MatrixXd::Zero(p, p).selfAdjointView(Eigen::Lower).rankUpdate(A.adjoint()));
	VectorXd     coef = Ch.solve(A.adjoint() * b);
	double         s2 = (b - A*coef).squaredNorm()/df;

	ArrayXd        se = (Ch.solve(MatrixXd::Identity(p, p)).diagonal().array() * s2).sqrt();
	NumericVector Rse(p);	// should define a wrap method for ArrayXd, ArrayXXd, etc.
	std::copy(se.data(), se.data() + p, Rse.begin());
			    
	return List::create(_["coefficients"] = coef,
			    _["se"]           = Rse,
			    _["df.residual"]  = df);
    } catch( std::exception &ex ) {
	forward_exception_to_r( ex );
    } catch(...) { 
	::Rf_error( "c++ exception (unknown reason)" ); 
    }
    return R_NilValue; // -Wall
}

extern "C" SEXP fastLmSVD(SEXP Xs, SEXP ys) {
    using namespace Eigen;
    using namespace Rcpp;
    try {
	NumericMatrix   X(Xs);
	NumericVector   y(ys);
	size_t          n = X.nrow(), p = X.ncol();
	ptrdiff_t      df = n - p;
	if ((size_t)y.size() != n)
	    throw std::invalid_argument("size mismatch");
	
	MatrixXd        A = Map<MatrixXd>(X.begin(), n, p); // shares storage
	VectorXd        b = Map<VectorXd>(y.begin(), n);

	VectorXd     coef = A.jacobiSvd(ComputeThinU|ComputeThinV).solve(b);
//	double         s2 = (b - A*coef).squaredNorm()/df;

	// ArrayXd        se = (Ch.solve(MatrixXd::Identity(p, p)).diagonal().array() * s2).sqrt();
	// NumericVector Rse(p);	// should define a wrap method for ArrayXd, ArrayXXd, etc.
	// std::copy(se.data(), se.data() + p, Rse.begin());
			    
	return List::create(_["coefficients"] = coef,
			    // _["se"]           = Rse,
			    _["df.residual"]  = df);
    } catch( std::exception &ex ) {
	forward_exception_to_r( ex );
    } catch(...) { 
	::Rf_error( "c++ exception (unknown reason)" ); 
    }
    return R_NilValue; // -Wall
}

#endif
