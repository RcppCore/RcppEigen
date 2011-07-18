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

lm::lm(const MMatrixXd &X, const MVectorXd &y)
    : m_n(X.rows()),
      m_p(X.cols()),
      m_coef(m_p),
      m_r(NA_INTEGER),
      m_df(m_n - m_p),
      m_perm(m_p),
      m_fitted(m_n),
      m_unsc(m_p, m_p),
      m_usePrescribedThreshold(false) {
}

lm& lm::setThreshold(const RealScalar& threshold)
{
    m_usePrescribedThreshold = true;
    m_prescribedThreshold = threshold;
    return *this;
}

/** Returns the threshold that will be used by certain methods such as rank().
 *
 */
RealScalar lm::threshold() const
{
    return m_usePrescribedThreshold ? m_prescribedThreshold
// this formula comes from experimenting (see "LU precision tuning" thread on the list)
// and turns out to be identical to Higham's formula used already in LDLt.
	: std::numeric_limits<double>::epsilon() * m_p; 
}

ColPivQR::ColPivQR(const MMatrixXd &X, const MVectorXd &y) : lm(X, y) {
    PivQRType           PQR(X);	// decompose the model matrix
    PermutationType    Pmat = PQR.colsPermutation();
    m_perm                  = Pmat.indices();
    m_r                     = PQR.rank();
    MatrixXd              R = PQR.matrixQR().topLeftCorner(m_p, m_p);

    if (m_r < (int)m_p) {		// The rank-deficient case
	m_df                = m_n - m_r;
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
    m_unsc        = m_unsc.setZero().selfadjointView<Eigen::Upper>().rankUpdate(Rinv);
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

// SVD method
MatrixXd pseudoInverse(const MatrixXd& X, double tolerance) {
    SVDType  UDV = X.jacobiSvd(Eigen::ComputeThinU|Eigen::ComputeThinV);
    VectorXd   D = UDV.singularValues();
    double   tol = D[0] * tolerance;
				// There are better ways of doing this
    double  *Dpt = D.data();
    for (int i = 0; i < D.size(); ++i)
	Dpt[i] = Dpt[i] < tol ? 0. : 1/Dpt[i];
// Eigen2 code
//UDV.matrixV() * (D.cwise() > tol).select(D.cwise().inverse(), 0).
//    asDiagonal() * UDV.matrixU().adjoint();
    return UDV.matrixV() * D.asDiagonal() * UDV.matrixU().adjoint();
}

SVD::SVD(const MMatrixXd &X, const MVectorXd &y) : lm(X, y) {
    SVDType  UDV = X.jacobiSvd(Eigen::ComputeThinU|Eigen::ComputeThinV);
    VectorXd   D = UDV.singularValues();
    m_r          = (D.array() > threshold() * D[0]).count();
    m_coef       = UDV.solve(y);
    m_fitted     = X * m_coef;
    MatrixXd VDi = UDV.matrixV() * DiagType(UDV.singularValues().array().inverse().matrix());
    m_unsc       = m_unsc.selfadjointView<Eigen::Lower>().rankUpdate(VDi);
    for (Index i = 0; i < m_p; ++i) m_perm[i] = i;
}


SymmEigen::SymmEigen(const MMatrixXd &X, const MVectorXd &y) : lm(X, y) {
    EigvType eig(m_unsc.setZero().selfadjointView<Eigen::Lower>().rankUpdate(X.adjoint()));
    ArrayXd    D = eig.eigenvalues().array().sqrt(); // singular values of X
    MatrixXd VDi = eig.eigenvectors() * DiagType(D.inverse().matrix());
    m_coef       = VDi * VDi.adjoint() * X.adjoint() * y;
    m_r          = std::count_if(D.data(), D.data() + m_p,
				 std::bind2nd(std::greater<double>(), D[0] * threshold()));
    m_fitted     = X * m_coef;
    m_unsc       = m_unsc.setZero().selfadjointView<Eigen::Lower>().rankUpdate(VDi);
    for (Index i = 0; i < m_p; ++i) m_perm[i] = i;
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
    case SVD_t:
     	return SVD(X, y);
    case SymmEigen_t:
     	return SymmEigen(X, y);
    }
    throw std::invalid_argument("invalid type");
    return ColPivQR(X, y);	// -Wall
}

extern "C" SEXP fastLm(SEXP Xs, SEXP ys, SEXP type) {
    try {
	const MMatrixXd      X(as<MMatrixXd>(Xs));
	const MVectorXd      y(as<MVectorXd>(ys));
	Index                n = X.rows(), p = X.cols();
	if ((Index)y.size() != n)
	    throw std::invalid_argument("size mismatch");

	lm                 ans = do_lm(X, y, ::Rf_asInteger(type));
				// Copy coefficients and install names, if available
	NumericVector     coef = wrap(ans.coef());
	List          dimnames = NumericMatrix(Xs).attr("dimnames");
	if (dimnames.size() > 1) {
	    RObject   colnames = dimnames[1];
	    if (!(colnames).isNULL())
		coef.attr("names") = clone(CharacterVector(colnames));
	}
	    
	VectorXd         resid = y - ans.fitted();
	double              s2 = resid.squaredNorm()/ans.df();
				// Create the standard errors
	PermutationType   Pmat = PermutationType(p);
	Pmat.indices()         = ans.perm();
	VectorXd            dd = Pmat * ans.unsc().diagonal();
	ArrayXd             se = (dd.array() * s2).sqrt();

	return List::create(_["coefficients"]  = coef,
			    _["se"]            = se,
			    _["rank"]          = ans.rank(),
			    _["df.residual"]   = ans.df(),
			    _["perm"]          = ans.perm(),
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

extern "C" SEXP crossprod(SEXP Xs) {
    try {
	const NumericMatrix X(Xs);
	Index               n = X.nrow(), p = X.ncol();
	const MMatrixXd    Xe(X.begin(), n, p);
	MatrixXd          XtX(p, p);
	XtX                   = XtX.setZero().selfadjointView<Eigen::Lower>().rankUpdate(Xe.adjoint());
	return                  wrap(XtX);
    } catch( std::exception &ex ) {
	forward_exception_to_r( ex );
    } catch(...) { 
	::Rf_error( "c++ exception (unknown reason)" ); 
    }
    return R_NilValue; // -Wall
}

extern "C" SEXP crossprod1(SEXP Xs) {
    try {
	const NumericMatrix X(Xs);
	Index               n = X.nrow(), p = X.ncol();
	const MMatrixXd    Xe(X.begin(), n, p);
	MatrixXd          XtX(p, p);
	XtX                   = XtX.setZero().selfadjointView<Eigen::Lower>().rankUpdate(Xe.adjoint());
	NumericMatrix     ans(p, p);
	std::copy(XtX.data(), XtX.data() + XtX.size(), ans.begin());
	return                   ans;
    } catch( std::exception &ex ) {
	forward_exception_to_r( ex );
    } catch(...) { 
	::Rf_error( "c++ exception (unknown reason)" ); 
    }
    return R_NilValue; // -Wall
}

extern "C" SEXP tcrossprod(SEXP Xs) {
    try {
	const NumericMatrix X(Xs);
	Index    n = X.nrow(), p = X.ncol();
	MatrixXd          XXt(n, n);
	XXt        = XXt.setZero().selfadjointView<Eigen::Lower>().rankUpdate(MMatrixXd(X.begin(), n, p));
	return wrap(XXt);
    } catch( std::exception &ex ) {
	forward_exception_to_r( ex );
    } catch(...) { 
	::Rf_error( "c++ exception (unknown reason)" ); 
    }
    return R_NilValue; // -Wall
}

extern "C" SEXP tcrossprod1(SEXP Xs) {
    try {
	const NumericMatrix X(Xs);
	Index    n = X.nrow(), p = X.ncol();
	MatrixXd          XXt(n, n);
	XXt        = XXt.setZero().selfadjointView<Eigen::Lower>().rankUpdate(MMatrixXd(X.begin(), n, p));
	NumericMatrix     ans(n, n);
	std::copy(XXt.data(), XXt.data() + XXt.size(), ans.begin());
	return   ans;
    } catch( std::exception &ex ) {
	forward_exception_to_r( ex );
    } catch(...) { 
	::Rf_error( "c++ exception (unknown reason)" ); 
    }
    return R_NilValue; // -Wall
}

