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

namespace lmsol {
    using Rcpp::_;
    using Rcpp::as;
    using Rcpp::CharacterVector;
    using Rcpp::clone;
    using Rcpp::List;
    using Rcpp::NumericMatrix;
    using Rcpp::NumericVector;
    using Rcpp::RObject;
    using Rcpp::wrap;

    using std::invalid_argument;
    using std::numeric_limits;

    lm::lm(const Map<MatrixXd> &X, const Map<VectorXd> &y)
	: m_X(X),
	  m_y(y),
	  m_n(X.rows()),
	  m_p(X.cols()),
	  m_coef(VectorXd::Constant(m_p, ::NA_REAL)),
	  m_r(::NA_INTEGER),
	  m_df(m_n - m_p),
	  m_fitted(m_n),
	  m_se(VectorXd::Constant(m_p, ::NA_REAL)),
	  m_usePrescribedThreshold(false) {
    }

    lm& lm::setThreshold(const RealScalar& threshold) {
	m_usePrescribedThreshold = true;
	m_prescribedThreshold = threshold;
	return *this;
    }

    SelfAdjointView<MatrixXd,Lower> lm::XtX() const {
	return MatrixXd(m_p, m_p).setZero().selfadjointView<Lower>().
	    rankUpdate(m_X.adjoint());
    }
// For some reason that function returning a SelfAdjointView encounters a bad_alloc error
// Use a macro for the time being
#define XtX MatrixXd(m_p, m_p).setZero().selfadjointView<Lower>().rankUpdate(m_X.adjoint())
    /** Returns the threshold that will be used by certain methods such as rank().
     * 
     *  The default value comes from experimenting (see "LU precision
     *  tuning" thread on the Eigen list) and turns out to be
     *  identical to Higham's formula used already in LDLt. 
     *
     *  @return The user-prescribed threshold or the default.
     */
    RealScalar lm::threshold() const {
	return m_usePrescribedThreshold ? m_prescribedThreshold
	    : numeric_limits<double>::epsilon() * m_p; 
    }

    ColPivQR::ColPivQR(const Map<MatrixXd> &X, const Map<VectorXd> &y)
	: lm(X, y) {
	ColPivHouseholderQR<MatrixXd> PQR(X); // decompose the model matrix
	Permutation                  Pmat(PQR.colsPermutation());
	m_r                               = PQR.rank();
	if (m_r == m_p) {	// full rank case
	    m_coef     = PQR.solve(y);
	    m_fitted   = X * m_coef;
	    m_se       = Pmat * PQR.matrixQR().topRows(m_p).
		triangularView<Upper>().solve(I_p()).rowwise().norm();
	    return;
	} 
	m_df                              = m_n - m_r;
	MatrixXd                     Rinv(PQR.matrixQR().topLeftCorner(m_r, m_r).
					  triangularView<Upper>().
					  solve(MatrixXd::Identity(m_r, m_r)));
	VectorXd                  effects(PQR.householderQ().adjoint() * y);
	m_coef.head(m_r)                  = Rinv * effects.head(m_r);
	m_coef                            = Pmat * m_coef;
				// create fitted values from effects
				// (can't use X*m_coef if X is rank-deficient)
	effects.tail(m_df).setZero();
	m_fitted                          = PQR.householderQ() * effects;
	m_se.head(m_r)                    = Rinv.rowwise().norm();
	m_se                              = Pmat * m_se;
    }
    
    QR::QR(const Map<MatrixXd> &X, const Map<VectorXd> &y) : lm(X, y) {
	HouseholderQR<MatrixXd> QR(X);
	m_coef                     = QR.solve(y);
	m_fitted                   = X * m_coef;
	m_se                       = QR.matrixQR().topRows(m_p).
	    triangularView<Upper>().solve(I_p()).rowwise().norm();
    }
    
    
    Llt::Llt(const Map<MatrixXd> &X, const Map<VectorXd> &y) : lm(X, y) {
	LLT<MatrixXd>  Ch(XtX);
	m_coef            = Ch.solve(X.adjoint() * y);
	m_fitted          = X * m_coef;
	m_se              = Ch.matrixL().solve(I_p()).colwise().norm();
    }
    
    inline DiagonalMatrix<double, Dynamic> Dplus(const ArrayXd& D,
						 Index r, bool rev=false) {
	VectorXd   Di(VectorXd::Constant(D.size(), 0.));
	if (rev) Di.tail(r)  = D.tail(r).inverse();
	else Di.head(r)      = D.head(r).inverse();
	return DiagonalMatrix<double, Dynamic>(Di);
    }

    Ldlt::Ldlt(const Map<MatrixXd> &X, const Map<VectorXd> &y) : lm(X, y) {
	LDLT<MatrixXd> Ch(XtX);
	ArrayXd         D(Ch.vectorD());
	m_r               = (D > D.maxCoeff() * threshold()).count();
	// FIXME: work out how to use Dplus with elements of D unsorted.
	m_coef            = Ch.solve(X.adjoint() * y);
	m_fitted          = X * m_coef;
	m_se              = Ch.solve(I_p()).diagonal().array().sqrt();
    }
    
    SVD::SVD(const Map<MatrixXd> &X, const Map<VectorXd> &y) : lm(X, y) {
	JacobiSVD<MatrixXd>  UDV(X.jacobiSvd(ComputeThinU|ComputeThinV));
	ArrayXd                D(UDV.singularValues());
	m_r                      = (D > D[0] * threshold()).count();
	m_df                     = m_n - m_r;
	MatrixXd             VDi(UDV.matrixV() * Dplus(D, m_r));
	m_coef                   = VDi * UDV.matrixU().adjoint() * y;
	m_fitted                 = X * m_coef;
	m_se                     = VDi.rowwise().norm();
    }

    SymmEigen::SymmEigen(const Map<MatrixXd> &X, const Map<VectorXd> &y)
	: lm(X, y) {
	SelfAdjointEigenSolver<MatrixXd> eig(XtX);
	ArrayXd                      D(eig.eigenvalues());
	m_r                            = (D > D[m_p - 1] * threshold()).count();
	D                              = D.sqrt();
	m_df                           = m_n - m_r;
	MatrixXd                   VDi(eig.eigenvectors() * Dplus(D, m_r, true));
	m_coef                         = VDi * VDi.adjoint() * X.adjoint() * y;
	m_fitted                       = X * m_coef;
	m_se                           = VDi.rowwise().norm();
    }

    enum {ColPivQR_t = 0, QR_t, LLT_t, LDLT_t, SVD_t, SymmEigen_t};

    static inline lm do_lm(const Map<MatrixXd> &X, const Map<VectorXd> &y, int type) {
	switch(type) {
	case ColPivQR_t:
	    return ColPivQR(X, y);
	case QR_t:
	    return QR(X, y);
	case LLT_t:
	    return Llt(X, y);
	case LDLT_t:
	    return Ldlt(X, y);
	case SVD_t:
	    return SVD(X, y);
	case SymmEigen_t:
	    return SymmEigen(X, y);
	}
	throw invalid_argument("invalid type");
	return ColPivQR(X, y);	// -Wall
    }

    extern "C" SEXP fastLm(SEXP Xs, SEXP ys, SEXP type) {
	try {
	    const Map<MatrixXd>  X(as<Map<MatrixXd> >(Xs));
	    const Map<VectorXd>  y(as<Map<VectorXd> >(ys));
	    Index                n = X.rows();
	    if ((Index)y.size() != n) throw invalid_argument("size mismatch");

				// Select and apply the least squares method
	    lm                 ans(do_lm(X, y, ::Rf_asInteger(type)));

				// Copy coefficients and install names, if any
	    NumericVector     coef(wrap(ans.coef()));
	    List          dimnames(NumericMatrix(Xs).attr("dimnames"));
	    if (dimnames.size() > 1) {
		RObject   colnames = dimnames[1];
		if (!(colnames).isNULL())
		    coef.attr("names") = clone(CharacterVector(colnames));
	    }
	    
	    VectorXd         resid = y - ans.fitted();
	    double               s = resid.norm() / std::sqrt(double(ans.df()));
				// Create the standard errors
	    VectorXd            se = s * ans.se();

	    return List::create(_["coefficients"]  = coef,
				_["se"]            = se,
				_["rank"]          = ans.rank(),
				_["df.residual"]   = ans.df(),
				_["residuals"]     = resid,
				_["s"]             = s,
				_["fitted.values"] = ans.fitted());

	} catch( std::exception &ex ) {
	    forward_exception_to_r( ex );
	} catch(...) { 
	    ::Rf_error( "c++ exception (unknown reason)" ); 
	}
	return R_NilValue; // -Wall
    }
}
