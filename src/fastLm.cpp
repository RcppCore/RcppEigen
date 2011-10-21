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
	: m_n(X.rows()),
	  m_p(X.cols()),
	  m_coef(m_p),
	  m_r(NA_INTEGER),
	  m_df(m_n - m_p),
	  m_perm(m_p),
	  m_fitted(m_n),
	  m_unsc(MatrixXd::Identity(m_p, m_p)),
	  m_usePrescribedThreshold(false) {
	for (Index i = 0; i < m_p; ++i) m_perm[i] = i;
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
	    : numeric_limits<double>::epsilon() * m_p; 
    }

    inline MatrixXd AAt(const MatrixXd& A) {
	return MatrixXd(MatrixXd(A.cols(), A.cols()).
			setZero().selfadjointView<Lower>().rankUpdate(A));
    }

    ColPivQR::ColPivQR(const Map<MatrixXd> &X, const Map<VectorXd> &y) : lm(X, y) {
	ColPivHouseholderQR<MatrixXd> PQR(X); // decompose the model matrix
	Permutation                  Pmat(PQR.colsPermutation());
	m_perm                            = Pmat.indices();
	m_r                               = PQR.rank();
	MatrixXd                        R(PQR.matrixQR().topRows(m_p).triangularView<Upper>());
	
	if (m_r < (int)m_p) {	// X is rank-deficient 
	    m_df                           = m_n - m_r;
	    int                      nsing(m_p - m_r);
	    MatrixXd                Atrunc((X * Pmat).leftCols(m_r));
	    HouseholderQR<MatrixXd>     QR(Atrunc);
	    VectorXd            coefTrunc(QR.solve(y));
	    m_fitted                      = Atrunc * coefTrunc;
	    m_coef.fill(::NA_REAL);
	    m_coef.topRows(m_r)           = coefTrunc;
	    m_coef                        = Pmat * m_coef;
	    MatrixXd               Rtrunc(MatrixXd(R).topLeftCorner(m_r, m_r));
	    MatrixXd            Rinvtrunc(Rtrunc.triangularView<Upper>().
					  solve(MatrixXd::Identity(m_r, m_r)));
	    m_unsc.topLeftCorner(m_r,m_r) = AAt(Rinvtrunc);
	    m_unsc.rightCols(nsing).fill(NA_REAL);
	    m_unsc.bottomRows(nsing).fill(NA_REAL);
	} else {		// full rank X
	    MatrixXd                 Rinv(R.triangularView<Upper>().solve(m_unsc));
	    m_unsc                        = AAt(Rinv);
	    m_coef                        = PQR.solve(y);
	    m_fitted                      = X * m_coef;
	}
    }
    
    QR::QR(const Map<MatrixXd> &X, const Map<VectorXd> &y) : lm(X, y) {
	HouseholderQR<MatrixXd> QR(X);
	MatrixXd              Rinv(QR.matrixQR().topRows(m_p).triangularView<Upper>().solve(m_unsc));
	m_unsc                     = AAt(Rinv);
	m_coef                     = QR.solve(y);
	m_fitted                   = X * m_coef;
    }
    
    Llt::Llt(const Map<MatrixXd> &X, const Map<VectorXd> &y) : lm(X, y) {
	LLT<MatrixXd>  Ch(m_unsc.setZero().selfadjointView<Lower>().rankUpdate(X.adjoint()));
	m_coef            = Ch.solve(X.adjoint() * y);
	m_fitted          = X * m_coef;
	m_unsc            = Ch.solve(MatrixXd::Identity(m_p, m_p));
    }
    
    Ldlt::Ldlt(const Map<MatrixXd> &X, const Map<VectorXd> &y) : lm(X, y) {
	LDLT<MatrixXd> Ch(m_unsc.setZero().selfadjointView<Lower>().rankUpdate(X.adjoint()));
	m_coef            = Ch.solve(X.adjoint() * y);
	m_fitted          = X * m_coef;
	m_unsc            = Ch.solve(MatrixXd::Identity(m_p, m_p));
    }
    
    typedef Eigen::DiagonalMatrix<double, Eigen::Dynamic>  DiagXd;

    inline DiagXd Dplus(const ArrayXd& D, Index r) {
	VectorXd   Di(VectorXd::Constant(D.size(), 0.));
	for (Index i = 0; i < r; ++i) Di[i] = 1. / D[i];
	return DiagXd(Di);
    }

    SVD::SVD(const Map<MatrixXd> &X, const Map<VectorXd> &y) : lm(X, y) {
	JacobiSVD<MatrixXd>  UDV(X.jacobiSvd(Eigen::ComputeThinU|Eigen::ComputeThinV));
	ArrayXd                D(UDV.singularValues());
	m_r                      = (D > D[0] * threshold()).count();
	MatrixXd             VDi(UDV.matrixV() * Dplus(D, m_r));
	m_coef                   = VDi * UDV.matrixU().adjoint() * y;
	m_fitted                 = X * m_coef;
	m_unsc                   = AAt(VDi);
    }

    SymmEigen::SymmEigen(const Map<MatrixXd> &X, const Map<VectorXd> &y) : lm(X, y) {
	SelfAdjointEigenSolver<MatrixXd>
	    eig(m_unsc.setZero().selfadjointView<Lower>().rankUpdate(X.adjoint()));
	ArrayXd                            D(eig.eigenvalues().array().sqrt());
	m_r                                  = (D > D[0] * threshold()).count();
	MatrixXd                         VDi(eig.eigenvectors() * Dplus(D, m_r));
	m_unsc                               = AAt(VDi);
	m_coef                               = m_unsc * X.adjoint() * y;
	m_fitted                             = X * m_coef;
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
	    Index                n = X.rows(), p = X.cols();
	    if ((Index)y.size() != n)
		throw invalid_argument("size mismatch");
				// Select and apply the least squares method
	    lm                 ans(do_lm(X, y, ::Rf_asInteger(type)));
				// Copy coefficients
	    NumericVector     coef(wrap(ans.coef()));
				// Install names, if available
	    List          dimnames(NumericMatrix(Xs).attr("dimnames"));
	    if (dimnames.size() > 1) {
		RObject   colnames = dimnames[1];
		if (!(colnames).isNULL())
		    coef.attr("names") = clone(CharacterVector(colnames));
	    }
	    
	    VectorXd         resid = y - ans.fitted();
	    double              s2 = resid.squaredNorm()/ans.df();
				// Create the standard errors
	    Permutation       Pmat = Permutation(p);
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
}
