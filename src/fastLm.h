// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// fastLm.h: Rcpp/Eigen example of a simple lm() alternative
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
#ifndef RCPPEIGEN_FASTLM_H
#define RCPPEIGEN_FASTLM_H

#include <RcppEigen.h>

namespace lmsol {
    using Eigen::ArrayXd;
    using Eigen::ColPivHouseholderQR;
    using Eigen::DiagonalMatrix;
    using Eigen::Dynamic;
    using Eigen::HouseholderQR;
    using Eigen::JacobiSVD;
    using Eigen::LDLT;
    using Eigen::LLT;
    using Eigen::Lower;
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::SelfAdjointEigenSolver;
    using Eigen::SelfAdjointView;
    using Eigen::TriangularView;
    using Eigen::VectorXd;
    using Eigen::Upper;

    typedef MatrixXd::Index                                 Index;
    typedef MatrixXd::Scalar                                Scalar;
    typedef MatrixXd::RealScalar                            RealScalar;
    typedef ColPivHouseholderQR<MatrixXd>::PermutationType  Permutation;
    typedef Permutation::IndicesType                        Indices;

    class lm {
    protected:
	Map<MatrixXd> m_X;	/**< model matrix */
	Map<VectorXd> m_y;	/**< response vector */
	Index         m_n;	/**< number of rows of X */
	Index         m_p;	/**< number of columns of X */
	VectorXd      m_coef;	/**< coefficient vector */
	int           m_r;	/**< computed rank or NA_INTEGER */
	int           m_df;	/**< residual degrees of freedom */
	VectorXd      m_fitted;	/**< vector of fitted values */
	VectorXd      m_se;	/**< standard errors  */
	RealScalar    m_prescribedThreshold; /**< user specified tolerance */
	bool          m_usePrescribedThreshold;
    public:
	lm(const Map<MatrixXd>&, const Map<VectorXd>&);

         // setThreshold and threshold are based on ColPivHouseholderQR methods
	MatrixXd                        I_p() const {
	    return MatrixXd::Identity(m_p, m_p);
	}
	RealScalar                threshold() const;
	SelfAdjointView<MatrixXd,Lower> XtX() const {
	    return MatrixXd(m_p, m_p).setZero().selfadjointView<Lower>().rankUpdate(m_X.adjoint());
	};
	const VectorXd&                  se() const {return m_se;}
	const VectorXd&                coef() const {return m_coef;}
	const VectorXd&              fitted() const {return m_fitted;}
	int                              df() const {return m_df;}
	int                            rank() const {return m_r;}
	lm&                    setThreshold(const RealScalar&); 
    };

    class ColPivQR : public lm {
    public:
	ColPivQR(const Map<MatrixXd>&, const Map<VectorXd>&);
    };

    class Llt : public lm {
    public:
	Llt(const Map<MatrixXd>&, const Map<VectorXd>&);
    };

    class Ldlt : public lm {
    public:
	Ldlt(const Map<MatrixXd>&, const Map<VectorXd>&); 
    };

    class QR : public lm {
    public:
	QR(const Map<MatrixXd>&, const Map<VectorXd>&);
    };

    class SVD : public lm {
    public:
	SVD(const Map<MatrixXd>&, const Map<VectorXd>&);
    };

    class SymmEigen : public lm {
    public:
	SymmEigen(const Map<MatrixXd>&, const Map<VectorXd>&);
    };
}

extern "C" SEXP fastLm(SEXP Xs, SEXP ys, SEXP types);

#endif

