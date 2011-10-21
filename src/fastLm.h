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
    using Eigen::HouseholderQR;
    using Eigen::JacobiSVD;
    using Eigen::LDLT;
    using Eigen::LLT;
    using Eigen::Lower;
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::SelfAdjointEigenSolver;
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
	Index      m_n;		/**< number of rows */
	Index      m_p;		/**< number of columns */
	VectorXd   m_coef;	/**< coefficient vector */
	int        m_r;		/**< computed rank or NA_INTEGER */
	int        m_df;	/**< residual degrees of freedom */
	Indices    m_perm;	/**< column permutation */
	VectorXd   m_fitted;	/**< vector of fitted values */
	MatrixXd   m_unsc;	/**< unscaled variance-covariance  */
	RealScalar m_prescribedThreshold; /**< user specified tolerance */
	bool       m_usePrescribedThreshold;
    public:
	lm(const Map<MatrixXd>&, const Map<VectorXd>&);
	lm&        setThreshold(const RealScalar&); // patterned after ColPivHouseholderQR code
	RealScalar    threshold() const;
	const VectorXd&    coef() const {return m_coef;}
	int                rank() const {return m_r;}
	int                  df() const {return m_df;}
	const Indices&     perm() const {return m_perm;}
	const VectorXd&  fitted() const {return m_fitted;}
	const MatrixXd&    unsc() const {return m_unsc;}
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

