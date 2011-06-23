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
				// Basic matrix, vector and array types for double
typedef Eigen::MatrixXd                                MatrixXd;
typedef Eigen::VectorXd                                VectorXd;
typedef Eigen::ArrayXd                                 ArrayXd;
typedef Eigen::ArrayXXd                                ArrayXXd;
				// integer
typedef Eigen::MatrixXi                                MatrixXi;
typedef Eigen::VectorXi                                VectorXi;
typedef Eigen::ArrayXi                                 ArrayXi;
typedef Eigen::ArrayXXi                                ArrayXXi;
				// complex
typedef Eigen::MatrixXcd                               MatrixXcd;
typedef Eigen::VectorXcd                               VectorXcd;
typedef Eigen::ArrayXcd                                ArrayXcd;
typedef Eigen::ArrayXXcd                               ArrayXXcd;
				// these should be defined for each base type - we use double
typedef MatrixXd::Index                                Index;
typedef MatrixXd::Scalar                               Scalar;
typedef MatrixXd::RealScalar                           RealScalar;
				// Mapped matrix and vector types (share storage)
typedef Eigen::Map<MatrixXd>                           MMatrixXd;
typedef Eigen::Map<VectorXd>                           MVectorXd;
				// Views
typedef Eigen::TriangularView<MatrixXd, Eigen::Upper>  UpperTri;
typedef Eigen::TriangularView<MatrixXd, Eigen::Lower>  LowerTri;
typedef Eigen::SelfAdjointView<MatrixXd, Eigen::Upper> UpperSym;
typedef Eigen::SelfAdjointView<MatrixXd, Eigen::Lower> LowerSym;
				// Decomposition types
typedef Eigen::LLT<MatrixXd>                           LLTType;
typedef Eigen::LDLT<MatrixXd>                          LDLTType;
typedef Eigen::ColPivHouseholderQR<MatrixXd>           PivQRType;
typedef Eigen::HouseholderQR<MatrixXd>                 QRType;
typedef Eigen::JacobiSVD<MatrixXd>                     SVDType;
typedef Eigen::DiagonalMatrix<double, Eigen::Dynamic>  DiagType;
typedef Eigen::SelfAdjointEigenSolver<MatrixXd>        EigvType;
				// Types derived from decompositions
typedef PivQRType::PermutationType                     PermutationType;
typedef PermutationType::IndicesType                   IndicesType;

class lm {
protected:
    Index        m_n;	   /**< number of rows */
    Index        m_p;	   /**< number of columns */
    VectorXd     m_coef;   /**< coefficient vector */
    int          m_r;	   /**< computed rank or NA_INTEGER */
    int          m_df;	   /**< residual degrees of freedom */
    IndicesType  m_perm;   /**< column permutation */
    VectorXd     m_fitted; /**< vector of fitted values */
    MatrixXd     m_unsc;   /**< unscaled variance-covariance matrix */
    RealScalar   m_prescribedThreshold;	/**< user specified tolerance */
    bool         m_usePrescribedThreshold;
public:
    lm(const MMatrixXd&, const MVectorXd&);
    lm&        setThreshold(const RealScalar&); // patterned after ColPivHouseholderQR
    RealScalar    threshold() const;
    const VectorXd&    coef() const {return m_coef;}
    int                rank() const {return m_r;}
    int                  df() const {return m_df;}
    const IndicesType& perm() const {return m_perm;}
    const VectorXd&  fitted() const {return m_fitted;}
    const MatrixXd&    unsc() const {return m_unsc;}
};

class ColPivQR : public lm {
public:
    ColPivQR(const MMatrixXd&, const MVectorXd&);
};

class LLT : public lm {
public:
    LLT(const MMatrixXd&, const MVectorXd&);
};

class LDLT : public lm {
public:
    LDLT(const MMatrixXd&, const MVectorXd&);
};

class QR : public lm {
public:
    QR(const MMatrixXd&, const MVectorXd&);
};

class SVD : public lm {
public:
    SVD(const MMatrixXd&, const MVectorXd&);
};

class SymmEigen : public lm {
public:
    SymmEigen(const MMatrixXd&, const MVectorXd&);
};

extern "C" SEXP fastLm(SEXP Xs, SEXP ys, SEXP types);

extern "C" SEXP crossprod(SEXP Xs);

extern "C" SEXP tcrossprod(SEXP Xs);

#endif

