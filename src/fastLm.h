// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
#ifndef RCPPEIGEN_FASTLM_H
#define RCPPEIGEN_FASTLM_H

#include <RcppEigen.h>
				// Basic matrix, vector and array types
typedef Eigen::MatrixXd                                MatrixXd;
typedef Eigen::VectorXd                                VectorXd;
typedef Eigen::VectorXi                                VectorXi;
typedef Eigen::ArrayXd                                 ArrayXd;
typedef Eigen::ArrayXXd                                ArrayXXd;
typedef Eigen::ArrayXi                                 ArrayXi;
typedef Eigen::ArrayXXi                                ArrayXXi;
typedef typename MatrixXd::Index                       Index;
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
typedef Eigen::DiagonalMatrix<double,Eigen::Dynamic>   DiagType;
				// Types derived from decompositions
typedef typename PivQRType::PermutationType            PermutationType;
typedef typename PermutationType::IndicesType          IndicesType;

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
public:
    lm(const MMatrixXd&, const MVectorXd&);
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

extern "C" SEXP fastLm(SEXP Xs, SEXP ys, SEXP types);

#endif

