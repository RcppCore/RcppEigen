
library(inline)
library(RcppEigen)


## section 3.1
(A <- matrix(1:6, ncol=2))
str(A)

transCpp <-'
using Eigen::Map;
using Eigen::MatrixXi;
                 // Map the integer matrix AA from R
const Map<MatrixXi>  A(as<Map<MatrixXi> >(AA));
                 // evaluate and return the transpose of A
const MatrixXi      At(A.transpose());
return wrap(At);
'

ftrans <- cxxfunction(signature(AA="matrix"), transCpp, plugin="RcppEigen")
(At <- ftrans(A))
stopifnot(all.equal(At, t(A)))



## section 3.2
prodCpp <- '
using Eigen::Map;
using Eigen::MatrixXi;
const Map<MatrixXi>    B(as<Map<MatrixXi> >(BB));
const Map<MatrixXi>    C(as<Map<MatrixXi> >(CC));
return List::create(_["B %*% C"]         = B * C,
                    _["crossprod(B, C)"] = B.adjoint() * C);
'

fprod <- cxxfunction(signature(BB = "matrix", CC = "matrix"), prodCpp, "RcppEigen")
B <- matrix(1:4, ncol=2)
C <- matrix(6:1, nrow=2)
str(fp <- fprod(B, C))
stopifnot(all.equal(fp[[1]], B %*% C), all.equal(fp[[2]], crossprod(B, C)))



## section 3.3

crossprodCpp <- '
using Eigen::Map;
using Eigen::MatrixXi;
using Eigen::Lower;

const Map<MatrixXi> A(as<Map<MatrixXi> >(AA));
const int           m(A.rows()), n(A.cols());
MatrixXi          AtA(MatrixXi(n, n).setZero().
                      selfadjointView<Lower>().rankUpdate(A.adjoint()));
MatrixXi          AAt(MatrixXi(m, m).setZero().
                      selfadjointView<Lower>().rankUpdate(A));

return List::create(_["crossprod(A)"]  = AtA,
                    _["tcrossprod(A)"] = AAt);
'
fcprd <- cxxfunction(signature(AA = "matrix"), crossprodCpp, "RcppEigen")
str(crp <- fcprd(A))
stopifnot(all.equal(crp[[1]], crossprod(A)),
          all.equal(crp[[2]], tcrossprod(A)))



## section 3.4

cholCpp <- '
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::LLT;
using Eigen::Lower;

const Map<MatrixXd>   A(as<Map<MatrixXd> >(AA));
const int             n(A.cols());
const LLT<MatrixXd> llt(MatrixXd(n, n).setZero().
                        selfadjointView<Lower>().rankUpdate(A.adjoint()));

return List::create(_["L"] = MatrixXd(llt.matrixL()),
                    _["R"] = MatrixXd(llt.matrixU()));
'

fchol <- cxxfunction(signature(AA = "matrix"), cholCpp, "RcppEigen")
(ll <- fchol(A))
stopifnot(all.equal(ll[[2]], chol(crossprod(A))))


# section 3.5

cholDetCpp <- '
using Eigen::Lower;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

const Map<MatrixXd>   A(as<Map<MatrixXd> >(AA));
const int             n(A.cols());
const MatrixXd      AtA(MatrixXd(n, n).setZero().
                        selfadjointView<Lower>().rankUpdate(A.adjoint()));
const MatrixXd     Lmat(AtA.llt().matrixL());
const double       detL(Lmat.diagonal().prod());
const VectorXd     Dvec(AtA.ldlt().vectorD());

return List::create(_["d1"] = detL * detL,
                    _["d2"] = Dvec.prod(),
                    _["ld"] = Dvec.array().log().sum());
'

fdet <- cxxfunction(signature(AA = "matrix"), cholDetCpp, "RcppEigen")
unlist(ll <- fdet(A))

