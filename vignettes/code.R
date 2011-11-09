
library(inline)
library(RcppEigen)


incl <- '
using   Eigen::LLT;
using   Eigen::Lower;
using   Eigen::Map;
using   Eigen::MatrixXd;
using   Eigen::MatrixXi;
using   Eigen::Upper;
using   Eigen::VectorXd;
typedef Map<MatrixXd>  MapMatd;
typedef Map<MatrixXi>  MapMati;
typedef Map<VectorXd>  MapVecd;
inline MatrixXd AtA(const MapMatd& A) {
    int    n(A.cols());
    return   MatrixXd(n,n).setZero().selfadjointView<Lower>()
             .rankUpdate(A.adjoint());
}
'


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
typedef Eigen::Map<Eigen::MatrixXi>   MapMati;
const MapMati    B(as<MapMati>(BB));
const MapMati    C(as<MapMati>(CC));
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

storage.mode(A) <- "double"

cholCpp <- '
const  LLT<MatrixXd> llt(AtA(as<MapMatd>(AA)));
return List::create(_["L"] = MatrixXd(llt.matrixL()),
                    _["R"] = MatrixXd(llt.matrixU()));
'

fchol <- cxxfunction(signature(AA = "matrix"), cholCpp, "RcppEigen", incl)
(ll <- fchol(A))
stopifnot(all.equal(ll[[2]], chol(crossprod(A))))


# section 3.5

cholDetCpp <- '
const MatrixXd      ata(AtA(as<MapMatd>(AA)));
const MatrixXd     Lmat(ata.llt().matrixL());
const double       detL(Lmat.diagonal().prod());
const VectorXd     Dvec(ata.ldlt().vectorD());
return List::create(_["d1"] = detL * detL,
                    _["d2"] = Dvec.prod(),
                    _["ld"] = Dvec.array().log().sum());
'

fdet <- cxxfunction(signature(AA = "matrix"), cholDetCpp, "RcppEigen", incl)
unlist(ll <- fdet(A))


## section 4.1
lltLSCpp <- '
const MapMatd         X(as<MapMatd>(XX));
const MapVecd         y(as<MapVecd>(yy));
const int             n(X.rows()), p(X.cols());
const LLT<MatrixXd> llt(AtA(X));
const VectorXd  betahat(llt.solve(X.adjoint() * y));
const VectorXd   fitted(X * betahat);
const VectorXd    resid(y - fitted);
const int            df(n - p);
const double          s(resid.norm() / std::sqrt(double(df)));
const VectorXd       se(s * llt.matrixL().solve(MatrixXd::Identity(p, p))
                        .colwise().norm());
return     List::create(_["coefficients"]   = betahat,
                        _["fitted.values"]  = fitted,
                        _["residuals"]      = resid,
                        _["s"]              = s,
                        _["df.residual"]    = df,
                        _["rank"]           = p,
                        _["Std. Error"]     = se);
'

lltLS <- cxxfunction(signature(XX = "matrix", yy = "numeric"),
                     lltLSCpp, "RcppEigen", incl)
data(trees, package="datasets")
str(lltFit <- with(trees, lltLS(cbind(1, log(Girth)), log(Volume))))
str(lmFit <- with(trees, lm.fit(cbind(1, log(Girth)), log(Volume))))
for (nm in c("coefficients", "residuals", "fitted.values", "rank"))
    stopifnot(all.equal(lltFit[[nm]], unname(lmFit[[nm]])))
stopifnot(all.equal(lltFit[["Std. Error"]],
                    unname(coef(summary(lm(log(Volume) ~ log(Girth), trees)))[,2])))


## section 4.3

dd <- data.frame(f1 = gl(4, 6, labels = LETTERS[1:4]),
                 f2 = gl(3, 2, labels = letters[1:3]))[-(7:8), ]
xtabs(~ f2 + f1, dd)                    # one missing cell
mm <- model.matrix(~ f1 * f2, dd)
kappa(mm)         # large condition number, indicating rank deficiency
rcond(mm)         # alternative evaluation, the reciprocal condition number
(c(rank=qr(mm)$rank, p=ncol(mm))) # rank as computed in R's qr function
set.seed(1)
dd$y <- mm %*% seq_len(ncol(mm)) + rnorm(nrow(mm), sd = 0.1)
                         # lm detects the rank deficiency
fm1 <- lm(y ~ f1 * f2, dd)
writeLines(capture.output(print(summary(fm1), signif.stars=FALSE))[9:22])


## section 4.6
print(summary(fmPQR <- fastLm(y ~ f1 * f2, dd)), signif.stars=FALSE)
all.equal(coef(fm1), coef(fmPQR))
all.equal(unname(fitted(fm1)), fitted(fmPQR))
all.equal(unname(residuals(fm1)), residuals(fmPQR))


print(summary(fmSVD <- fastLm(y ~ f1 * f2, dd, method=4L)), signif.stars=FALSE)
all.equal(coef(fm1), coef(fmSVD))
all.equal(unname(fitted(fm1)), fitted(fmSVD))
all.equal(unname(residuals(fm1)), residuals(fmSVD))


print(summary(fmVLV <- fastLm(y ~ f1 * f2, dd, method=5L)), signif.stars=FALSE)
all.equal(coef(fmSVD), coef(fmVLV))
all.equal(unname(fitted(fm1)), fitted(fmSVD))
all.equal(unname(residuals(fm1)), residuals(fmSVD))



## section 5

badtransCpp <- '
const MapMati  A(as<MapMati>(AA));
return wrap(A.transpose());
'

Ai <- matrix(1:6, ncol=2L)
ftrans2 <- cxxfunction(signature(AA = "matrix"), badtransCpp, "RcppEigen")
(At <- ftrans2(Ai))
all.equal(At, t(Ai))



## section 6
sparseProdCpp <- '
using Eigen::Map;
using Eigen::MappedSparseMatrix;
using Eigen::SparseMatrix;
using Eigen::VectorXd;

const MappedSparseMatrix<double>  A(as<MappedSparseMatrix<double> >(AA));
const Map<VectorXd>               y(as<Map<VectorXd> >(yy));
const SparseMatrix<double>       At(A.adjoint());
return List::create(_["At"]  = At,
                    _["Aty"] = At * y);
'

sparse1 <- cxxfunction(signature(AA = "dgCMatrix", yy = "numeric"),
                       sparseProdCpp, "RcppEigen")
data(KNex, package="Matrix")
rr <- sparse1(KNex$mm, KNex$y)
stopifnot(all.equal(rr$At, t(KNex$mm)),
          all.equal(rr$Aty, as.vector(crossprod(KNex$mm, KNex$y))))


sparseLSCpp <- '
using   Eigen::Lower;
using   Eigen::VectorXd;
typedef Eigen::Map<VectorXd>               MapVec;
typedef Eigen::MappedSparseMatrix<double>  MSpMat;
typedef Eigen::SparseMatrix<double>         SpMat;
typedef Eigen::SimplicialLDLt<SpMat>       SpChol;
typedef Eigen::CholmodDecomposition<SpMat> CholMD;

const SpMat      At(as<MSpMat>(AA).adjoint());
const VectorXd  Aty(At * as<MapVec>(yy));
const SpChol     Ch(At * At.adjoint());
if (Ch.info() != Eigen::Success)
   return R_NilValue;
const CholMD      L(At);
if (L.info() != Eigen::Success)
   return R_NilValue;
return List::create(_["L"]        = wrap(L),
                    _["betahatS"] = Ch.solve(Aty),
                    _["betahatC"] = L.solve(Aty),
                    _["perm"]     = Ch.permutationP().indices());
'


sparse2 <- cxxfunction(signature(AA = "dgCMatrix", yy = "numeric"),
                       sparseLSCpp, "RcppEigen")
str(rr <-  sparse2(KNex$mm, KNex$y))
res <- as.vector(solve(Ch <- Cholesky(crossprod(KNex$mm)),
                       crossprod(KNex$mm, KNex$y)))
stopifnot(all.equal(rr$betahatS, res), all.equal(rr$betahatC, res))
all.equal(rr$L, Ch)   # not sure yet why these are different.  The one from Eigen is smaller
## It's because the Ch was created from mm'mm and L was created directly from mm'
all(rr$perm == Ch@perm) # fill-reducing permutations are different
