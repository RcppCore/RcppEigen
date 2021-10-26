#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP RcppEigen_Eigen_SSE();
extern SEXP RcppEigen_eigen_version(SEXP);
extern SEXP RcppEigen_fastLm_Impl(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"RcppEigen_Eigen_SSE",     (DL_FUNC) &RcppEigen_Eigen_SSE,     0},
    {"RcppEigen_eigen_version", (DL_FUNC) &RcppEigen_eigen_version, 1},
    {"RcppEigen_fastLm_Impl",   (DL_FUNC) &RcppEigen_fastLm_Impl,   3},
    {NULL, NULL, 0}
};

void R_init_RcppEigen(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
