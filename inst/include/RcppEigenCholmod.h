// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// RcppEigenCholmod.h: Provide access to the Matrix API and in turn
// to Eigen's CholmodSupport module.  Use of this header relies on
// compilation of ../../src/RcppEigenStubs.cpp and LinkingTo: Matrix.
//
// Copyright (C)      2011 Douglas Bates, Martin Maechler, Dirk Eddelbuettel and Romain Francois
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
// along with RcppEigen.  If not, see <http://www.gnu.org/licenses/>.

#ifndef RcppEigen__RcppEigenCholmod__h
#define RcppEigen__RcppEigenCholmod__h

#include <Matrix.h>
#ifndef R_MATRIX_CHOLMOD /* Matrix <= 1.6-1.1 */
# define R_MATRIX_CHOLMOD(_NAME_) M_cholmod_ ## _NAME_
# define M_cholmod_start M_R_cholmod_start /* sigh */
#endif

#include <Eigen/CholmodSupport>

namespace Rcpp {

	namespace RcppEigen {

	template<typename T>
	SEXP Eigen_cholmod_wrap(const Eigen::CholmodDecomposition<Eigen::SparseMatrix<T>> &obj) {
		const cholmod_factor *f = obj.factor();
		if (f->minor < f->n)
			throw std::runtime_error("CHOLMOD factorization was unsuccessful");
		S4 ans(std::string((f->is_super) ? "dCHMsuper" : "dCHMsimpl"));
		IntegerVector dd(2);
		IntegerVector tt((f->is_super) ? 6 : 4);
		dd[0] = dd[1] = f->n;
		tt[0] = f->ordering;
		tt[1] = f->is_ll;
		tt[2] = f->is_super;
		tt[3] = f->is_monotonic;
		ans.slot("Dim") = dd;
		ans.slot("type") = tt;
		ans.slot("colcount") = ::Rcpp::wrap((int *) f->ColCount, (int *) f->ColCount + f->n);
		ans.slot("perm") = ::Rcpp::wrap((int *) f->Perm, (int *) f->Perm + f->n);
		if (f->is_super) {
		tt[4] = f->maxcsize;
		tt[5] = f->maxesize;
		ans.slot("super") = ::Rcpp::wrap((int *) f->super, (int *) f->super + f->nsuper + 1);
		ans.slot("pi") = ::Rcpp::wrap((int *) f->pi, (int *) f->pi + f->nsuper + 1);
		ans.slot("px") = ::Rcpp::wrap((int *) f->px, (int *) f->px + f->nsuper + 1);
		ans.slot("s") = ::Rcpp::wrap((int *) f->s, (int *) f->s + f->ssize);
		ans.slot("x") = ::Rcpp::wrap((T *) f->x, (T *) f->x + f->xsize);
		} else {
		ans.slot("nxt") = ::Rcpp::wrap((int *) f->next, (int *) f->next + f->n + 2);
		ans.slot("prv") = ::Rcpp::wrap((int *) f->prev, (int *) f->prev + f->n + 2);
		ans.slot("nz") = ::Rcpp::wrap((int *) f->nz, (int *) f->nz + f->n);
		ans.slot("p") = ::Rcpp::wrap((int *) f->p, (int *) f->p + f->n + 1);
		ans.slot("i") = ::Rcpp::wrap((int *) f->i, (int *) f->i + f->nzmax);
		ans.slot("x") = ::Rcpp::wrap((T *) f->x, (T *) f->x + f->nzmax);
		}
		return ::Rcpp::wrap(ans);
	}

    } /* namespace RcppEigen */

    template<typename T>
    SEXP wrap(const Eigen::CholmodDecomposition<Eigen::SparseMatrix<T>> &obj) {
        return RcppEigen::Eigen_cholmod_wrap(obj);
    }

} /* namespace Rcpp */

#endif
