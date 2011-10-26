// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// RcppEigen.cpp: Rcpp/Eigen glue
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
// along with RcppEigen.  If not, see <http://www.gnu.org/licenses/>.

#include <RcppEigen.h>

extern "C" {
    SEXP eigen_version(SEXP single_){
	using Rcpp::_;
	using Rcpp::IntegerVector;
	
	BEGIN_RCPP;
	bool single = Rcpp::as<bool>(single_) ;
	if( single ){
	    return Rcpp::wrap( 10000 * EIGEN_WORLD_VERSION +
			       100 * EIGEN_MAJOR_VERSION + 
			       EIGEN_MINOR_VERSION ) ;
	}
	
	return IntegerVector::create(_["major"] = EIGEN_WORLD_VERSION,
				     _["minor"] = EIGEN_MAJOR_VERSION,
				     _["patch"] = EIGEN_MINOR_VERSION);
	END_RCPP;
    }

    SEXP Eigen_SSE() {
	BEGIN_RCPP;
	return Rcpp::wrap(Eigen::SimdInstructionSetsInUse());
	END_RCPP;
    }
}




