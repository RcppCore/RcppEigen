// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// RcppEigen.cpp: Rcpp/Armadillo glue
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

using namespace Rcpp;
extern "C" SEXP eigen_version(SEXP single_){

    bool single = as<bool>( single_) ;
    int major = 3, minor = 0, patch = 1;
    if( single ){
	return wrap( 10000*major +
		     100*minor + 
		     patch ) ;
    }

    IntegerVector version = 
	IntegerVector::create(_["major"] = major,
			      _["minor"] = minor,
			      _["patch"] = patch);

   return version ;

}


