// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-
//
// RcppEigenWrap.h: Rcpp wrap methods for Eigen matrices, vectors and arrays
//
// Copyright (C)      2011 Douglas Bates, Dirk Eddelbuettel and Romain Francois
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

#ifndef RcppEigen__RcppEigenWrap__h
#define RcppEigen__RcppEigenWrap__h

namespace Rcpp{

    namespace RcppEigen{

		template <typename T>
		SEXP Eigen_wrap( const T& object, const ::Rcpp::Dimension& dim){
			::Rcpp::RObject x = ::Rcpp::wrap( object.data() , object.data() + object.size() ) ;
			x.attr( "dim" ) = dim ;
			return x; 
		}

    } /* namespace RcppEigen */
	
    /* wrap */

    template <typename T>
	SEXP wrap(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& data) {
		return RcppEigen::Eigen_wrap(data, Dimension(data.rows(), data.cols()));
	}
    
	template <typename T>
	SEXP wrap(const Eigen::Matrix<T, Eigen::Dynamic, 1>& object ){
		return ::Rcpp::wrap(object.data(), object.data() + object.size());
    }

    template <typename T>
	SEXP wrap( const Eigen::Matrix<T, 1, Eigen::Dynamic>& data ){
		return RcppEigen::Eigen_wrap(data, Dimension(1, data.size()));
    }

    template <typename T>
	SEXP wrap(const Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>& data) {
		return RcppEigen::Eigen_wrap(data, Dimension(data.rows(), data.cols()));
	}
    
	template <typename T>
	SEXP wrap(const Eigen::Array<T, Eigen::Dynamic, 1>& object ){
		return ::Rcpp::wrap(object.data(), object.data() + object.size());
    }

#if 0
    /* support for Rcpp::as */
	
    namespace traits {
		
		template <typename T> 
		class Exporter< Eigen::Matrix<T, Eigen::Dynamic, 1> >
			: public IndexingExporter< Eigen::Matrix<T, Eigen::Dynamic, 1>, T > {
		public: 
			Exporter(SEXP x) : IndexingExporter< Eigen::Matrix<T, Eigen::Dynamic, 1>, T >(x){}
		}; 
		
		template <typename T> 
		class Exporter< Eigen::Matrix<T, 1, Eigen::Dynamic> >
			: public IndexingExporter< Eigen::Matrix<T, 1, Eigen::Dynamic>, T > {
		public:
			Exporter(SEXP x) : IndexingExporter< Eigen::Matrix<T, 1, Eigen::Dynamic>, T >(x){}
		}; 
		
		template <typename T> 
		class Exporter< Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >
			: public MatrixExporter< Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>, T > {
		public:
			Exporter(SEXP x) :
				MatrixExporter< Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>, T >(x){}
		}; 
		
    } // namespace traits
#endif
}

#endif

