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
        
        // helper trait to identify if T is a plain object type
        // TODO: perhaps move this to its own file
        template <typename T> struct is_plain : Rcpp::traits::same_type<T,typename T::PlainObject>{} ;
        
        // plain object, so we can assume data() and size()
        template <typename T>
        SEXP eigen_wrap_is_plain( const T& obj, ::Rcpp::traits::true_type ){
            // FIXME: deal with RowMajor, etc ...
            const int RTYPE = Rcpp::traits::r_sexptype_traits<typename T::Scalar>::rtype ;
            if( obj.cols() == 1 ) {
                return wrap( obj.data(), obj.data() + obj.size() ) ;
            } else {
                Rcpp::Matrix<RTYPE> x( obj.rows(), obj.cols(), obj.data() ) ;
                return x; 
            }
        }
       
        // when the object is not plain, we need to eval()uate it
        template <typename T>
        inline SEXP eigen_wrap_is_plain( const T& obj, ::Rcpp::traits::false_type ){
            return eigen_wrap_is_plain( obj.eval(), Rcpp::traits::true_type() ) ;
        }
        
        
        // at that point we know that T derives from EigenBase
        // so it is either a plain object (Matrix, etc ...) or an expression
        // that eval()uates into a plain object
        //
        // so the first thing we need to do is to find out so that we don't evaluate if we don't need to
        template <typename T>
        inline SEXP eigen_wrap( const T& obj ){
            return eigen_wrap_is_plain( obj, 
                typename is_plain<T>::type() 
                ) ;
        }
        
        
		template <typename T>
		SEXP Eigen_wrap( const T& object, const ::Rcpp::Dimension& dim){
			::Rcpp::RObject x = ::Rcpp::wrap(object.data(), object.data() + object.size());
			x.attr( "dim" ) = dim ;
			return x; 
		}

    } /* namespace RcppEigen */

    // /* wrap */
    // [romain] : no longer necessary
    // template <typename Derived>
    // SEXP wrap(const Eigen::EigenBase<Derived>& object) {
    //     //FIXME: Check IsRowMajor and transpose if needed
    // 	::Rcpp::RObject x = ::Rcpp::wrap(object.data(), object.data() + object.size());
    // 	if (object.ColsAtCompileTime == 1) return x; // represented as a vector
    // 	x.attr("dim") = ::Rcpp::Dimension(object.rows(), object.cols());
    // }
    // 
	// template <typename T>
	// SEXP wrap(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& data) {
	// 	return RcppEigen::Eigen_wrap(data, Dimension(data.rows(), data.cols()));
	// }
    // 
	// template <typename T>
	// SEXP wrap(const Eigen::Matrix<T, Eigen::Dynamic, 1>& object ){
	// 	return ::Rcpp::wrap(object.data(), object.data() + object.size());
    // }
    // 
    // template <typename T>
	// SEXP wrap( const Eigen::Matrix<T, 1, Eigen::Dynamic>& data ){
	// 	return RcppEigen::Eigen_wrap(data, Dimension(1, data.size()));
    // }
    // 
    // template <typename T>
	// SEXP wrap(const Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>& data) {
	// 	return RcppEigen::Eigen_wrap(data, Dimension(data.rows(), data.cols()));
	// }
    // 
	// template <typename T>
	// SEXP wrap(const Eigen::Array<T, Eigen::Dynamic, 1>& object ){
	// 	return ::Rcpp::wrap(object.data(), object.data() + object.size());
    // }
    // 
    // template <typename T>
	// SEXP wrap(const Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& data) {
	// 	return RcppEigen::Eigen_wrap(data, Dimension(data.rows(), data.cols()));
	// }
    // 
	// template <typename T>
	// SEXP wrap(const Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1> >& object ){
	// 	return ::Rcpp::wrap(object.data(), object.data() + object.size());
    // }
    // 
    // template <typename T>
	// SEXP wrap(const Eigen::Map<Eigen::Matrix<T, 1, Eigen::Dynamic> >& data ){
	// 	return RcppEigen::Eigen_wrap(data, Dimension(1, data.size()));
    // }
    // 
    // template <typename T>
	// SEXP wrap(const Eigen::Map<Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> >& data) {
	// 	return RcppEigen::Eigen_wrap(data, Dimension(data.rows(), data.cols()));
	// }
    // 
	// template <typename T>
	// SEXP wrap(const Eigen::Map<Eigen::Array<T, Eigen::Dynamic, 1> >& object ){
	// 	return ::Rcpp::wrap(object.data(), object.data() + object.size());
    // }
               
    // we can probably deal with sparse stuff more generically
	template <typename T>
    SEXP wrap(const Eigen::Map<Eigen::SparseMatrix<T> >& object ) {
		int          nnz = object.nonZeros(), p = object.outerSize();
		Dimension    dim(object.innerSize(), p);
		const int    *ip = object._innerIndexPtr(), *pp = object._outerIndexPtr();
		const T      *xp = object._valuePtr();
		IntegerVector iv(ip, ip + nnz), pv(pp, pp + p + 1);
		NumericVector xv(xp, xp + nnz);
		
		return ::Rcpp::wrap(List::create(_["Dim"] = dim,
										 _["i"]   = iv,
										 _["p"]   = pv,
										 _["x"]   = xv));
	}

    /* support for Rcpp::as */
	
    namespace traits {

		template<typename T>
		class Exporter<Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1> > > {
		public:
			typedef typename Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1> >  MVType;
			Exporter(SEXP x) : d_size(::Rf_length(x)) {
				const int RTYPE = ::Rcpp::traits::r_sexptype_traits<T>::rtype ;
				if (TYPEOF(x) != RTYPE)
					throw std::invalid_argument("Wrong R type for mapped vector");
				typedef typename ::Rcpp::traits::storage_type<RTYPE>::type STORAGE;
				d_start         = ::Rcpp::internal::r_vector_start<RTYPE,STORAGE>(x);
			}
			MVType get() {return MVType(d_start, d_size);}
		protected:
			const int d_size;
			T*        d_start;
		};

		template<typename T>
		class Exporter<Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > > {
		public:
			typedef typename Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >  MMType;
			Exporter(SEXP x) : d_nrow(::Rf_length(x)), d_ncol(1) {
				const int RTYPE = ::Rcpp::traits::r_sexptype_traits<T>::rtype ;
				if (TYPEOF(x) != RTYPE)
					throw std::invalid_argument("Wrong R type for mapped vector");
				typedef typename ::Rcpp::traits::storage_type<RTYPE>::type STORAGE;
				d_start         = ::Rcpp::internal::r_vector_start<RTYPE,STORAGE>(x);
				if (::Rf_isMatrix(x)) {
					int *dims = INTEGER(::Rf_getAttrib(x, R_DimSymbol));
					d_nrow = dims[0];
					d_ncol = dims[1];
				}
			}
			MMType get() {return MMType(d_start, d_nrow, d_ncol);}
		protected:
			int   d_nrow, d_ncol;
			T*    d_start;
		};

		template <typename T> 
		class Exporter<Eigen::Matrix<T, Eigen::Dynamic, 1> >
			: public IndexingExporter<Eigen::Matrix<T, Eigen::Dynamic, 1>, T> {
		public: 
			Exporter(SEXP x) : IndexingExporter<Eigen::Matrix<T, Eigen::Dynamic, 1>, T >(x){}
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
}

#endif
