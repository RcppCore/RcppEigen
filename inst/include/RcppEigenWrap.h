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
          
        // helper trait to identify if the object has dense storage
        template <typename T> struct is_dense : Rcpp::traits::same_type<typename T::StorageKind,Eigen::Dense>{} ;
        
        // for plain dense objects
        template <typename T> 
        SEXP eigen_wrap_plain_dense( const T& obj, Rcpp::traits::true_type ){
            // FIXME: deal with RowMajor, etc ...
            const int RTYPE = Rcpp::traits::r_sexptype_traits<typename T::Scalar>::rtype ;
            if( T::ColsAtCompileTime == 1 ) {
                return wrap( obj.data(), obj.data() + obj.size() ) ;
            } else {
                Rcpp::Matrix<RTYPE> x( obj.rows(), obj.cols(), obj.data() ) ;
                return x; 
            }   
        }
        
        // for plain sparse objects
        template <typename T> 
        SEXP eigen_wrap_plain_dense( const T& object, Rcpp::traits::false_type ){
            typedef typename T::Scalar Scalar ;
            const int RTYPE = Rcpp::traits::r_sexptype_traits<Scalar>::rtype  ;  
            int          nnz = object.nonZeros(), p = object.outerSize();
	        Dimension    dim(object.innerSize(), p);
	        const int    *ip = object._innerIndexPtr(), *pp = object._outerIndexPtr();
	        const Scalar      *xp = object._valuePtr();
	        IntegerVector iv(ip, ip + nnz), pv(pp, pp + p + 1);
	        Vector<RTYPE> xv(xp, xp + nnz);
	        
	        return List::create(_["Dim"] = dim,
	        								 _["i"]   = iv,
	        								 _["p"]   = pv,
	        								 _["x"]   = xv);
	    } 
        
        
        // plain object, so we can assume data() and size()
        template <typename T>
        inline SEXP eigen_wrap_is_plain( const T& obj, ::Rcpp::traits::true_type ){
            return eigen_wrap_plain_dense( obj, typename is_dense<T>::type() ) ;
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

		template<typename T>
		class Exporter<Eigen::MappedSparseMatrix<T> > {
		public:
			Exporter(SEXP x)
				: d_x(x), d_dims(d_x.slot("Dim")), d_i(d_x.slot("i")), d_p(d_x.slot("p")) {
				if (!d_x.is("CsparseMatrix")) 
					throw std::invalid_argument("Need S4 class CsparseMatrix for an mapped sparse matrix");
				const int RTYPE = ::Rcpp::traits::r_sexptype_traits<T>::rtype ;
				SEXP xx = d_x.slot("x");
				if (TYPEOF(xx) != RTYPE)
					throw std::invalid_argument("Wrong R type for mapped sparse matrix");
				typedef typename ::Rcpp::traits::storage_type<RTYPE>::type STORAGE;
				d_start         = ::Rcpp::internal::r_vector_start<RTYPE,STORAGE>(xx);
			}
			Eigen::MappedSparseMatrix<T> get() {
				return Eigen::MappedSparseMatrix<T>(d_dims[0], d_dims[1], d_p[d_dims[1]],
													d_p.begin(), d_i.begin(), d_start);
			}
		protected:
			S4            d_x;
			T*            d_start;
			IntegerVector d_dims, d_i, d_p;
		};

		template<typename T>
		class Exporter<Eigen::SparseMatrix<T> > {
		public:
			Exporter(SEXP x)
				: d_x(x), d_dims(d_x.slot("Dim")), d_i(d_x.slot("i")), d_p(d_x.slot("p")) {
				if (!d_x.is("CsparseMatrix"))
					throw std::invalid_argument("Need S4 class CsparseMatrix for an mapped sparse matrix");
				const int RTYPE = ::Rcpp::traits::r_sexptype_traits<T>::rtype ;
				SEXP xx = d_x.slot("x");
				if (TYPEOF(xx) != RTYPE) // should coerce instead - see Rcpp/inst/include/Rcpp/internal/export.h
					throw std::invalid_argument("Wrong R type for sparse matrix");
				typedef typename ::Rcpp::traits::storage_type<RTYPE>::type STORAGE;
				d_start         = ::Rcpp::internal::r_vector_start<RTYPE,STORAGE>(xx);
			}
			Eigen::SparseMatrix<T> get() {
				Eigen::SparseMatrix<T>  ans(d_dims[0], d_dims[1]);
				ans.reserve(d_p[d_dims[1]]);
				for(int j = 0; j < d_dims[1]; ++j) {
					ans.startVec(j);
					for (int k = d_p[j]; k < d_p[j + 1]; ++k) ans.insertBack(d_i[k], j) = d_start[k];
				}
				ans.finalize();  
				return ans;
			}
		protected:
			S4            d_x;
			T*            d_start;
			IntegerVector d_dims, d_i, d_p;
		};
				
    } // namespace traits
}

#endif
