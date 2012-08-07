#include "rcppeigen_hello_world.h"

using namespace Rcpp ;

SEXP rcppeigen_hello_world(){
	
    Eigen::MatrixXd m1 = Eigen::MatrixXd::Identity( 3, 3 ) ;
    Eigen::MatrixXd m2 = Eigen::MatrixXd::Random( 3, 3 ) ;
    
    return List::create(Named("m1")   = m1,
			Named("m2")   = m2,
			Named("comb") = m1 + 3 *(m1 + m2));
}

