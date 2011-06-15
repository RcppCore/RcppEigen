#include "rcppeigen_hello_world.h"

using namespace Rcpp ;

SEXP rcppeigen_hello_world(){
	
    Eigen::MatrixXd m1 = Eigen::MatrixXd::Identity( 3, 3 ) ;
    Eigen::MatrixXd m2 = Eigen::MatrixXd::Random( 3, 3 ) ;
				// need more flexible wrap methods
    Eigen::MatrixXd comb = m1 + 3 * (m1 + m2);
    
    return List::create(_["m1"]   = m1,
			_["m2"]   = m2,
			_["comb"] = comb);
}

