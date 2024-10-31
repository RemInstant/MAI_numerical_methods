#ifndef _NUM_METHODS_UTIL_ALGORITHMS_
#define _NUM_METHODS_UTIL_ALGORITHMS_

#include <checks.h>
#include <vecN.h>
#include <matrixNxN.h>

namespace algorithms
{
    std::vector<vecN> solve_linear_equation(
        matrixNxN coefs,
        std::vector<vecN> constant_terms,
        double eps);
        
    vecN solve_linear_equation(
        matrixNxN coefs,
        vecN constant_terms,
        double eps);
}

#endif // _NUM_METHODS_UTIL_ALGORITHMS_