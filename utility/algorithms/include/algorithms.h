#ifndef _NUM_METHODS_UTIL_ALGORITHMS_H_
#define _NUM_METHODS_UTIL_ALGORITHMS_H_

#include <checks.h>
#include <vecN.h>
#include <matrixNxN.h>
#include <polynome.h>

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
    
    polynome interpolate_with_lagrange(
        std::vector<double> points,
        std::vector<double> values);
    
    polynome interpolate_with_newton(
        std::vector<double> points,
        std::vector<double> values);
}

#endif // _NUM_METHODS_UTIL_ALGORITHMS_H_