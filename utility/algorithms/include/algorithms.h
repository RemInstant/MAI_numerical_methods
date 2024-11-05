#ifndef _NUM_METHODS_UTIL_ALGORITHMS_H_
#define _NUM_METHODS_UTIL_ALGORITHMS_H_

#include <checks.h>
#include <vecN.h>
#include <matrixNxN.h>
#include <polynomial.h>
#include <piecewise_polynomial.h>

namespace algorithms
{
    // TODO: rename to system
    std::vector<vecN> solve_linear_equation(
        matrixNxN coefs,
        std::vector<vecN> constant_terms);
        
    vecN solve_linear_equation(
        matrixNxN coefs,
        vecN constant_terms);
    
    polynomial interpolate_with_lagrange(
        std::vector<double> const &points,
        std::vector<double> const &values);
    
    polynomial interpolate_with_newton(
        std::vector<double> const &points,
        std::vector<double> const &values);
    
    piecewise_polynomial build_spline(
        std::vector<double> const &points,
        std::vector<double> const &values);
}

#endif // _NUM_METHODS_UTIL_ALGORITHMS_H_