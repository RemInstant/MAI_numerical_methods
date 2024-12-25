#ifndef _NUM_METHODS_UTIL_ALGORITHMS_H_
#define _NUM_METHODS_UTIL_ALGORITHMS_H_

#include <checks.h>
#include <vecN.h>
#include <matrixNxN.h>
#include <polynomial.h>
#include <piecewise_polynomial.h>

namespace algorithms
{
    std::vector<vecN> solve_linear_equation_system(
        matrixNxN coefs,
        std::vector<vecN> constant_terms);
        
    vecN solve_linear_equation_system(
        matrixNxN coefs,
        vecN constant_terms);
    
    polynomial interpolate_with_lagrange(
        std::vector<double> const &points,
        std::vector<double> const &values,
        bool trace = false,
        size_t trace_precision = 3);
    
    polynomial interpolate_with_newton(
        std::vector<double> const &points,
        std::vector<double> const &values,
        bool trace = false,
        size_t trace_precision = 3);
    
    piecewise_polynomial build_spline(
        std::vector<double> const &points,
        std::vector<double> const &values);
    
    polynomial approximate_with_least_squares_method(
        std::vector<double> const &points,
        std::vector<double> const &values,
        size_t degree);
}

#endif // _NUM_METHODS_UTIL_ALGORITHMS_H_