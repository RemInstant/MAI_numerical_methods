#include <stdexcept>
#include <limits>
#include <algorithm>

#include "../include/algorithms.h"

std::vector<vecN> solve_tridiagonal_linear_equation_system(
    matrixNxN const &coefs,
    std::vector<vecN> const &constant_terms)
{
    size_t equation_dimension = coefs.size();
    size_t term_group_count = constant_terms.size();
    
    vecN coefs_a(equation_dimension);
    vecN coefs_b(equation_dimension);
    vecN coefs_c(equation_dimension);
    std::vector<vecN> const &coefs_d = constant_terms;
    
    for (size_t i = 0; i < equation_dimension; ++i)
    {
        coefs_a[i] = i > 0 ? coefs[i][i-1] : 0;
        coefs_b[i] = coefs[i][i];
        coefs_c[i] = i + 1 < equation_dimension ? coefs[i][i+1] : 0;
    }
    
    // First phase of Thomas algorithm
    vecN coefs_p(equation_dimension, 0);
    std::vector<vecN> coefs_q(term_group_count, vecN(equation_dimension, 0));
    
    if (abs(coefs_b[0]) < std::numeric_limits<double>::epsilon())
    {
        throw std::invalid_argument("Incorrect system");
    }
    
    coefs_p[0] = -coefs_c[0] / coefs_b[0];
    coefs_q[0] = coefs_d[0] / coefs_b[0];
    
    for (size_t i = 1; i < equation_dimension; ++i)
    {
        double div = coefs_b[i] + coefs_a[i] * coefs_p[i-1];
        
        if (abs(div) < std::numeric_limits<double>::epsilon())
        {
            throw std::invalid_argument("Incorrect system");
        }
        
        coefs_p[i] = -coefs_c[i] / div;
        
        for (size_t k = 0; k < term_group_count; ++k)
        {
            coefs_q[k][i] = (coefs_d[k][i] - coefs_a[i] * coefs_q[k][i-1]) / div;
        }
    }
    
    // Second phase of Thomas algorithm
    std::vector<vecN> roots(term_group_count, vecN(equation_dimension, 0));
    
    for (size_t k = 0; k < term_group_count; ++k)
    {
        roots[k][equation_dimension - 1] = coefs_q[k][equation_dimension - 1];
        
        for (size_t i = equation_dimension - 1; i > 0; --i)
        {
            roots[k][i-1] = coefs_p[i-1] * roots[k][i] + coefs_q[k][i-1];
        }
    }
    
    return roots;
}

vecN solve_tridiagonal_linear_equation_system(
    matrixNxN coefs,
    vecN constant_terms)
{
    return solve_tridiagonal_linear_equation_system(coefs, std::vector<vecN>(1, constant_terms))[0];
}

std::vector<vecN> algorithms::solve_linear_equation_system(
    matrixNxN coefs,
    std::vector<vecN> constant_terms)
{
    if (constant_terms.size() == 0)
    {
        return std::vector<vecN>(0);
    }
    
    if (coefs.size() != constant_terms[0].size())
    {
        throw std::invalid_argument("Dimensions does not correspond");
    }
    
    size_t equation_dimension = coefs.size();
    size_t term_group_count = constant_terms.size();
    std::vector<vecN> roots(term_group_count, vecN(coefs.size()));
    
    // First phase of gauss elimination
    for (size_t i = 0; i < equation_dimension; ++i)
    {
        // Search of leading element
        size_t idx_to_swap = i;
        
        for (size_t j = idx_to_swap + 1; j < equation_dimension; ++j)
        {
            if (std::abs(coefs[j][i]) > std::abs(coefs[idx_to_swap][i]))
            {
                idx_to_swap = j;
            }
        }
        
        std::swap(coefs[i], coefs[idx_to_swap]);
        for (size_t k = 0; k < term_group_count; ++k)
        {
            std::swap(constant_terms[k][i], constant_terms[k][idx_to_swap]);
        }
        
        if (std::abs(coefs[i][i]) < std::numeric_limits<double>::epsilon())
        {
            throw std::invalid_argument("Equation cannot be solved");
        }
        
        for (size_t j = i+1; j < equation_dimension; ++j)
        {
            if (coefs[j][i] < std::numeric_limits<double>::epsilon())
            {
                continue;
            }
            
            double mult = coefs[j][i] / coefs[i][i];
            
            for (size_t k = i; k < equation_dimension; ++k)
            {
                coefs[j][k] -= coefs[i][k] * mult;
            }
            for (size_t k = 0; k < term_group_count; ++k)
            {
                constant_terms[k][j] -= constant_terms[k][i] * mult;
            }
        }
    }
    
    // Second phase of gauss elimination
    for (size_t k = 0; k < term_group_count; ++k)
    {
        for (size_t i = equation_dimension; i > 0; --i)
        {
            roots[k][i-1] = constant_terms[k][i-1];
            
            for (size_t j = i; j < equation_dimension; ++j)
            {
                roots[k][i-1] -= coefs[i-1][j] * roots[k][j];
            }
            
            roots[k][i-1] /= coefs[i-1][i-1];
        }
        
    }
    
    return roots;
}

vecN algorithms::solve_linear_equation_system(
    matrixNxN coefs,
    vecN constant_terms)
{
    return solve_linear_equation_system(coefs, std::vector<vecN>(1, constant_terms))[0];
}

polynomial algorithms::interpolate_with_lagrange(
    std::vector<double> const &points,
    std::vector<double> const &values)
{
    if (points.size() != values.size())
    {
        throw std::invalid_argument("Point and values count are not equal");
    }
    if (points.size() == 0)
    {
        return polynomial();
    }
    
    polynomial interpolation;
    
    for (size_t i = 0; i < points.size(); ++i)
    {
        polynomial li({1});
        
        for (size_t j = 0; j < points.size(); ++j)
        {
            if (i == j)
            {
                continue;
            }
            
            li *= polynomial({-points[j], 1});
            li /= points[i] - points[j];
        }
        
        interpolation += li * values[i];
    }
    
    return interpolation;
}

polynomial algorithms::interpolate_with_newton(
    std::vector<double> const &points,
    std::vector<double> const &values)
{
    if (points.size() != values.size())
    {
        throw std::invalid_argument("Point and values count are not equal");
    }
    if (points.size() == 0)
    {
        return polynomial();
    }
    
    std::vector<std::vector<double>> div_diff(points.size(), std::vector<double>(points.size()));
       
    for (size_t i = 0; i < points.size(); ++i)
    {
        div_diff[i][0] = values[i]; 
    }
    
    for (size_t j = 1; j < points.size(); ++j)
    {
        for (size_t i = 0; i+j < points.size(); ++i)
        {
            div_diff[i][j] = (div_diff[i+1][j-1] - div_diff[i][j-1]) / (points[i+j] - points[i]);
        }
    }
    
    polynomial interpolation, mult({1});
    
    for (size_t i = 0; i < points.size(); ++i)
    {
        interpolation += mult * div_diff[0][i];
        mult *= polynomial({-points[i], 1});
    }
    
    return interpolation;
}

piecewise_polynomial algorithms::build_spline(
    std::vector<double> const &points,
    std::vector<double> const &values)
{
    if (points.size() != values.size())
    {
        throw std::invalid_argument("Point and values count are not equal");
    }
    if (points.size() == 0)
    {
        return piecewise_polynomial();
    }
    
    std::vector<std::pair<double, double>> pv(points.size());
    for (size_t i = 0; i < pv.size(); ++i)
    {
        pv[i].first = points[i];
        pv[i].second = values[i];
    }
    
    std::sort(pv.begin(), pv.end(),
            [](std::pair<double, double> const &a, std::pair<double, double> const &b)
            {
                return a.first < b.first;
            });
    
    std::vector<double> h(points.size() - 1);
    matrixNxN eq_coefs(points.size() - 2);
    vecN eq_const_terms(points.size() - 2);
    
    for (size_t i = 0; i < h.size(); ++i)
    {
        h[i] = pv[i+1].first - pv[i].first;
    }
    
    for (size_t i = 0; i < eq_coefs.size(); ++i)
    {
        if (i > 0)
        {
            eq_coefs[i][i-1] = h[i];
        }
        
        eq_coefs[i][i] = 2 * (h[i] + h[i+1]);
        
        if (i + 1 < eq_coefs.size())
        {
            eq_coefs[i][i+1] = h[i+1];
        }
        
        eq_const_terms[i] = 3 * ((pv[i+2].second - pv[i+1].second) / h[i+1] -
                (pv[i+1].second - pv[i].second) / h[i]);
    }
    
    vecN roots = solve_tridiagonal_linear_equation_system(eq_coefs, eq_const_terms);
    
    std::vector<double> c(points.size() - 1, 0);
    for (size_t i = 0; i < roots.size(); ++i)
    {
        c[i+1] = roots[i];
    }
    
    piecewise_polynomial spline;
    for (size_t i = 0; i < c.size(); ++i)
    {
        double a = pv[i].second;
        double b = (pv[i+1].second - pv[i].second) / h[i] -
                (i + 1 < c.size()
                ? h[i] * (c[i+1] + 2*c[i]) / 3
                : 2 * h[i] * c[i] / 3);
        double d = ((i + 1 < c.size() ? c[i+1] : 0) - c[i]) / (3 * h[i]);
        
        polynomial p({-pv[i].first, 1});
        polynomial q = polynomial({a}) + p*b + p*p*c[i] + p*p*p*d;
        
        if (i == 0)
        {
            spline = piecewise_polynomial(pv[i].first, pv[i+1].first, q);
        }
        else
        {
            spline.push_back(pv[i+1].first, q);
        }
    }
    
    return spline;
}