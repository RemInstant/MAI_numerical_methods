#include <stdexcept>
#include <limits>

#include "../include/algorithms.h"

std::vector<vecN> algorithms::solve_linear_equation(
    matrixNxN coefs,
    std::vector<vecN> constant_terms,
    double eps)
{
    checks::throw_if_invalid_eps(eps);
    
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

vecN algorithms::solve_linear_equation(
    matrixNxN coefs,
    vecN constant_terms,
    double eps)
{
    return solve_linear_equation(coefs, std::vector<vecN>(1, constant_terms), eps)[0];
}