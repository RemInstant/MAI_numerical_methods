#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>

#include <vecN.h>
#include <matrixNxN.h>

std::pair<matrixNxN, vecN>
jacobi_transform(
    matrixNxN const &matrix,
    vecN const &const_terms)
{
    matrixNxN tr_matrix(matrix.size());
    vecN tr_const_terms(matrix.size());
    
    for (size_t i = 0; i < matrix.size(); ++i)
    {
        tr_const_terms[i] = const_terms[i] / matrix[i][i];
        
        for (size_t j = 0; j < matrix.size(); ++j)
        {
            if (i == j)
            {
                tr_matrix[i][j] = 0;
            }
            else
            {
                tr_matrix[i][j] = -matrix[i][j] / matrix[i][i]; 
            }
        }
    }
    
    return std::make_pair(tr_matrix, tr_const_terms);
}

void validate_iterations(
    matrixNxN const &tr_matrix)
{
    double tr_matrix_l1_norm = tr_matrix.l1_norm();
    double tr_matrix_c_norm = tr_matrix.continuous_norm();
    double eps = std::numeric_limits<double>::epsilon();
    
    if (((tr_matrix_c_norm + 1 <= -eps) || (tr_matrix_c_norm - 1 >= -eps)) &&
            ((tr_matrix_l1_norm + 1 <= -eps) || (tr_matrix_l1_norm - 1 >= -eps)))
    {
        throw std::invalid_argument("The method does not converges for the system");
    }
}

std::pair<vecN, size_t>
solve_with_iterations(
    matrixNxN const &matrix,
    vecN const &const_terms,
    double eps)
{
    if (eps < 0)
    {
        throw std::invalid_argument("Invalid epsilon");
    }
    
    auto [tr_matrix, tr_const_terms] = jacobi_transform(matrix, const_terms);
    
    validate_iterations(tr_matrix);
    
    // Iterations
    vecN prev_roots(0);
    vecN roots = tr_const_terms;
    
    size_t k = 0;
    double tr_matrix_c_norm = tr_matrix.continuous_norm();
    double eps_k_coef = tr_matrix_c_norm / (1 - tr_matrix_c_norm);
    double eps_k = 2*eps;
    
    while (eps_k > eps)
    {
        ++k;
        prev_roots = std::move(roots);
        roots = tr_const_terms + tr_matrix * prev_roots;   
        
        eps_k = eps_k_coef * (roots - prev_roots).continuous_norm();
    }   
    
    return std::make_pair(roots, k);
}

std::pair<vecN, size_t>
solve_with_seidel(
    matrixNxN const &matrix,
    vecN const &const_terms,
    double eps)
{
    if (eps < 0)
    {
        throw std::invalid_argument("Invalid epsilon");
    }
    
    auto [tr_matrix, tr_const_terms] = jacobi_transform(matrix, const_terms);
    
    matrixNxN b(matrix.size(), 0);
    matrixNxN c(matrix.size(), 0);
    
    for (size_t i = 0; i < matrix.size(); ++i)
    {
        for (size_t j = 0; j < matrix.size(); ++j)
        {
            if (i < j)
            {
                b[i][j] = tr_matrix[i][j];
            }
            else
            {
                c[i][j] = tr_matrix[i][j];
            }
        }
    }
    
    matrixNxN tmp = (matrixNxN::identical(matrix.size()) - b).inversed();
    matrixNxN modified_matrix = tmp * c;
    vecN modified_const_terms = tmp * tr_const_terms;
    validate_iterations(modified_matrix);
    
    // Iterations
    vecN prev_roots(0);
    vecN roots = tr_const_terms;
    
    size_t k = 0;
    double eps_k_coef = c.continuous_norm() / (1 - modified_matrix.continuous_norm());
    double eps_k = 2*eps;
    
    // while (eps_k > eps)
    // {
    //     ++k;
        
    //     prev_roots = std::move(roots);
    //     roots = modified_const_terms + modified_matrix * prev_roots;   
        
    //     eps_k = eps_k_coef * (roots - prev_roots).l1_norm();
    // }
    
    while (eps_k > eps)
    {
        ++k;
        
        prev_roots = roots;
        
        for (size_t i = 0; i < roots.size(); ++i)
        {
            roots[i] = tr_const_terms[i];
            
            for (size_t j = 0; j < i; ++j)
            {
                roots[i] += tr_matrix[i][j] * roots[j];
            }
            for (size_t j = i; j < roots.size(); ++j)
            {
                roots[i] += tr_matrix[i][j] * prev_roots[j];
            }
        }
        
        eps_k = eps_k_coef * (roots - prev_roots).continuous_norm();
    } 
    
    return std::make_pair(roots, k);
}

int main()
{
    std::freopen("input.txt", "r", stdin);
       
    size_t n = 0;
    double eps;
    
    std::cin >> n >> eps;
    
    matrixNxN coefs(n);
    vecN const_terms(n);
    double max_row_abs_sum = 0;
    
    for (size_t i = 0; i < n; ++i)
    {
        double tmp_abs_sum = 0;
        
        for (size_t j = 0; j < n; ++j)
        {
            std::cin >> coefs[i][j];
            tmp_abs_sum += std::abs(coefs[i][j]);
        }
        
        max_row_abs_sum = std::max(max_row_abs_sum, tmp_abs_sum);
        
        std::cin >> const_terms[i];
    }
    
    std::cout << std::setprecision(15) << std::endl;
    
    try
    {
        auto [roots1, k1] = solve_with_iterations(coefs, const_terms, eps);
        
        std::cout << "Iterations count: " << k1 << std::endl;
        
        for (size_t i = 0; i < n; ++i)
        {
            std::cout << "x" << (i+1) << " = " << roots1[i] << std::endl;
        }
        
        if (vecN(n, 0).equals(coefs * roots1 - const_terms, eps * max_row_abs_sum))
        {
            std::cout << "Solution verified" << std::endl << std::endl;
        }
        else
        {
            std::cout << "Solution is wrong" << std::endl << std::endl;
        }
    }
    catch (std::invalid_argument const &e)
    {
        std::cout << e.what() << std::endl;
    }
    
    try
    {
        auto [roots2, k2] = solve_with_seidel(coefs, const_terms, eps);
        
        std::cout << "Iterations count: " << k2 << std::endl;
        
        for (size_t i = 0; i < n; ++i)
        {
            std::cout << "x" << (i+1) << " = " << roots2[i] << std::endl;
        }
        
        if (vecN(n, 0).equals(coefs * roots2 - const_terms, eps * max_row_abs_sum))
        {
            std::cout << "Solution verified" << std::endl << std::endl;
        }
        else
        {
            std::cout << "Solution is wrong" << std::endl << std::endl;
        }
    }
    catch (std::invalid_argument const &e)
    {
        std::cout << e.what() << std::endl;
    }
    
    try
    {
        std::vector<size_t> iter_cnt;
        std::vector<size_t> seidel_cnt;
        
        eps = 1;
        
        for (size_t i = 0; i < 16; ++i)
        {
            iter_cnt.push_back(solve_with_iterations(coefs, const_terms, eps).second);
            seidel_cnt.push_back(solve_with_seidel(coefs, const_terms, eps).second);
            eps /= 10;
        }
        
        for (size_t i = 0; i < iter_cnt.size(); ++i)
        {
            std::cout << std::setw(5) << ("e-" + std::to_string(i)) << ' ';
        }
        std::cout << std::endl;
        
        for (size_t i = 0; i < iter_cnt.size(); ++i)
        {
            std::cout << std::setw(5) << iter_cnt[i] << ' ';
        }
        std::cout << std::endl;
        
        for (size_t i = 0; i < iter_cnt.size(); ++i)
        {
            std::cout << std::setw(5) << seidel_cnt[i] << ' ';
        }
        std::cout << std::endl;
    }
    catch (std::invalid_argument const &e)
    {
        std::cout << e.what() << std::endl;
    }
}