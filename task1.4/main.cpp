#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <cmath>
#include <numbers>

#include <vecN.h>
#include <matrixNxN.h>

double calc_check_sum(
    matrixNxN const &matrix)
{
    double check_sum = 0;
    
    for (size_t i = 1; i < matrix.size(); ++i)
    {
        for (size_t j = 0; j < i; ++j)
        {
            check_sum += matrix[i][j] * matrix[i][j];
        }
    }
    
    return std::sqrt(check_sum);
}

std::pair< std::vector< std::pair<double, vecN> >, size_t>
calc_eigen_with_rotation_iterations(
    matrixNxN matrix,
    double eps)
{
    if (eps < 0)
    {
        throw std::invalid_argument("Invalid epsilon");
    }
    
    if (matrix.size() == 0)
    {
        throw std::invalid_argument("Empty matrix");
    }
    if (matrix.size() == 1)
    {
        return std::make_pair(std::vector<std::pair<double, vecN>>(
            1, std::make_pair(matrix[0][0], vecN(1, 1))), 0);
    }
    
    if (!matrix.is_symmetric(eps))
    {
        throw std::invalid_argument("Matrix is not symmetric");
    }
    
    matrixNxN transformation_matrix = matrixNxN::identical(matrix.size());
    size_t iter = 0;
    
    while (calc_check_sum(matrix) > eps)
    {
        size_t mx_i = 0, mx_j = 1;
        
        for (size_t i = 1; i < matrix.size(); ++i)
        {
            for (size_t j = 0; j < i; ++j)
            {
                if (std::abs(matrix[i][j]) > std::abs(matrix[mx_i][mx_j]))
                {
                    mx_i = i;
                    mx_j = j;
                }
            }
        }
        
        double phi = 0.5 * std::atan(2 * matrix[mx_i][mx_j] / (matrix[mx_i][mx_i] - matrix[mx_j][mx_j]));
        
        matrixNxN rotate_matrix = matrixNxN::identical(matrix.size());
        
        rotate_matrix[mx_i][mx_i] =  std::cos(phi);
        rotate_matrix[mx_i][mx_j] = -std::sin(phi);
        rotate_matrix[mx_j][mx_i] =  std::sin(phi);
        rotate_matrix[mx_j][mx_j] =  std::cos(phi);
    
        transformation_matrix *= rotate_matrix;
        matrix = rotate_matrix.transposed() * matrix * rotate_matrix; 
        ++iter;
    }
    
    transformation_matrix = transformation_matrix.transposed();
    
    std::vector< std::pair<double, vecN> > ans(matrix.size());
    
    for (size_t i = 0; i < matrix.size(); ++i)
    {
        ans[i] = std::make_pair(matrix[i][i], transformation_matrix[i]);
    }
    
    return std::make_pair(ans, iter);
}

std::pair< std::pair<double, vecN> , size_t>
calc_eigen_radius_with_power_iterations(
    matrixNxN const &matrix,
    double eps)
{
    if (eps < 0)
    {
        throw std::invalid_argument("Invalid epsilon");
    }
    
    if (matrix.size() == 0)
    {
        throw std::invalid_argument("Empty matrix");
    }
    if (matrix.size() == 1)
    {
        return std::make_pair(std::make_pair(matrix[0][0], vecN(1, 1)), 0);
    }
    
    size_t iter = 0;
    double x, x_prev;
    
    vecN h(matrix.size(), 1);
    h /= h.continuous_norm();
    
    do
    {
        double tmp = h[0];
        h = matrix * h;
        x_prev = x;
        x = h[0] / tmp;
        h /= h.continuous_norm();
        ++iter;
    } while (std::abs(x - x_prev) > eps);
    
    return std::make_pair(std::make_pair(x, h), iter);
}

bool validate_eigen(
    matrixNxN const &matrix,
    std::pair<double, vecN> const &eigen,
    double eps)
{
    double mx = 1;
    
    for (size_t i = 0; i < matrix.size(); ++i)
    {
        for (size_t j = 0; j < matrix.size(); ++j)
        {
            mx = std::max(mx, std::abs(matrix[i][j]));
        }
    }
    
    return (matrix * eigen.second - eigen.first * eigen.second).continuous_norm() <= eps * mx;
}



int main()
{
    std::freopen("input.txt", "r", stdin);
       
    size_t n = 0;
    double eps;
    
    std::cin >> n >> eps;
    
    matrixNxN matrix(n);
    
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            std::cin >> matrix[i][j];
        }
    }
    
    std::cout << std::setprecision(8) << std::endl;
    
    try
    {
        auto [eigen, k] = calc_eigen_with_rotation_iterations(matrix, eps);
        
        std::cout << "Rotation iterations count: " << k << std::endl;
        
        for (size_t i = 0; i < n; ++i)
        {
            std::cout << "x" << (i+1) << " = " << eigen[i].first << " : h = (";
            
            for (size_t j = 0; j < n; ++j)
            {
                std::cout << eigen[i].second[j];
                if (j + 1 < n)
                {
                    std::cout << ", ";
                }
            }
            
            std::cout << ")^T" << std::endl;
        }
        
        bool flag = true;
        
        for (size_t i = 0; i < n && flag; ++i)
        {
            flag = validate_eigen(matrix, eigen[i], eps);
        }
        
        for (size_t i = 0; i < n && flag; ++i)
        {
            for (size_t j = i + 1; j < n && flag; ++j)
            {
                flag = std::abs(vecN::scalar_prod(eigen[i].second, eigen[j].second)) <
                        std::numeric_limits<double>::epsilon();
            }
        }
        
        if (flag)
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
        auto [eigen_radius, k] = calc_eigen_radius_with_power_iterations(matrix, eps);
        
        std::cout << "Rotation iterations count: " << k << std::endl;
        std::cout << "x_rad" << " = " << eigen_radius.first << " : h = (";
        
        for (size_t i = 0; i < n; ++i)
        {
            std::cout << eigen_radius.second[i];
            if (i + 1 < n)
            {
                std::cout << ", ";
            }
        }
        
        std::cout << ")^T" << std::endl;

        
        if (validate_eigen(matrix, eigen_radius, eps))
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
            iter_cnt.push_back(calc_eigen_with_rotation_iterations(matrix, eps).second);
            seidel_cnt.push_back(calc_eigen_radius_with_power_iterations(matrix, eps).second);
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