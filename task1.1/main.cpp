#include <iostream>
#include <iomanip>
#include <vector>

void print_augmented_matrix(
        std::vector<std::vector<double>> &augmented_matrix)
{
    for (size_t i = 0; i < augmented_matrix.size(); ++i)
    {
        for (size_t j = 0; j < augmented_matrix[0].size(); ++j)
        {
            std::cout << augmented_matrix[i][j] << ' ';
        }
        std::cout << std::endl;
    }
}

size_t search_leading_element(
    std::vector<std::vector<double>> const &augmented_matrix,
    size_t column)
{
    size_t n = augmented_matrix.size();
    size_t idx = column;
    
    for (size_t i = idx + 1; i < n; ++i)
    {
        if (std::abs(augmented_matrix[i][column]) > std::abs(augmented_matrix[idx][column]))
        {
            idx = i;
        }
    }
    
    return idx;
}

void process_gaussian_elimination_first_phase(
        std::vector<std::vector<double>> &augmented_matrix,
        double eps)
{
    size_t n = augmented_matrix.size();
    size_t m = augmented_matrix[0].size();
    size_t swap_count = 0;
    
    for (size_t i = 0; i < n; ++i)
    {
        size_t idx_to_swap = search_leading_element(augmented_matrix, i);
        swap_count += (idx_to_swap > i) ? (idx_to_swap - i) : (i - idx_to_swap);
        
        std::swap(augmented_matrix[i], augmented_matrix[idx_to_swap]);
        
        if (std::abs(augmented_matrix[i][i]) < eps)
        {
            throw std::invalid_argument("Invalid matrix");
        }
        
        for (size_t j = i+1; j < n; ++j)
        {
            double mult = augmented_matrix[j][i] / augmented_matrix[i][i];
            
            for (size_t k = i; k < m; ++k)
            {
                augmented_matrix[j][k] -= augmented_matrix[i][k] * mult;
            }
        }
    }
    
    if (swap_count & 1)
    {
        for (size_t i = 0; i < n+1; ++i)
        {
            augmented_matrix[0][i] *= -1;
        }
    }
}

std::vector<std::vector<double>>
process_gaussian_elimination_second_phase(
        std::vector<std::vector<double>> &augmented_matrix)
{
    size_t n = augmented_matrix.size();
    size_t m = augmented_matrix[0].size();
    
    std::vector<std::vector<double>> roots(n, std::vector<double>(m - n));
    
    for (size_t k = 0; k < m - n; ++k)
    {
        for (size_t i = n; i > 0; --i)
        {
            roots[i-1][k] = augmented_matrix[i-1][n+k];
            
            for (size_t j = i; j < n; ++j)
            {
                roots[i-1][k] -= augmented_matrix[i-1][j] * roots[j][k];
            }
            
            roots[i-1][k] /= augmented_matrix[i-1][i-1];
        }
        
    }
    
    return roots;
}

bool validate_solution(
        std::vector<std::vector<double>> const &augmented_matrix,
        std::vector<std::vector<double>> const &roots,
        double eps)
{
    size_t n = augmented_matrix.size();
    size_t m = augmented_matrix[0].size();
    
    for (size_t k = 0; k < m - n; ++k)
    {
        for (size_t i = 0; i < n; ++i)
        {
            double sum = 0;
            
            for (size_t j = 0; j < n; ++j)
            {
                sum += augmented_matrix[i][j] * roots[j][k];
            }
            
            if (std::abs(sum - augmented_matrix[i][n+k]) > eps)
            {
                return false;
            }
        }
    }
    
    return true;
}

double calc_determinant(
        std::vector<std::vector<double>> &augmented_matrix)
{
    size_t n = augmented_matrix.size();
    double det = 1;
    
    for (size_t i = 0; i < n; ++i)
    {
        det *= augmented_matrix[i][i];
    }
    
    return det;
}

std::vector<std::vector<double>>
inverse_matrix(
    std::vector<std::vector<double>> matrix,
    double eps)
{
    size_t n = matrix.size();
    
    std::vector<std::vector<double>> inv_matrix(n, std::vector<double>(n));
    std::vector<std::vector<double>> augmented_matrix(n, std::vector<double>(2*n));
    
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            augmented_matrix[i][j] = matrix[i][j];
            augmented_matrix[i][n+j] = (i == j) ? 1 : 0;
        }
    }
    
    process_gaussian_elimination_first_phase(augmented_matrix, eps);
    auto roots = process_gaussian_elimination_second_phase(augmented_matrix);
    
    return roots;
}

int main()
{
    std::freopen("input.txt", "r", stdin);
    std::cout << std::fixed << std::setprecision(6);
    
    double const EPS = 1e-8;
    
    size_t n = 0;
    std::cin >> n;
    
    std::vector<std::vector<double>> matrix(n, std::vector<double>(n));
    std::vector<std::vector<double>> augmented_matrix(n, std::vector<double>(n+1));
    
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            std::cin >> matrix[i][j];
            augmented_matrix[i][j] = matrix[i][j];
        }
        
        std::cin >> augmented_matrix[i][n];
    }
    
    try
    {
        process_gaussian_elimination_first_phase(augmented_matrix, EPS);
        auto roots = process_gaussian_elimination_second_phase(augmented_matrix);
        
        for (size_t i = 0; i < n; ++i)
        {
            std::cout << "x" << (i+1) << " = " << roots[i][0] << std::endl;
        }
        
        std::cout << std::endl;
        
        if (validate_solution(augmented_matrix, roots, EPS))
        {
            std::cout << "Solution is verified" << std::endl << std::endl;
        }
        else
        {
            std::cout << "Solution is incorrect" << std::endl << std::endl;
        }
    }
    catch (std::invalid_argument const &ex)
    {
        std::cout << "Inconsistent equation system" << std::endl;
        return 1;
    }
    
    double determinant = calc_determinant(augmented_matrix);
    std::cout << "Determinant = " << determinant << std::endl << std::endl;
    
    try
    {
        auto inv_matrix = inverse_matrix(matrix, EPS);
        
        std::cout << "Inverse matrix: " << std::endl;
        for (auto line : inv_matrix)
        {
            for (auto elem : line)
            {
                std::cout << elem << ' ';
            }
            
            std::cout << std::endl;
        }
        
        std::cout << std::endl;
        // verifying inverse matrix
        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = 0; j < n; ++j)
            {
                double prod_elem = 0;
                
                for (size_t k = 0; k < n; ++k)
                {
                    prod_elem += matrix[i][k] * inv_matrix[k][j];
                }
                
                std::cout << std::setprecision(20) << prod_elem << ' ';
                
                if ((i == j && std::abs(prod_elem - 1) > EPS) ||
                        (i != j && std::abs(prod_elem) > EPS))
                {
                    std::cout << "The inverse matrix is incorrect!" << std::endl;
                    return 1;
                }
            }
            
            std::cout << std::endl;
        }
        
        std::cout << std::endl << "The inverse matrix is verified" << std::endl;
    }
    catch (std::invalid_argument const &ex)
    {
        std::cout << "The matrix has no inverse" << std::endl;
        return 1;
    }
}