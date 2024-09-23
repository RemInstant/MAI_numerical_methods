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
            
            for (size_t k = i; k < n+1; ++k)
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

std::vector<double>
process_gaussian_elimination_second_phase(
        std::vector<std::vector<double>> &augmented_matrix)
{
    size_t n = augmented_matrix.size();
    
    std::vector<double> roots(n);
    
    for (size_t i = n; i > 0; --i)
    {
        roots[i-1] = augmented_matrix[i-1][n];
        
        for (size_t j = i; j < n; ++j)
        {
            roots[i-1] -= augmented_matrix[i-1][j] * roots[j];
        }
        
        roots[i-1] /= augmented_matrix[i-1][i-1];
    }
    
    return roots;
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
    
    for (size_t i = 0; i < n; ++i)
    {
        std::vector<std::vector<double>> augmented_matrix(n, std::vector<double>(n+1));
        for (size_t j = 0; j < n; ++j)
        {
            for (size_t k = 0; k < n; ++k)
            {
                augmented_matrix[j][k] = matrix[j][k];
            }
            
            augmented_matrix[j][n] = i == j ? 1 : 0;
        }
        
        process_gaussian_elimination_first_phase(augmented_matrix, eps);
        auto roots = process_gaussian_elimination_second_phase(augmented_matrix);
        
        for (size_t j = 0; j < n; ++j)
        {
            inv_matrix[j][i] = roots[j];
        }
    }
    
    return inv_matrix;
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
            std::cout << "x" << (i+1) << " = " << roots[i] << std::endl;
        }
        
        std::cout << std::endl;
    }
    catch (std::invalid_argument const &ex)
    {
        std::cout << "Invalid equation system" << std::endl;
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
                
                //std::cout << prod_elem << ' ';
                
                if ((i == j && std::abs(prod_elem - 1) > EPS) ||
                        (i != j && std::abs(prod_elem) > EPS))
                {
                    std::cout << "The inverse matrix is incorrect!" << std::endl;
                    return 1;
                }
            }
            
            //std::cout << std::endl;
        }
        
        std::cout << std::endl << "The inverse matrix is verified" << std::endl;
    }
    catch (std::invalid_argument const &ex)
    {
        std::cout << "The matrix has no inverse" << std::endl;
        return 1;
    }
}