#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <cmath>
#include <numbers>
#include <format>
#include <algorithm>

#include <checks.h>
#include <polynomial.h>
#include <algorithms.h>

double calc_quadratic_error(
    polynomial approximate_polynomial,
    std::vector<double> points,
    std::vector<double> values)
{
    if (points.size() != values.size())
    {
        throw std::invalid_argument("Point and values count are not equal");
    }
    
    double error = 0;
    
    for (size_t i = 0; i < points.size(); ++i)
    {
        error += std::pow(approximate_polynomial.valueAt(points[i]) - values[i], 2);
    }
    
    return error;
}


int main()
{
    std::freopen("input.txt", "r", stdin);
    
    size_t n;
    std::cin >> n;
    
    std::vector<double> points(n), values(n);
    
    for (size_t i = 0; i < n; ++i)
    {
        std::cin >> points[i] >> values[i];
    }
    
    try
    {
        
        polynomial approximate_polynomial =
                algorithms::approximate_with_least_squares_method(points, values, 1);
        
        std::cout << "P1(x) = ";
        approximate_polynomial.println(std::cout, 4);
        std::cout << "Error = " << calc_quadratic_error(approximate_polynomial, points, values) <<
                std::endl << std::endl;
    }
    catch (std::invalid_argument const &e)
    {
        std::cout << e.what() << std::endl;
    }
    
    try
    {
        
        polynomial approximate_polynomial =
                algorithms::approximate_with_least_squares_method(points, values, 2);
        
        std::cout << "P2(x) = ";
        approximate_polynomial.println(std::cout, 4);
        std::cout << "Error = " << calc_quadratic_error(approximate_polynomial, points, values) <<
                std::endl << std::endl;
    }
    catch (std::invalid_argument const &e)
    {
        std::cout << e.what() << std::endl;
    }
    
    
    
}