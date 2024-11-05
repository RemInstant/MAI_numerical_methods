#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <cmath>
#include <numbers>
#include <format>
#include <algorithm>

#include <checks.h>
#include <vecN.h>
#include <matrixNxN.h>
#include <algorithms.h>

int get_interpolation_accuracy(
    polynomial interpolation,
    std::vector<double> points,
    std::vector<double> values)
{
    if (points.size() != values.size())
    {
        throw std::invalid_argument("Point and values count are not equal");
    }
    
    double eps = 1.0;
    size_t accuracy = 0;
    
    for (; eps > std::numeric_limits<double>::epsilon(); ++accuracy, eps /= 10)
    {
        for (size_t i = 0; i < points.size(); ++i)
        {
            if (std::abs(interpolation.valueAt(points[i]) - values[i]) > eps)
            {
                return accuracy - 1;
            }
        }
    }
    
    return accuracy;
}


int main()
{
    std::freopen("input.txt", "r", stdin);
    double test_point = 3.0 / 16 * std::numbers::pi;
    double test_value = std::tan(test_point);
    
    std::cout << std::fixed << std::setw(8);
    
    try
    {
        std::cout << "Lagrange interpolation" << std::endl;
        
        std::vector<double> points =
        {
            0,
            1.0 / 8 * std::numbers::pi,
            2.0 / 8 * std::numbers::pi,
            3.0 / 8 * std::numbers::pi,
        };
    
        std::vector<double> values(points);
        std::for_each(values.begin(), values.end(), [](double &x){ x = std::tan(x); });
        
        polynomial interpolation = algorithms::interpolate_with_lagrange(points, values);
        std::cout << "P(x) = ";
        interpolation.println();
        
        int accuracy = get_interpolation_accuracy(interpolation, points, values);
        std::cout << "Interpolation is " << (accuracy != -1
                ? std::format("verified (accurate to {} decimal places)", accuracy)
                : "wrong") << std::endl;
        
        double approx_test_value = interpolation.valueAt(test_point);
        
        std::cout << "tan(x*) = " << test_value << std::endl;
        std::cout << "P(x*)   = " << approx_test_value << std::endl;
        std::cout << "error   = " << std::abs(test_value - approx_test_value) << std::endl;
        std::cout << std::endl;
    }
    catch (std::invalid_argument e)
    {
        std::cout << e.what() << std::endl;
    }
    
    try
    {
        std::cout << "Newton interpolation" << std::endl;
        
        std::vector<double> points =
        {
            0,
            1.0 / 8 * std::numbers::pi,
            1.0 / 3 * std::numbers::pi,
            3.0 / 8 * std::numbers::pi,
        };
    
        std::vector<double> values(points);
        std::for_each(values.begin(), values.end(), [](double &x){ x = std::tan(x); });
        
        polynomial interpolation = algorithms::interpolate_with_newton(points, values);
        std::cout << "P(x) = ";
        interpolation.println();
        
        int accuracy = get_interpolation_accuracy(interpolation, points, values);
        std::cout << "Interpolation is " << (accuracy != -1
                ? std::format("verified (accurate to {} decimal places)", accuracy)
                : "wrong") << std::endl;
        
        double approx_test_value = interpolation.valueAt(test_point);
        
        std::cout << "tan(x*) = " << test_value << std::endl;
        std::cout << "P(x*)   = " << approx_test_value << std::endl;
        std::cout << "error   = " << std::abs(test_value - approx_test_value) << std::endl;
        std::cout << std::endl;
    }
    catch (std::invalid_argument e)
    {
        std::cout << e.what() << std::endl;
    }
}