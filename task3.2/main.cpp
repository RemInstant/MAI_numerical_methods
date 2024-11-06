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
#include <piecewise_polynomial.h>
#include <algorithms.h>

int get_spline_accuracy(
    piecewise_polynomial spline,
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
            if (std::abs(spline.valueAt(points[i]) - values[i]) > eps)
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
    
    size_t n;
    std::cin >> n;
    
    std::vector<double> points(n), values(n);
    double test_point;
    
    for (size_t i = 0; i < n; ++i)
    {
        std::cin >> points[i] >> values[i];
    }
    std::cin >> test_point;
    
    try
    {
        piecewise_polynomial spline = algorithms::build_spline(points, values);
        spline.println(std::cout, 6);
        
        int accuracy = get_spline_accuracy(spline, points, values);
        std::cout << "Spline is " << (accuracy != -1
                ? std::format("verified (accurate to {} decimal places)", accuracy)
                : "wrong") << std::endl << std::endl;
        
        std::cout << std::format("F({}) = {}", test_point, spline.valueAt(test_point)) << std::endl;
    }
    catch (std::invalid_argument const &e)
    {
        std::cout << e.what() << std::endl;
    }
    
    
    piecewise_polynomial spline = algorithms::build_spline(points, values);
    
    
}