#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <cmath>
#include <numbers>
#include <format>
#include <algorithm>

#include <checks.h>
#include <piecewise_polynomial.h>
#include <polynomial.h>
#include <algorithms.h>

piecewise_polynomial calc_first_derivative1(
    std::vector<double> const &points,
    std::vector<double> const &values)
{
    if (points.size() != values.size())
    {
        throw std::invalid_argument("Point and values count are not equal");
    }
    if (points.size() < 2)
    {
        return piecewise_polynomial(0, 0, polynomial({0}));
    }
    
    piecewise_polynomial deriv(points[0], points[1], polynomial(
            {(values[1] - values[0]) / (points[1] - points[0])}));
    
    for (size_t i = 2; i < points.size(); ++i)
    {
        deriv.push_back(points[i], polynomial(
                {(values[i] - values[i-1]) / (points[i] - points[i-1])}));
    }
    
    return deriv;
}

piecewise_polynomial calc_first_derivative2(
    std::vector<double> const &points,
    std::vector<double> const &values)
{
    if (points.size() != values.size())
    {
        throw std::invalid_argument("Point and values count are not equal");
    }
    if (points.size() < 3)
    {
        return piecewise_polynomial(0, 0, polynomial({0}));
    }
    
    piecewise_polynomial deriv;
    
    for (size_t i = 2; i < points.size(); ++i)
    {
        
        double coef_A = (values[i  ] - values[i-1]) / (points[i  ] - points[i-1]);
        double coef_B = (values[i-1] - values[i-2]) / (points[i-1] - points[i-2]);
        double coef_C = points[i] - points[i-2];
        
        double coef = (coef_A - coef_B) / coef_C;
        
        polynomial p({-points[i-2] - points[i-1], 2});
        p *= polynomial({coef});
        p += polynomial({coef_B});
        
        if (i == 2)
        {
            deriv = piecewise_polynomial(points[i-2], points[i-1], p);
        }
        else
        {
            deriv.push_back(points[i-1], p);
        }
    }
    
    return deriv;
}

piecewise_polynomial calc_second_derivative(
    std::vector<double> const &points,
    std::vector<double> const &values)
{
    if (points.size() != values.size())
    {
        throw std::invalid_argument("Point and values count are not equal");
    }
    if (points.size() < 3)
    {
        return piecewise_polynomial(0, 0, polynomial({0}));
    }
    
    piecewise_polynomial deriv(points[0], points[1], polynomial(
            {(values[1] - values[0]) / (points[1] - points[0])}));
    
    for (size_t i = 2; i < points.size(); ++i)
    {
        
        double coef_A = (values[i  ] - values[i-1]) / (points[i  ] - points[i-1]);
        double coef_B = (values[i-1] - values[i-2]) / (points[i-1] - points[i-2]);
        double coef_C = points[i] - points[i-2];
        
        double coef = (coef_A - coef_B) / coef_C;
        polynomial p({2 * coef});
        
        if (i == 2)
        {
            deriv = piecewise_polynomial(points[i-2], points[i-1], p);
        }
        else
        {
            deriv.push_back(points[i-1], p);
        }
    }
    
    return deriv;
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
    std::cout << "TEST POINT X* = " << test_point << std::endl << std::endl;
    
    try
    {
        
        piecewise_polynomial deriv = calc_first_derivative1(points, values);
        std::cout << "1th-order accurate first derivative" << std::endl;
        deriv.println(std::cout, 6, 1);
        std::cout << "F(X*) = " << deriv.valueAt(test_point) << " / " <<
                deriv(test_point + 2*std::numeric_limits<double>::epsilon()) << std::endl;
        std::cout << std::endl;
    }
    catch (std::invalid_argument const &e)
    {
        std::cout << e.what() << std::endl;
    }
    
    try
    {
        
        piecewise_polynomial deriv = calc_first_derivative2(points, values);
        std::cout << "2th-order accurate first derivative" << std::endl;
        deriv.println(std::cout, 6, 1);
        std::cout << "F(X*) = " << deriv.valueAt(test_point) << " / " <<
                deriv(test_point + 2*std::numeric_limits<double>::epsilon()) << std::endl;
        std::cout << std::endl;
    }
    catch (std::invalid_argument const &e)
    {
        std::cout << e.what() << std::endl;
    }
    
    try
    {
        
        piecewise_polynomial deriv = calc_second_derivative(points, values);
        std::cout << "1th-order accurate second derivative" << std::endl;
        deriv.println(std::cout, 6, 1);
        std::cout << "F(X*) = " << deriv.valueAt(test_point) << " / " <<
                deriv(test_point + 2*std::numeric_limits<double>::epsilon()) << std::endl;
        std::cout << std::endl;
    }
    catch (std::invalid_argument const &e)
    {
        std::cout << e.what() << std::endl;
    }
}