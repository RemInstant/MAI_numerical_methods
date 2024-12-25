#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <cmath>
#include <numbers>
#include <format>
#include <algorithm>

#include <checks.h>
#include <algorithms.h>

double func(double x)
{
    return x / std::pow(3*x + 4, 3);
}

double integrate_with_rectangle_method(double (func)(double), double a, double b, double h)
{
    size_t n = std::round((b - a) / h);
    
    double sum = 0;
    
    for (size_t i = 0; i < n; ++i)
    {
        double l = a + h*i;
        double r = a + h*(i+1);
        
        sum += func((l + r) / 2);
    }
    
    return sum * h;
}

double integrate_with_trapezium_method(double (func)(double), double a, double b, double h)
{
    size_t n = std::round((b - a) / h);
    
    double sum = 0;
    
    for (size_t i = 0; i < n; ++i)
    {
        double l = a + h*i;
        double r = a + h*(i+1);
        
        sum += func(l) + func(r);
    }
    
    return sum * h / 2;
}

double integrate_with_simpson_method(double (func)(double), double a, double b, double h)
{
    size_t n = std::round((b - a) / h);
    
    double sum = 0;
    
    for (size_t i = 0; i < n; i += 2)
    {
        double l = a + h*i;
        double m = a + h*(i+1);
        double r = a + h*(i+2);
        
        sum += func(l) + 4*func(m) + func(r);
    }
    
    return sum * h / 3;
}

double integrate_with_rrr_method(double f_h, double f_kh, double k, double p)
{
    return f_h + (f_h - f_kh) / (std::pow(k, p) - 1);
}


int main()
{
    std::freopen("input.txt", "r", stdin);
    
    double a, b, h1, h2;
    std::cin >> a >> b >> h1 >> h2;
    
    double solution = -0.122448979592;
    
    double integr_rectangle_1 = integrate_with_rectangle_method(func, a, b, h1);
    double integr_trapezium_1 = integrate_with_trapezium_method(func, a, b, h1);
    double integr_simpson_1 = integrate_with_simpson_method(func, a, b, h1);
    
    double integr_rectangle_2 = integrate_with_rectangle_method(func, a, b, h2);
    double integr_trapezium_2 = integrate_with_trapezium_method(func, a, b, h2);
    double integr_simpson_2 = integrate_with_simpson_method(func, a, b, h2);
    
    double integr_rrr_rectangle = integrate_with_rrr_method(
            integr_rectangle_1, integr_rectangle_2, h2/h1, 2);
    double integr_rrr_trapezium = integrate_with_rrr_method(
            integr_trapezium_1, integr_trapezium_2, h2/h1, 2);
    double integr_rrr_simpson = integrate_with_rrr_method(
            integr_simpson_1, integr_simpson_2, h2/h1, 2);
    
    std::cout.precision(6);
    std::cout << std::fixed;
    
    std::cout << "Solution = " << solution << std::endl << std::endl;
    
    std::cout << "METHOD     " << std::setw(9) << "VALUE" << "  " << 
        std::setw(9) << "ERROR" << std::endl;
    
    std::cout << "h = " << h1 << std::endl;
    std::cout << "Rectangle  " << std::setw(9) << integr_rectangle_1 << "  " << std::endl;
    std::cout << "Trapezium  " << std::setw(9) << integr_trapezium_1 << "  " << std::endl;
    std::cout << "Simpson    " << std::setw(9) << integr_simpson_1 << "  " << std::endl;
            
    std::cout << "h = " << h2 << std::endl;
    std::cout << "Rectangle  " << std::setw(9) << integr_rectangle_2 << "  " << std::endl;
    std::cout << "Trapezium  " << std::setw(9) << integr_trapezium_2 << "  " << std::endl;
    std::cout << "Simpson    " << std::setw(9) << integr_simpson_2 << "  " << std::endl;
    
    std::cout << "RRR improvement" << std::endl;
    std::cout << "Rectangle  " << std::setw(9) << integr_rrr_rectangle << "  " <<
                std::abs(solution - integr_rrr_rectangle) << std::endl;
    std::cout << "Trapezium  " << std::setw(9) << integr_rrr_trapezium << "  " <<
                std::abs(solution - integr_rrr_trapezium) << std::endl;
    std::cout << "Simpson    " << std::setw(9) << integr_rrr_simpson << "  " <<
                std::abs(solution - integr_rrr_simpson) << std::endl;
    
}