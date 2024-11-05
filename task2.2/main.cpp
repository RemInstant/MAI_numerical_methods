#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <cmath>
#include <numbers>
#include <format>

#include <checks.h>
#include <vecN.h>
#include <matrixNxN.h>
#include <algorithms.h>

double f1(double x1, double x2)
{
    return (x1*x1 + 16) * x2 - 64;
}

double f2(double x1, double x2)
{
    return std::pow(x1 - 2, 2) + std::pow(x2 - 2, 2) - 16;
}


double f1_x1_deriv(double x1, double x2)
{
    return 2*x1*x2;
}

double f1_x2_deriv(double x1, double x2)
{
    return x1*x1 + 16;
}

double f2_x1_deriv(double x1, double x2)
{
    return 2 * (x1 - 2);
}

double f2_x2_deriv(double x1, double x2)
{
    return 2 * (x2 - 2);
}


double phi1(double x1, double x2)
{
    return std::sqrt(16 - std::pow(x2 - 2, 2)) + 2;
}

double phi2(double x1, double x2)
{
    return 64 / (x1*x1 + 16);
}

std::pair<std::pair<double, double>, size_t>
solve_equation_system_with_iterations(
    double (*phi1)(double, double),
    double (*phi2)(double, double),
    double q, double x1_0, double x2_0, double eps)
{
    checks::throw_if_invalid_eps(eps);
    
    double x1 = x1_0;
    double x2 = x2_0;
    double x1_prev;
    double x2_prev;
    double eps_mult = q * (1 - q);
    double check;
    size_t iter = 0;
    
    do
    {
        x1_prev = x1;
        x2_prev = x2;
        double x1_tmp = phi1(x1, x2);
        double x2_tmp = phi2(x1, x2);
        x1 = x1_tmp;
        x2 = x2_tmp;
        
        check = std::max(std::abs(x1 - x1_prev), std::abs(x2 - x2_prev));
        ++iter;
    } while (eps_mult * check > eps);
    
    return std::make_pair(std::make_pair(x1, x2), iter);
}

std::pair<std::pair<double, double>, size_t>
solve_equation_system_with_newton(
    double (*f1)(double, double),
    double (*f2)(double, double),
    double (*f1_x1_deriv)(double, double),
    double (*f1_x2_deriv)(double, double),
    double (*f2_x1_deriv)(double, double),
    double (*f2_x2_deriv)(double, double),
    double x1_0, double x2_0, double eps)
{
    checks::throw_if_invalid_eps(eps);
    
    double x1 = x1_0;
    double x2 = x2_0;
    double x1_prev;
    double x2_prev;
    double check;
    size_t iter = 0;
    
    do
    {
        x1_prev = x1;
        x2_prev = x2;
        
        matrixNxN coefs(2);
        vecN constant_terms(2);
        coefs[0][0] = f1_x1_deriv(x1, x2);
        coefs[0][1] = f1_x2_deriv(x1, x2);
        coefs[1][0] = f2_x1_deriv(x1, x2);
        coefs[1][1] = f2_x2_deriv(x1, x2);
        constant_terms[0] = -f1(x1, x2);
        constant_terms[1] = -f2(x1, x2);
        
        vecN delta = algorithms::solve_linear_equation(coefs, constant_terms);
        
        x1 += delta[0];
        x2 += delta[1];
        
        check = std::max(std::abs(x1 - x1_prev), std::abs(x2 - x2_prev));
        ++iter;
    } while (check > eps);
    
    return std::make_pair(std::make_pair(x1, x2), iter);
}

int get_solution_accuracy(
    double (*f1)(double, double),
    double (*f2)(double, double),
    std::pair<double, double> x)
{
    double x1 = x.first;
    double x2 = x.second;
    double eps = 1.0;
    size_t accuracy = 0;
    
    for (; eps > std::numeric_limits<double>::epsilon(); ++accuracy, eps /= 10)
    {
        if (std::abs(f1(x1, x2)) > eps || std::abs(f2(x1, x2)) > eps)
        {
            return accuracy - 1;
        }
    }
    
    return accuracy;
}

void print_verifictation(int accuracy)
{
    std::cout << (accuracy != -1
            ? std::format("Verified (accurate to {} decimal places)", accuracy)
            : "Wrong") << std::endl;
}

int main()
{
    std::freopen("input.txt", "r", stdin);
       
    size_t n = 0;
    double eps;
    
    std::cin >> eps;
    std::cout.precision(15);
    std::cout.setf(std::ios_base::fixed);
    
    std::cout << "           eps = " << eps << "  " << eps << std::endl;
    
    double x1_0_iter = 6;
    double x2_0_iter = 2;
    double q_iter = 0.6;
    
    double x1_0_newton = 6;
    double x2_0_newton = 2;
    
    try
    {
        auto [x, k_iter] = solve_equation_system_with_iterations(
                phi1, phi2, q_iter, x1_0_iter, x2_0_iter, eps);
        std::cout << "Iterations: x = (" << x.first << ", " << x.second << ")  ("
                << std::setw(2) << k_iter << ")  ";
        print_verifictation(get_solution_accuracy(f1, f2, x));
    }
    catch (std::invalid_argument e)
    {
        std::cout << "Iterations: " << e.what() << std::endl;
    }
    
    try
    {
        auto [x, k_iter] = solve_equation_system_with_newton(
                f1, f2, f1_x1_deriv, f1_x2_deriv, f2_x1_deriv, f2_x2_deriv,
                x1_0_newton, x2_0_newton, eps);
        std::cout << "Newton:     x = (" << x.first << ", " << x.second << ")  ("
                << std::setw(2) << k_iter << ")  ";
        print_verifictation(get_solution_accuracy(f1, f2, x));
    }
    catch (std::invalid_argument e)
    {
        std::cout << "Newton: " << e.what() << std::endl;
    }
    
    
    try
    {
        std::vector<size_t> iter_cnt, newton_cnt;
        eps = 1;
        
        std::cout << std::endl;
        std::cout << "EPS        | ";
        for (size_t i = 0; i < 16; ++i)
        {
            iter_cnt.push_back(solve_equation_system_with_iterations(
                phi1, phi2, q_iter, x1_0_iter, x2_0_iter, eps).second);
            newton_cnt.push_back(solve_equation_system_with_newton(
                    f1, f2, f1_x1_deriv, f1_x2_deriv, f2_x1_deriv, f2_x2_deriv,
                    x1_0_newton, x2_0_newton, eps).second);
            eps /= 10;
        }
        
        for (size_t i = 0; i < iter_cnt.size(); ++i)
        {
            std::cout << std::setw(4) << ("e-" + std::to_string(i)) << ' ';
        }
        std::cout << std::endl;
        
        void (*print_vec)(std::string, std::vector<size_t>) =
                [](std::string str, std::vector<size_t> data)
                {
                    std::cout << str << ' ';
                    for (auto elem : data)
                    {
                        std::cout << std::setw(4) << elem << ' ';
                    }
                    std::cout << std::endl;
                };
        
        print_vec("Iterations |", iter_cnt);
        print_vec("Newton:    |", newton_cnt);
    }
    catch (std::invalid_argument e)
    {
        std::cout << "Exception occurred in one of the methods" << std::endl << e.what() << std::endl;
    }
}