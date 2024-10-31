#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <cmath>
#include <numbers>

#include <checks.h>

double func(double x)
{
    return std::sqrt(1 - x*x) - std::exp(x) + 0.1;
}

double func_first_deriv(double x)
{
    return -x / std::sqrt(1 - x*x) - std::exp(x);
}

double func_second_deriv(double x)
{
    return -std::pow(1 - x*x, -3.0 / 2) - std::exp(x);
}

double phi(double x)
{
    return std::log(std::sqrt(1 - x*x) + 0.1);
}

std::pair<double, size_t>
solve_equation_with_iterations(double (*phi)(double), double q, double x0, double eps)
{
    checks::throw_if_invalid_eps(eps);
    
    double x = x0;
    double x_prev;
    double eps_mult = q * (1 - q);
    size_t iter = 0;
    
    do
    {
        x_prev = x;
        x = phi(x);
        ++iter;
    } while (eps_mult * std::abs(x - x_prev) > eps);
    
    return std::make_pair(x, iter);
}

std::pair<double, size_t>
solve_equation_with_dichotomy(double (*f)(double), double a, double b, double eps)
{
    checks::throw_if_invalid_eps(eps);
    checks::throw_if_invalid_interval(a, b, eps);
    
    size_t iter = 0;
    
    do
    {
        double mid = (a + b) / 2;
        double a_sign = f(a) > 0 ? 1 : -1;
        double b_sign = f(b) > 0 ? 1 : -1;
        double mid_sign = f(mid) > 0 ? 1 : -1;
        
        if (a_sign != mid_sign && mid_sign == b_sign)
        {
            b = mid;
        }
        else if (b_sign != mid_sign && mid_sign == a_sign)
        {
            a = mid;
        }
        else if (a_sign != mid_sign && mid_sign != b_sign)
        {
            throw std::invalid_argument("There are TWO roots on the interval");
        }
        else if (a_sign == mid_sign && mid_sign == b_sign)
        {
            throw std::invalid_argument("There are NO roots on the interval");
        }
        
        ++iter;
    } while (b - a > 2 * eps);
    
    return std::make_pair((a + b) / 2, iter);
}

std::pair<double, size_t>
solve_equation_with_newton(
    double (*func)(double),
    double (*func_first_deriv)(double),
    double (*func_second_deriv)(double),
    double x0, double eps)
{
    checks::throw_if_invalid_eps(eps);
    
    if (func(x0) * func_second_deriv(x0) < std::numeric_limits<double>::epsilon())
    {
        throw std::invalid_argument("Newton method will not converge with the given approximation");
    }
    
    double x = x0;
    double x_prev;
    size_t iter = 0;
    
    do
    {
        x_prev = x;
        x -= func(x) / func_first_deriv(x);
        ++iter;
    } while (std::abs(x - x_prev) > eps);
    
    return std::make_pair(x, iter);
}

std::pair<double, size_t>
solve_equation_with_secant(
    double (*func)(double),
    double (*func_second_deriv)(double),
    double x0, double x1, double eps)
{
    checks::throw_if_invalid_eps(eps);
    
    if (func(x0) * func_second_deriv(x0) < std::numeric_limits<double>::epsilon() ||
            func(x1) * func_second_deriv(x1) < std::numeric_limits<double>::epsilon())
    {
        throw std::invalid_argument("Newton method will not converge with the given approximation");
    }
    
    double x = x1;
    double x_prev = x0;
    size_t iter = 0;
    
    do
    {
        double tmp = x;
        x -= func(x) * (x - x_prev) / (func(x) - func(x_prev));
        x_prev = tmp;
        ++iter;
    } while (std::abs(x - x_prev) > eps);
    
    return std::make_pair(x, iter);
}

int main()
{
    std::freopen("input.txt", "r", stdin);
       
    size_t n = 0;
    double eps;
    
    std::cin >> eps;
    std::cout.precision(15);
    std::cout.setf(std::ios_base::fixed);
    
    std::cout << "          eps = " << eps << std::endl;
    
    double a_iter = 0.05;
    double b_iter = 0.1;
    double x0_iter = (a_iter + b_iter) / 2;
    double q_iter = 0.1;
    
    double a_dich = 0;
    double b_dich = 1;
    
    double x0_newton = 0.12;
    
    double x0_secant = 0.12;
    double x1_secant = 0.1;
    
    try
    {
        auto [x_iter, k_iter] = solve_equation_with_iterations(phi, q_iter, x0_iter, eps);
        std::cout << "Iterations: x = " << x_iter << " (" << k_iter << ")" << std::endl;
    }
    catch (std::invalid_argument e)
    {
        std::cout << "Iterations: " << e.what() << std::endl;
    }
    
    try
    {
        auto [x_dich, k_dich] = solve_equation_with_dichotomy(func, a_dich, b_dich, eps);
        std::cout << "Dichotomy:  x = " << x_dich << " (" << k_dich << ")" << std::endl;
    }
    catch (std::invalid_argument e)
    {
        std::cout << "Dichotomy: " << e.what() << std::endl;
    }
    
    try
    {
        auto [x_newton, k_newton] = solve_equation_with_newton(
                func, func_first_deriv, func_second_deriv, x0_newton, eps);
        std::cout << "Newton:     x = " << x_newton << " (" << k_newton << ")" << std::endl;
    }
    catch (std::invalid_argument e)
    {
        std::cout << "Newton: " << e.what() << std::endl;
    }
    
    try
    {
        auto [x_secant, k_secant] = solve_equation_with_secant(
                func, func_second_deriv, x0_secant, x1_secant, eps);
        std::cout << "Secant:     x = " << x_secant << " (" << k_secant << ")" << std::endl;
    }
    catch (std::invalid_argument e)
    {
        std::cout << "Secant: " << e.what() << std::endl;
    }
    
    try
    {
        std::vector<size_t> iter_cnt, dich_cnt, newton_cnt, secant_cnt;
        eps = 1;
        
        std::cout << std::endl;
        std::cout << "EPS        | ";
        for (size_t i = 0; i < 16; ++i)
        {
            iter_cnt.push_back(solve_equation_with_iterations(phi, q_iter, x0_iter, eps).second);
            dich_cnt.push_back(solve_equation_with_dichotomy(func, a_dich, b_dich, eps).second);
            newton_cnt.push_back(solve_equation_with_newton(
                    func, func_first_deriv, func_second_deriv, x0_newton, eps).second);
            secant_cnt.push_back(solve_equation_with_secant(
                    func, func_second_deriv, x0_secant, x1_secant, eps).second);
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
        print_vec("Dichotomy  |", dich_cnt);
        print_vec("Newton:    |", newton_cnt);
        print_vec("Secant:    |", secant_cnt);
    }
    catch (std::invalid_argument e)
    {
        std::cout << "Exception occurred in one of the methods" << std::endl << e.what() << std::endl;
    }
}