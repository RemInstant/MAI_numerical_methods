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
#include <algorithms.h>

double calc_solution(double x)
{
    return std::exp(x*x) + std::exp(x*std::sqrt(2)) + std::exp(-x*std::sqrt(2));
}

std::pair<double, double> calc_equation(double x, double y, double z)
{
    return std::make_pair(z, 2*y + 4*x*x * std::exp(x*x));
}

std::vector<std::tuple<double, double, double>>
solve_cauchy_with_euler(
    double a, double b, double h,
    double y0, double z0,
    std::pair<double, double> (func)(double, double, double))
{
    size_t n = std::round((b - a) / h) + 1;
    
    std::vector<double> vx(n, a);
    std::vector<double> vy(n, y0);
    std::vector<double> vz(n, z0);
    
    for (size_t i = 1; i < n; ++i)
    {
        vx[i] = a + i*h;
        
        auto p = func(vx[i-1], vy[i-1], vz[i-1]);
        vy[i] = vy[i-1] + h * p.first;
        vz[i] = vz[i-1] + h * p.second;
    }
    
    std::vector<std::tuple<double, double, double>> pv(n);
    
    for (size_t i = 0; i < n; ++i)
    {
        pv[i] = std::make_tuple(vx[i], vy[i], vz[i]);
    }
    
    return pv;
}

std::vector<std::tuple<double, double, double>>
solve_cauchy_with_runge_kutta4(
    double a, double b, double h,
    double y0, double z0,
    std::pair<double, double> (func)(double, double, double))
{
    size_t n = std::round((b - a) / h) + 1;
    
    std::vector<double> vx(n, a);
    std::vector<double> vy(n, y0);
    std::vector<double> vz(n, z0);
    
    for (size_t i = 1; i < n; ++i)
    {
        double const &x = vx[i-1];
        double const &y = vy[i-1];
        double const &z = vz[i-1];
        
        auto [k1, l1] = func(x, y, z);
        auto [k2, l2] = func(x + h/2, y + h*k1/2, z + h*l1/2);
        auto [k3, l3] = func(x + h/2, y + h*k2/2, z + h*l2/2);
        auto [k4, l4] = func(x + h, y + h*k3, z + h*l3);
        
        vx[i] = a + i*h;
        vy[i] = vy[i-1] + h/6 * (k1 + 2*k2 + 2*k3 + k4);
        vz[i] = vz[i-1] + h/6 * (l1 + 2*l2 + 2*l3 + l4);
    }
    
    std::vector<std::tuple<double, double, double>> pv(n);
    
    for (size_t i = 0; i < n; ++i)
    {
        pv[i] = std::make_tuple(vx[i], vy[i], vz[i]);
    }
    
    return pv;
}

std::vector<std::tuple<double, double, double>>
solve_cauchy_with_adams4(
    double a, double b, double h,
    double y0, double y1, double y2, double y3,
    double z0, double z1, double z2, double z3,
    std::pair<double, double> (func)(double, double, double))
{
    size_t n = std::round((b - a) / h) + 1;
    
    std::vector<double> vx(n, a);
    std::vector<double> vy = { y0, y1, y2, y3 };
    std::vector<double> vz = { z0, z1, z2, z3 };
    
    vy.resize(n, 0);
    vz.resize(n, 0);
    
    for (size_t i = 1; i < 4; ++i)
    {
        vx[i] = a + i*h;
    }
    
    for (size_t i = 4; i < n; ++i)
    {
        auto [f0, g0] = func(vx[i-1], vy[i-1], vz[i-1]);
        auto [f1, g1] = func(vx[i-2], vy[i-2], vz[i-2]);
        auto [f2, g2] = func(vx[i-3], vy[i-3], vz[i-3]);
        auto [f3, g3] = func(vx[i-4], vy[i-4], vz[i-4]);
        
        vx[i] = a + i*h;
        vy[i] = vy[i-1] + h/24 * (55*f0 - 59*f1 + 37*f2 - 9*f3);
        vz[i] = vz[i-1] + h/24 * (55*g0 - 59*g1 + 37*g2 - 9*g3);
    }
    
    std::vector<std::tuple<double, double, double>> pv(n);
    
    for (size_t i = 0; i < n; ++i)
    {
        pv[i] = std::make_tuple(vx[i], vy[i], vz[i]);
    }
    
    return pv;
}

double calc_error_with_exact(double x, double y)
{
    return std::abs(calc_solution(x) - y);
}

double calc_error_with_rrr(double yh, double y2h, double p)
{
    return std::abs((yh - y2h) / (std::pow(2, p) - 1));
}


int main()
{
    std::freopen("input.txt", "r", stdin);
    
    double a, b, h, y0, z0;
    std::cin >> a >> b >> h >> y0 >> z0;
    
    size_t n = std::round((b-a) / h) + 1;
    
    auto pv_euler = solve_cauchy_with_euler(a, b, h, y0, z0, calc_equation);
    auto pv_runge = solve_cauchy_with_runge_kutta4(a, b, h, y0, z0, calc_equation);
    
    auto pv_adams = solve_cauchy_with_adams4(a, b, h,
            y0, std::get<1>(pv_runge[1]), std::get<1>(pv_runge[2]), std::get<1>(pv_runge[3]),
            z0, std::get<2>(pv_runge[1]), std::get<2>(pv_runge[2]), std::get<2>(pv_runge[3]),
            calc_equation);
    
    auto pv_euler2 = solve_cauchy_with_euler(a, b, h/2, y0, z0, calc_equation);
    auto pv_runge2 = solve_cauchy_with_runge_kutta4(a, b, h/2, y0, z0, calc_equation);
    
    auto pv_adams2 = solve_cauchy_with_adams4(a, b, h/2,
            y0, std::get<1>(pv_runge2[1]), std::get<1>(pv_runge2[2]), std::get<1>(pv_runge2[3]),
            z0, std::get<2>(pv_runge2[1]), std::get<2>(pv_runge2[2]), std::get<2>(pv_runge2[3]),
            calc_equation);
    
    
    std::cout << std::fixed;
    
    std::cout << std::setw(3) << "x" << "  " << std::setw(8) << "exact" << " | " <<
            std::setw(8) << "euler" << "  " << "exact error  " << "rrr error" << " | " <<
            std::setw(8) << "runge" << "  " << "exact error  " << "rrr error" << " | " <<
            std::setw(8) << "adams" << "  " << "exact error  " << "rrr error" << std::endl;
    
    for (size_t i = 0; i < n; ++i)
    {
        double x = a + i*h;
        double y = calc_solution(x);
        
        
        std::cout.precision(1);
        std::cout << x << "  ";
        std::cout.precision(6);
        std::cout << y << " | ";
        
        std::cout << std::get<1>(pv_euler[i]) << "  ";
        std::cout << std::setw(11) << calc_error_with_exact(x, std::get<1>(pv_euler[i])) << "  ";
        std::cout << std::setw(9) << (calc_error_with_rrr(std::get<1>(pv_euler2[2*i]),
                std::get<1>(pv_euler[i]), 2)) << " | ";
                
        std::cout << std::get<1>(pv_runge[i]) << "  ";
        std::cout << std::setw(11) << calc_error_with_exact(x, std::get<1>(pv_runge[i])) << "  ";
        std::cout << std::setw(9) << (calc_error_with_rrr(std::get<1>(pv_runge2[2*i]),
                std::get<1>(pv_runge[i]), 4)) << " | ";
        
        std::cout << std::get<1>(pv_adams[i]) << "  ";
        std::cout << std::setw(11) << calc_error_with_exact(x, std::get<1>(pv_adams[i])) << "  ";
        std::cout << std::setw(9) << (calc_error_with_rrr(std::get<1>(pv_adams2[2*i]),
                std::get<1>(pv_adams[i]), 4)) << std::endl;
    }
}