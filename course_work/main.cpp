#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <limits>
#include <cmath>
#include <numbers>
#include <format>
#include <algorithm>
#include <functional>

#include <sys/wait.h>
#include <unistd.h>

#include <checks.h>
#include <vecN.h>
#include <algorithms.h>

struct gas_params {
    double x;
    double speed;
    double pressure;
    double density;
    double mach_number;
};

double calc_sound_speed(double pressure, double density, double adiabatic_index) {
    return std::sqrt(adiabatic_index * pressure / density);
}

double calc_mach_number(double speed, double pressure, double density, double adiabatic_index) {
    return speed / calc_sound_speed(pressure, density, adiabatic_index);
}


std::vector<gas_params> calc_gas_params_with_runge_kutta4(
    double x0, double x1, double h,
    double speed0, double pressure0, double density0,
    double adiabatic_index,
    std::function<double (double)> calc_square,
    std::function<double (double)> calc_square_deriv)
{
    auto calc_step = [calc_square, calc_square_deriv, adiabatic_index](
        double x, double speed, double pressure, double density)
    {
        double m = calc_mach_number(speed, pressure, density, adiabatic_index);
        
        double f = calc_square(x);
        double fd = calc_square_deriv(x);
        
        double k_speed = (speed * fd) / ((m*m - 1) * f);
        double k_pressure = -(density * speed * speed * fd) / ((m*m - 1) * f);
        double k_density = -(m * m * density * fd) / ((m*m - 1) * f);
        
        return std::make_tuple(k_speed, k_pressure, k_density);
    };
    
    size_t n = (x1 - x0) / h + 1;
    if (calc_square_deriv(x0) < 0)
    {
        n = std::numeric_limits<size_t>::max();
    }
    
    std::vector<gas_params> params(0);
    
    double mach0 = calc_mach_number(speed0, pressure0, density0, adiabatic_index);
    params.push_back( { x0, speed0, pressure0, density0, mach0 } );
    
    for (size_t i = 1; i < n; ++i)
    {
        double const &x = params[i-1].x;
        double const &speed = params[i-1].speed;
        double const &pressure = params[i-1].pressure;
        double const &density = params[i-1].density;
        double const &mach = params[i-1].mach_number;
        
        auto [s_k1, p_k1, d_k1] = calc_step(x, speed, pressure, density);
        auto [s_k2, p_k2, d_k2] = calc_step(x + h/2, speed + h*s_k1/2, pressure + h*p_k1/2, density + h*d_k1/2);
        auto [s_k3, p_k3, d_k3] = calc_step(x + h/2, speed + h*s_k2/2, pressure + h*p_k2/2, density + h*d_k2/2);
        auto [s_k4, p_k4, d_k4] = calc_step(x + h, speed + h*s_k3, pressure + h*p_k3, density + h*d_k3);
        double next_x        = x0 + i*h;
        double next_speed    = params[i-1].speed    + h/6 * (s_k1 + 2*s_k2 + 2*s_k3 + s_k4);
        double next_pressure = params[i-1].pressure + h/6 * (p_k1 + 2*p_k2 + 2*p_k3 + p_k4);
        double next_density  = params[i-1].density  + h/6 * (d_k1 + 2*d_k2 + 2*d_k3 + d_k4);
        double next_mach     = calc_mach_number(next_speed, next_pressure, next_density, adiabatic_index);
        
        if ((calc_square_deriv(x0) < 1) && (
                (mach0 < 1 && next_mach > 1 - std::numeric_limits<double>::epsilon()) ||
                (mach0 > 1 && next_mach < 1 + std::numeric_limits<double>::epsilon()) ||
                (mach0 < 1 && next_mach < mach) ||
                (mach0 > 1 && next_mach > mach)))
        {
            return params;
        }
        
        params.push_back( { next_x, next_speed, next_pressure, next_density, next_mach } );
    }
    
    return params;
}


int main()
{
    std::freopen("input.txt", "r", stdin);
    // x0 x1 r0 a h
    // speed0
    // pressure0
    // density0
    // adiabatic index
    
    double x0, x1, r0, a, h;
    double speed0, pressure0, density0;
    double adiabatic_index;
    
    std::cin >> x0 >> x1 >> r0 >> a >> h >> speed0 >> pressure0 >> density0 >> adiabatic_index;
    
    double mach0 = calc_mach_number(speed0, pressure0, density0, adiabatic_index);
    
    std::cout << a << std::endl;
    
    std::function<double (double)> calc_radius = [r0, a, x0](double x)
    {
        return r0 + a * (x - x0);
    };
    
    std::function<double (double)> calc_square = [calc_radius](double x) -> double
    {
        return std::numbers::pi * std::pow(calc_radius(x), 2);
    };
    
    std::function<double (double)> calc_square_deriv = [a, calc_radius](double x) -> double
    {
        return 2 * std::numbers::pi * a * calc_radius(x);
    };
    
    auto params = calc_gas_params_with_runge_kutta4(x0, x1, h, speed0, pressure0, density0, adiabatic_index,
                                                    calc_square, calc_square_deriv);
    
    double eq_const1 = pressure0 / std::pow(density0, adiabatic_index);
    double eq_const2 = density0 * speed0 * calc_square(x0);
    
    std::cout << 
        std::endl <<
        "* - величины полученные из скорости с помощью алгебраических соотношений" <<
        std::endl << std::endl;
    
    std::cout << 
            std::fixed <<
            std::setw(8) << "x" << " | " <<
            std::setw(8) << "Скорость" << " | " <<
            std::setw(8) << "Давление" << " | " <<
            std::setw(9) << "Плотность" << " | " <<
            std::setw(10) << "Число Маха" << " | " <<
            " | " <<
            std::setw(10) << "Давление*" << " | " <<
            std::setw(10) << "Плотность*" << " | " <<
            std::endl;
    
    for (size_t i = 0; i < params.size(); ++i)
    {
        if (i == 50 && params.size() > 150) {
            i = params.size() - 50;
        }
        
        double check_pressure = eq_const1 * std::pow(params[i].density, 1.4);
        double check_density = eq_const2 / params[i].speed / calc_square(params[i].x);
        
        double temperature = params[i].pressure * 0.02898 / params[i].density / 8.314;
        double u = 5/2 * 8.314 * temperature;
        double h = u + params[i].pressure;
        
        double test = params[i].density * params[i].speed * calc_square(params[i].x);
        
        std::cout <<
            std::setprecision(6) << std::setw(5)  << params[i].x           << " | " <<
            std::setprecision(2) << std::setw(8)  << params[i].speed       << " | " <<
            std::setprecision(1) << std::setw(8)  << params[i].pressure    << " | " <<
            std::setprecision(6) << std::setw(9)  << params[i].density     << " | " <<
            std::setprecision(6) << std::setw(10) << params[i].mach_number << " | " <<
            " | " <<
            std::setprecision(1) << std::setw(9)  << check_pressure << " | " <<
            std::setprecision(6) << std::setw(10)  << check_density << " | " <<
            std::setprecision(6) << std::setw(10)  << test << " | " <<
            std::endl;
    }
    
    std::cout << "Итераций: " << params.size() << std::endl;
    
    std::ofstream out("buf");
    
    
    out << r0 << ' ' << a << std::endl;
    out << params.size() << std::endl;
    for (size_t i = 0; i < params.size(); ++i) out << params[i].x << ' '; out << std::endl;
    for (size_t i = 0; i < params.size(); ++i) out << params[i].speed << ' '; out << std::endl;
    for (size_t i = 0; i < params.size(); ++i) out << params[i].pressure << ' '; out << std::endl;
    for (size_t i = 0; i < params.size(); ++i) out << params[i].density << ' '; out << std::endl;
    for (size_t i = 0; i < params.size(); ++i) out << params[i].mach_number << ' '; out << std::endl;
    
    out.flush();
    
    system("python3 graphs.py");
}