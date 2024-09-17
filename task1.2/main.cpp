#include <iostream>
#include <iomanip>
#include <vector>

void enter_coefs(
        std::vector<double> &coefs,
        char identifier,
        size_t from,
        size_t to)
{
    std::cout << "Enter coefficients from " << identifier << (from + 1) << " to " << identifier << (to + 1) << ": ";
    
    for (size_t i = from; i <= to; ++i)
    {
        std::cin >> coefs[i];
    }
}

// void print_system(std::vector<double> &a, std::vector<double> &b, std::vector<double> &c, std::vector<double> &d, )
// {
//     std::cout << std::setprecision(0);
//     std::cout << b[0] << "x1"
    
// }

bool check_stability(
        std::vector<double> const &coefs_a,
        std::vector<double> const &coefs_b,
        std::vector<double> const &coefs_c)
{
    size_t n = coefs_a.size();
    
    for (size_t i = 0; i < n; ++i)
    {
        if (abs(coefs_b[i]) < abs(coefs_a[i]) + abs(coefs_c[i]))
        {
            return false;
        }
    }
    
    return true;
}

std::pair<std::vector<double>, std::vector<double>>
process_thomas_alg_first_phase(
        std::vector<double> const &coefs_a,
        std::vector<double> const &coefs_b,
        std::vector<double> const &coefs_c,
        std::vector<double> const &coefs_d,
        double eps)
{
    size_t n = coefs_a.size();
    std::vector<double> coefs_p(n, 0), coefs_q(n, 0);
    
    if (abs(coefs_b[0]) < eps)
    {
        throw std::invalid_argument("Incorrect system");
    }
    
    coefs_p[0] = -coefs_c[0] / coefs_b[0];
    coefs_q[0] = coefs_d[0] / coefs_b[0];
    
    for (size_t i = 1; i < n; ++i)
    {
        double div = coefs_b[i] + coefs_a[i] * coefs_p[i - 1];
        
        if (abs(div) < eps)
        {
            throw std::invalid_argument("Incorrect system");
        }
        
        coefs_p[i] = -coefs_c[i] / div;
        coefs_q[i] = (coefs_d[i] - coefs_a[i] * coefs_q[i - 1]) / div;
    }
    
    return std::make_pair(coefs_p, coefs_q);
}

std::vector<double>
process_thomas_alg_second_phase(
        std::vector<double> const &coefs_p,
        std::vector<double> const &coefs_q)
{
    size_t n = coefs_p.size();
    std::vector<double> roots(n, 0);
    
    roots[n - 1] = coefs_q[n - 1];
    
    for (size_t i = n - 1; i > 0; --i)
    {
        roots[i - 1] = coefs_p[i - 1] * roots[i] + coefs_q[i - 1];
    }
    
    return roots;
}

double calc_determinant(
        std::vector<double> const &coefs_a,
        std::vector<double> const &coefs_b,
        std::vector<double> const &coefs_p)
{
    size_t n = coefs_a.size();
    double det = coefs_b[0];
    
    for (size_t i = 1; i < n; ++i)
    {
        det *= coefs_b[i] + coefs_a[i] * coefs_p[i - 1];
    }
    
    return det;
}

int main()
{
    double const EPS = 1e-8;
    
    size_t n = 0;
    std::cout << "Enter dimension: ";
    std::cin >> n;
    
    std::vector<double> coefs_a(n, 0), coefs_b(n, 0), coefs_c(n, 0), coefs_d(n, 0);
    
    enter_coefs(coefs_a, 'a', 1, n - 1);
    enter_coefs(coefs_b, 'b', 0, n - 1);
    enter_coefs(coefs_c, 'c', 0, n - 2);
    enter_coefs(coefs_d, 'd', 0, n - 1);

    bool stability = check_stability(coefs_a, coefs_b, coefs_c);
    
    if (stability)
    {
        std::cout << "The method is stable for the given system" << std::endl;
    }
    else
    {
        std::cout << "The method is not stable for the given system" << std::endl;
    }
    
    std::vector<double> roots(0);
    double determinant = 0;    
    
    try
    {
        auto [coefs_p, coefs_q] = process_thomas_alg_first_phase(coefs_a, coefs_b, coefs_c, coefs_d, EPS);   
        roots = process_thomas_alg_second_phase(coefs_p, coefs_q);
        determinant = calc_determinant(coefs_a, coefs_b, coefs_p);
    }
    catch (const std::invalid_argument &e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }
    
    for (size_t i = 0; i < n; ++i)
    {
        std::cout << std::fixed << std::setprecision(8);
        std::cout << 'x' << (i+1) << " = " << roots[i] << std::endl;
    }
    
    std::cout << "Determinant = " << determinant << std::endl;
}