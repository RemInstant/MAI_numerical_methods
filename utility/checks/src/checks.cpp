#include <stdexcept>

#include "../include/checks.h"

void checks::throw_if_invalid_eps(double eps)
{
    if (eps < 0)
    {
        throw std::invalid_argument("Invalid epsilon");
    }
}

void checks::throw_if_invalid_interval(double a, double b, double eps)
{
    if (a - b > -eps)
    {
        throw std::invalid_argument("Invalid interval");
    }
}