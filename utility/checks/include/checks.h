#ifndef _NUM_METHODS_UTIL_CHECKS_H_
#define _NUM_METHODS_UTIL_CHECKS_H_

namespace checks
{
    void throw_if_zero_divisor(double divisor);
    void throw_if_invalid_eps(double eps);
    void throw_if_invalid_interval(double a, double b, double eps);
}

#endif // _NUM_METHODS_UTIL_CHECKS_H_