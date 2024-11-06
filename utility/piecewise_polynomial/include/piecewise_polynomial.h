#ifndef _NUM_METHODS_UTIL_PIECEWISE_POLYNOMIAL_H_
#define _NUM_METHODS_UTIL_PIECEWISE_POLYNOMIAL_H_

#include <deque>
#include <vector>
#include <iostream>

#include <polynomial.h>

class piecewise_polynomial 
{

private:

    std::deque<double> _bounds;
    std::deque<polynomial> _polynomials;
    size_t _degree;
    
public:
    piecewise_polynomial();
    piecewise_polynomial(double left_bound, double right_bound, polynomial const &poly);

public:
    void print(std::ostream &stream = std::cout, size_t precision = 3, size_t dot_precision = 3) const;
    void println(std::ostream &stream = std::cout, size_t precision = 3, size_t dot_precision = 3) const;

public:
    double operator()(double x) const;
    double operator[](size_t idx) const;

public:
    piecewise_polynomial &push_front(double delim, polynomial const &poly);
    piecewise_polynomial &push_back(double delim, polynomial const &poly);
    
    double get_bound(size_t idx) const;
    void set_bound(size_t idx, double x);

public:
    size_t bound_count() const;
    size_t degree() const;
    double valueAt(double x) const;

};

#endif // _NUM_METHODS_UTIL_PIECEWISE_POLYNOMIAL_H_