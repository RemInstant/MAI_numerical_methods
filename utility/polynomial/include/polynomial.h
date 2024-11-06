#ifndef _NUM_METHODS_UTIL_POLYNOMIAL_H_
#define _NUM_METHODS_UTIL_POLYNOMIAL_H_

#include <deque>
#include <vector>
#include <iostream>

class polynomial 
{

private:
    std::deque<double> _coefs;
    
public:
    // TODO: reverse constructor
    polynomial();
    polynomial(std::vector<double> const &coefs);

public:
    void print(std::ostream &stream = std::cout, size_t precision = 3) const;
    void println(std::ostream &stream = std::cout, size_t precision = 3) const;

public:
    bool equals(polynomial const &other, double eps) const;

public:
    polynomial &operator+=(polynomial const &other);
    polynomial &operator-=(polynomial const &other);
    polynomial &operator*=(double term);
    polynomial &operator/=(double term);
    polynomial &operator*=(polynomial const &other);
    
    polynomial operator+(polynomial const &other) const;
    polynomial operator-(polynomial const &other) const;
    polynomial operator*(double term) const;
    polynomial operator/(double term) const;
    polynomial operator*(polynomial const &other) const;
    
    double &operator[](size_t idx);
    double operator[](size_t idx) const;
    
    double operator()(double x) const;

public:
    size_t degree() const;
    double valueAt(double x) const;

private:
    void extend_to_size_of(polynomial const &other);
    void shrink();

};

#endif // _NUM_METHODS_UTIL_POLYNOMIAL_H_