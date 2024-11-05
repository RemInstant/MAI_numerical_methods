#ifndef _NUM_METHODS_UTIL_POLYNOME_H_
#define _NUM_METHODS_UTIL_POLYNOME_H_

#include <deque>
#include <vector>
#include <iostream>

class polynome 
{

private:
    std::deque<double> _coefs;
    
public:
    polynome();
    polynome(std::vector<double> const &coefs);

public:
    void print(std::ostream &stream = std::cout, size_t precision = 3);
    void println(std::ostream &stream = std::cout, size_t precision = 3);

public:
    bool equals(polynome const &other, double eps) const;

public:
    polynome &operator+=(polynome const &other);
    polynome &operator-=(polynome const &other);
    polynome &operator*=(double term);
    polynome &operator/=(double term);
    polynome &operator*=(polynome const &other);
    
    polynome operator+(polynome const &other) const;
    polynome operator-(polynome const &other) const;
    polynome operator*(double term) const;
    polynome operator/(double term) const;
    polynome operator*(polynome const &other) const;
    
    double &operator[](size_t idx);
    double operator[](size_t idx) const;
    
    double operator()(double x) const;

public:
    size_t degree() const;
    
    double valueAt(double x) const;

private:
    void extend_to_size_of(polynome const &other);
    void shrink();

};

#endif // _NUM_METHODS_UTIL_POLYNOME_H_