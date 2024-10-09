#ifndef _NUM_METHODS_UTIL_MATRIXNXN_
#define _NUM_METHODS_UTIL_MATRIXNXN_

#include <vector>
#include <vecN.h>

class matrixNxN 
{

private:
    std::vector<vecN> _data;

public:
    matrixNxN(size_t dim);
    matrixNxN(size_t dim, double value);

public:
    static matrixNxN identical(size_t dim);

public:
    bool equals(matrixNxN const &other, double eps) const;

public:
    matrixNxN &operator+=(matrixNxN const &other);
    matrixNxN &operator-=(matrixNxN const &other);
    matrixNxN &operator*=(double term);
    matrixNxN &operator/=(double term);
    matrixNxN &operator*=(matrixNxN const &other);
    
    matrixNxN operator+(matrixNxN const &other) const;
    matrixNxN operator-(matrixNxN const &other) const;
    matrixNxN operator*(double term) const;
    matrixNxN operator/(double term) const;
    matrixNxN operator*(matrixNxN const &other) const;
    
    friend vecN operator*(vecN const &vec, matrixNxN const &matrix); 
    friend vecN operator*(matrixNxN const &matrix, vecN const &vec); 
    
    vecN &operator[](size_t idx);
    vecN operator[](size_t idx) const;

public:
    size_t size() const;
    
    double l1_norm() const;
    double continuous_norm() const;
    
    matrixNxN inverse() const;

private:
    void throw_if_other_dim(vecN const &other) const;
    void throw_if_other_dim(matrixNxN const &other) const;

};

#endif // _NUM_METHODS_UTIL_MATRIXNXN_