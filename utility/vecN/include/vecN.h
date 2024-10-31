#ifndef _NUM_METHODS_UTIL_VECN_
#define _NUM_METHODS_UTIL_VECN_

#include <vector>

class vecN 
{
    
private:
    std::vector<double> _data;

public:
    vecN();
    vecN(size_t dim);
    vecN(size_t dim, double value);
    vecN(std::vector<double> const &data);

public:
    void print(std::ostream &stream, size_t precision);
    void println(std::ostream &stream, size_t precision);

public:
    bool equals(vecN const &other, double eps) const;

public:
    vecN &operator+=(vecN const &other);
    vecN &operator-=(vecN const &other);
    vecN &operator*=(double term);
    vecN &operator/=(double term);
    
    vecN operator+(vecN const &other) const;
    vecN operator-(vecN const &other) const;
    vecN operator*(double term) const;
    vecN operator/(double term) const;
    
    friend vecN operator*(double term, vecN const &vec);
    
    static double scalar_prod(vecN const &lhs, vecN const &rhs);
    
    double &operator[](size_t idx);
    double operator[](size_t idx) const;

public:
    size_t size() const;
    
    double l1_norm() const;
    double l2_norm() const;
    double continuous_norm() const;

private:
    void throw_if_other_dim(vecN const &other) const;

};

#endif // _NUM_METHODS_UTIL_VECN_