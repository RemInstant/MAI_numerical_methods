#include <stdexcept>
#include <cmath>
#include <limits>
#include <iomanip>

#include <checks.h>

#include "../include/polynomial.h"

polynomial::polynomial():
    _coefs(1, 0)
{ }


polynomial::polynomial(std::vector<double> const &coefs):
    _coefs(coefs.size())
{
    if (coefs.empty())
    {
        _coefs.push_back(0);
    }
    else
    {
        for (size_t i = 0; i < coefs.size(); ++i)
        {
            _coefs[i] = coefs[i];
        }
    }
}


void polynomial::print(std::ostream &stream, size_t precision)
{
    std::ios_base::fmtflags former_flags = stream.flags();
    std::streamsize former_precision = stream.precision();
    
    stream.precision(precision);
    stream.setf(std::ios_base::fixed);
    
    double eps = std::pow(10.0, -static_cast<long long>(precision));
    
    std::string positive_start = "";
    std::string negative_start = "-";
    std::string spacer = "";
    
    for (size_t i = _coefs.size(); i > 0; --i)
    {
        if (_coefs.size() > 1 && std::abs(_coefs[i-1]) < eps)
        {
            continue;
        }
        
        stream << (_coefs[i-1] >= 0 ? positive_start : negative_start) << spacer;
        
        if (i == 1 || std::abs(_coefs[i-1] - 1) > eps)
        {
            stream << std::abs(_coefs[i-1]);
        }
        
        switch (i)
        {
        case 1:
            break;
        case 2:
            std::cout << 'x';
            break;
        default: 
            std::cout << "x^" << (i-1);
            break;
        }
        
        positive_start = " +";
        negative_start = " -";
        spacer = " ";
    }
    
    stream.setf(former_flags);
    stream.precision(former_precision);
}

void polynomial::println(std::ostream &stream, size_t precision)
{
    print(stream, precision);
    stream << std::endl;
}


bool polynomial::equals(polynomial const &other, double eps) const
{
    size_t mn = std::min(_coefs.size(), other._coefs.size());
    
    for (size_t i = 0; i < mn; ++i)
    {
        if (std::abs(_coefs[i] - other[i]) > eps)
        {
            return false;
        }
    }
    for (size_t i = mn; _coefs.size(); ++i)
    {
        if (std::abs(_coefs[i]) > eps)
        {
            return false;
        }
    }
    for (size_t i = mn; other._coefs.size(); ++i)
    {
        if (std::abs(other._coefs[i]) > eps)
        {
            return false;
        }
    }
    
    return true;
}


polynomial &polynomial::operator+=(polynomial const &other)
{
    extend_to_size_of(other);
    
    for (size_t i = 0; i < _coefs.size(); ++i)
    {
        _coefs[i] += other[i];
    }
    
    shrink();
    
    return *this;
}

polynomial &polynomial::operator-=(polynomial const &other)
{
    extend_to_size_of(other);
    
    for (size_t i = 0; i < _coefs.size(); ++i)
    {
        _coefs[i] -= other[i];
    }
    
    shrink();
    
    return *this;
}

polynomial &polynomial::operator*=(double term)
{
    for (size_t i = 0; i < _coefs.size(); ++i)
    {
        _coefs[i] *= term;
    }
    
    shrink();
    
    return *this;
}

polynomial &polynomial::operator/=(double term)
{
    checks::throw_if_zero_divisor(term);
    
    for (size_t i = 0; i < _coefs.size(); ++i)
    {
        _coefs[i] /= term;
    }
    
    return *this;
}

polynomial &polynomial::operator*=(polynomial const &other)
{
    std::deque tmp(_coefs);
    
    for (size_t i = 0; i < _coefs.size(); ++i)
    {
        _coefs[i] *= other._coefs.back();
    }
    
    for (size_t i = other._coefs.size() - 1; i > 0; --i)
    {
        double coef = other._coefs[i-1];
        _coefs.push_front(0);
        
        for (size_t j = 0; j < tmp.size(); ++j)
        {
            _coefs[j] += tmp[j] * coef;
        }
    }
    
    shrink();
    
    return *this;
}

polynomial polynomial::operator+(polynomial const &other) const
{
    return polynomial(*this) += other;
}

polynomial polynomial::operator-(polynomial const &other) const
{
    return polynomial(*this) -= other;
}

polynomial polynomial::operator*(double term) const
{
    return polynomial(*this) *= term;
}

polynomial polynomial::operator/(double term) const
{
    return polynomial(*this) /= term;
}

polynomial polynomial::operator*(polynomial const &other) const
{
    return polynomial(*this) *= other;
}

double &polynomial::operator[](size_t idx)
{
    return _coefs[idx];
}

double polynomial::operator[](size_t idx) const
{
    return _coefs[idx];
}

double polynomial::operator()(double x) const
{
    return valueAt(x);
}



size_t polynomial::degree() const
{
    return _coefs.size() - 1;
}

double polynomial::valueAt(double x) const
{
    double ans = 0;
    
    for (size_t i = _coefs.size(); i > 0; --i)
    {
        ans *= x;
        ans += _coefs[i-1];
    }
    
    return ans;
}


void polynomial::extend_to_size_of(polynomial const &other)
{
    if (other._coefs.size() > _coefs.size())
    {
        _coefs.resize(other._coefs.size(), 0);
    }
}

void polynomial::shrink()
{
    while (_coefs.size() > 1 && std::abs(_coefs.back()) < std::numeric_limits<double>::epsilon())
    {
        _coefs.pop_back();
    }
    
    _coefs.shrink_to_fit();
}