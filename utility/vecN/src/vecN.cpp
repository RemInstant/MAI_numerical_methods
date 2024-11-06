#include <stdexcept>
#include <cmath>
#include <iomanip>

#include <checks.h>

#include "../include/vecN.h"

vecN::vecN():
        _data()
{ }

vecN::vecN(size_t dim):
        _data(std::vector<double>(dim, 0))
{ }

vecN::vecN(size_t dim, double value):
        _data(std::vector<double>(dim, value))
{ }

vecN::vecN(std::vector<double> const &data):
        _data(data)
{ }

// vecN::vecN(std::vector<double> data):
//         _data(std::move(data))
// { }


void vecN::print(std::ostream &stream, size_t precision) const
{
    std::ios_base::fmtflags former_flags = stream.flags();
    std::streamsize former_precision = stream.precision();
    std::streamsize former_width = stream.width();
    
    stream.precision(precision);
    stream.setf(std::ios_base::fixed);
    
    for (size_t i = 0; i < _data.size(); ++i)
    {
        stream << _data[i] << ' ';
    }
    
    stream.setf(former_flags);
    stream.precision(former_precision);
    stream.width(former_width);
}

void vecN::println(std::ostream &stream, size_t precision) const
{
    print(stream, precision);
    stream << std::endl;
}


bool vecN::equals(vecN const &other, double eps) const
{
    throw_if_other_dim(other);
    checks::throw_if_invalid_eps(eps);
    
    for (size_t i = 0; i < _data.size(); ++i)
    {
        if (std::abs(_data[i] - other._data[i]) > eps)
        {
            return false;
        }
    }
    
    return true;
}


vecN &vecN::operator+=(vecN const &other)
{
    throw_if_other_dim(other);
    
    for (size_t i = 0; i < _data.size(); ++i)
    {
        _data[i] += other._data[i];
    }
    
    return *this;
}

vecN &vecN::operator-=(vecN const &other)
{
    throw_if_other_dim(other);
    
    for (size_t i = 0; i < _data.size(); ++i)
    {
        _data[i] -= other._data[i];
    }
    
    return *this;
}

vecN &vecN::operator*=(double term)
{
    for (size_t i = 0; i < _data.size(); ++i)
    {
        _data[i] *= term;
    }
    
    return *this;
}

vecN &vecN::operator/=(double term)
{
    checks::throw_if_zero_divisor(term);
    
    for (size_t i = 0; i < _data.size(); ++i)
    {
        _data[i] /= term;
    }
    
    return *this;
}


vecN vecN::operator+(vecN const &other) const
{
    return vecN(*this) += other;
}

vecN vecN::operator-(vecN const &other) const
{
    return vecN(*this) -= other;
}

vecN vecN::operator*(double term) const
{
    return vecN(*this) *= term;
}

vecN vecN::operator/(double term) const
{
    return vecN(*this) /= term;
}


vecN operator*(double term, vecN const &vec)
{
    return vec * term;
}


double vecN::scalar_prod(vecN const &lhs, vecN const &rhs)
{
    if (lhs.size() != rhs.size())
    {
        throw std::invalid_argument("Dimensions does not correspond");
    }
    
    double ans = 0;
    
    for (size_t i = 0; i < lhs.size(); ++i)
    {
        ans += lhs[i] * rhs[i];
    }
    
    return ans;
}


double &vecN::operator[](size_t idx)
{
    return _data[idx];
}

double vecN::operator[](size_t idx) const
{
    return _data[idx];
}


size_t vecN::size() const
{
    return _data.size();
}


double vecN::l1_norm() const
{
    double norm = 0;
    
    for (size_t i = 0; i < _data.size(); ++i)
    {
        norm += _data[i];
    }
    
    return norm;
}

double vecN::l2_norm() const
{
    double norm = 0;
    
    for (size_t i = 0; i < _data.size(); ++i)
    {
        norm += _data[i] * _data[i];
    }
    
    return std::sqrt(norm);
}

double vecN::continuous_norm() const
{
    double norm = 0;
    
    for (size_t i = 0; i < _data.size(); ++i)
    {
        norm = std::max(norm, std::abs(_data[i]));
    }
    
    return norm;
}

void vecN::throw_if_other_dim(vecN const &other) const
{
    if (_data.size() != other.size())
    {
        throw std::invalid_argument("Dimensions does not correspond");
    }
}


