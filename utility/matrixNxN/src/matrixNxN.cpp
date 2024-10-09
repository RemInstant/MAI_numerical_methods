#include <stdexcept>
#include <cmath>
#include <limits>

#include "../include/matrixNxN.h"

matrixNxN::matrixNxN(size_t dim):
    _data(std::vector<vecN>(dim, vecN(dim)))
{ }

matrixNxN::matrixNxN(size_t dim, double value):
    _data(std::vector<vecN>(dim, vecN(dim, value)))
{ }

matrixNxN matrixNxN::identical(size_t dim)
{
    matrixNxN identical(dim, 0);
    
    for (size_t i = 0; i < dim; ++i)
    {
        identical[i][i] = 1;
    }
    
    return identical;
}

bool matrixNxN::equals(matrixNxN const &other, double eps) const
{
    throw_if_other_dim(other);
    
    for (size_t i = 0; i < _data.size(); ++i)
    {
        if (!_data[i].equals(other._data[i], eps))
        {
            return false;
        }
    }
    
    return true;
}


matrixNxN &matrixNxN::operator+=(matrixNxN const &other)
{
    throw_if_other_dim(other);
    
    for (size_t i = 0; i < _data.size(); ++i)
    {
        _data[i] += other._data[i];
    }
    
    return *this;
}

matrixNxN &matrixNxN::operator-=(matrixNxN const &other)
{
    throw_if_other_dim(other);
    
    for (size_t i = 0; i < _data.size(); ++i)
    {
        _data[i] -= other._data[i];
    }
    
    return *this;
}

matrixNxN &matrixNxN::operator*=(double term)
{
    for (size_t i = 0; i < _data.size(); ++i)
    {
        _data[i] *= term;
    }
    
    return *this;
}

matrixNxN &matrixNxN::operator/=(double term)
{
    for (size_t i = 0; i < _data.size(); ++i)
    {
        _data[i] /= term;
    }
    
    return *this;
}

matrixNxN &matrixNxN::operator*=(matrixNxN const &other)
{
    throw_if_other_dim(other);
    
    vecN tmp(_data.size());
    
    for (size_t i = 0; i < _data.size(); ++i)
    {
        for (size_t j = 0; j < _data.size(); ++j)
        {
            double cell = 0;
            
            for (size_t k = 0; k < _data.size(); ++k)
            {
                cell += _data[i][k] * other[k][j];
            }
            
            tmp[j] = cell;
        }
        
        _data[i] = tmp;
    }
    
    return *this;
}


matrixNxN matrixNxN::operator+(matrixNxN const &other) const
{
    return matrixNxN(*this) += other;
}

matrixNxN matrixNxN::operator-(matrixNxN const &other) const
{
    return matrixNxN(*this) -= other;
}

matrixNxN matrixNxN::operator*(double term) const
{
    return matrixNxN(*this) *= term;
}

matrixNxN matrixNxN::operator/(double term) const
{
    return matrixNxN(*this) /= term;
}

matrixNxN matrixNxN::operator*(matrixNxN const& other) const
{
    return matrixNxN(*this) *= other;
}

vecN operator*(vecN const &vec, matrixNxN const &matrix)
{
    matrix.throw_if_other_dim(vec);
    
    vecN res(matrix.size());
    
    for (size_t j = 0; j < matrix.size(); ++j)
    {
        res[j] = 0;
        
        for (size_t k = 0; k < matrix.size(); ++k)
        {
            res[j] += vec[k] * matrix[k][j];
        }
    }
    
    return res;
}

vecN operator*(matrixNxN const &matrix, vecN const &vec)
{
    matrix.throw_if_other_dim(vec);
    
    vecN res(matrix.size());
    
    for (size_t i = 0; i < matrix.size(); ++i)
    {
        res[i] = 0;
        
        for (size_t k = 0; k < matrix.size(); ++k)
        {
            res[i] += matrix[i][k] * vec[k];
        }
    }
    
    return res;
}

vecN &matrixNxN::operator[](size_t idx)
{
    return _data[idx];
}

vecN matrixNxN::operator[](size_t idx) const
{
    return _data[idx];
}


size_t matrixNxN::size() const
{
    return _data.size();
}


double matrixNxN::l1_norm() const
{
    double norm = 0;
    
    for (size_t j = 0; j < _data.size(); ++j)
    {
        double tmp = 0;
        
        for (size_t i = 0; i < _data.size(); ++i)
        {
            tmp += std::abs(_data[i][j]);
        }
        
        norm = std::max(norm, tmp);
    }
    
    return norm;
}

double matrixNxN::continuous_norm() const
{
    double norm = 0;
    
    for (size_t i = 0; i < _data.size(); ++i)
    {
        double tmp = 0;
        
        for (size_t j = 0; j < _data.size(); ++j)
        {
            tmp += std::abs(_data[i][j]);
        }
        
        norm = std::max(norm, tmp);
    }
    
    return norm;
}


matrixNxN matrixNxN::inverse() const
{
    matrixNxN coefs = *this;
    matrixNxN constant_terms = identical(_data.size());
    matrixNxN inversed(_data.size());
    size_t swap_count = 0;
    
    // First phase of gauss elimination
    for (size_t i = 0; i < _data.size(); ++i)
    {
        // Search of leading element
        size_t idx_to_swap = i;
        
        for (size_t j = idx_to_swap + 1; j < _data.size(); ++j)
        {
            if (std::abs(coefs[j][i]) > std::abs(coefs[idx_to_swap][i]))
            {
                idx_to_swap = j;
            }
        }
        
        swap_count += (idx_to_swap > i) ? (idx_to_swap - i) : (i - idx_to_swap);
        std::swap(coefs[i], coefs[idx_to_swap]);
        std::swap(constant_terms[i], constant_terms[idx_to_swap]);
        
        if (std::abs(coefs[i][i]) < std::numeric_limits<double>::epsilon())
        {
            throw std::invalid_argument("Matrix cannot be inversed");
        }
        
        for (size_t j = i+1; j < _data.size(); ++j)
        {
            if (coefs[j][i] < std::numeric_limits<double>::epsilon())
            {
                continue;
            }
            
            double mult = coefs[j][i] / coefs[i][i];
            
            for (size_t k = i; k < _data.size(); ++k)
            {
                if (k >= i)
                {
                    coefs[j][k] -= coefs[i][k] * mult;
                }
                constant_terms[j][k] -= coefs[i][k] * mult;
            }
        }
    }
    
    // Second phase of gauss elimination
    for (size_t k = 0; k < _data.size(); ++k)
    {
        for (size_t i = _data.size(); i > 0; --i)
        {
            inversed[i-1][k] = constant_terms[i-1][k];
            
            for (size_t j = i; j < _data.size(); ++j)
            {
                inversed[i-1][k] -= coefs[i-1][j] * inversed[j][k];
            }
            
            inversed[i-1][k] /= coefs[i-1][i-1];
        }
        
    }
    
    return inversed;
}


void matrixNxN::throw_if_other_dim(vecN const &other) const
{
    if (_data.size() != other.size())
    {
        throw std::invalid_argument("Dimensions of vectors does not correspond");
    }
}

void matrixNxN::throw_if_other_dim(matrixNxN const &other) const
{
    if (_data.size() != other.size())
    {
        throw std::invalid_argument("Dimensions of vectors does not correspond");
    }
}


