#include <stdexcept>
#include <cmath>
#include <limits>
#include <iomanip>

#include <checks.h>

#include "../include/piecewise_polynomial.h"

piecewise_polynomial::piecewise_polynomial():
        _bounds(1, 0),
        _polynomials(0),
        _degree(0)
{ }

piecewise_polynomial::piecewise_polynomial(
    double left_bound, double right_bound, polynomial const &poly):
        _bounds({left_bound, right_bound}),
        _polynomials(1, poly),
        _degree(poly.degree())
{ }


void piecewise_polynomial::print(std::ostream &stream, size_t precision, size_t dot_precision)
{
    std::ios_base::fmtflags former_flags = stream.flags();
    std::streamsize former_width = stream.width();
    std::streamsize former_precision = stream.precision();
    size_t bound_max_len = 0;
    double eps = std::pow(10.0, -static_cast<long long>(precision));
    
    for (size_t i = 0; i < _bounds.size(); ++i)
    {
        bound_max_len = std::max(bound_max_len, std::to_string(_bounds[i]).find('.') +
                (_bounds[i] < 0 ? 1 : 0));
    }
    
    bound_max_len += dot_precision;
    stream.precision(dot_precision);
    stream.setf(std::ios_base::fixed);
    
    for (size_t i = 0; i < _polynomials.size(); ++i)
    {
        if (i > 0)
        {
            stream << std::endl;
        }
        
        stream << std::setw(bound_max_len) << _bounds[i] << " < x <= " <<
                std::setw(bound_max_len) << _bounds[i+1] << " : ";
        
        if (_polynomials[i][_polynomials[i].degree()] > 0)
        {
            stream << ' ';
        }
        
        _polynomials[i].print(stream, precision);
    }
    
    stream.setf(former_flags);
    stream.width(former_width);
    stream.precision(former_precision);
}

void piecewise_polynomial::println(std::ostream &stream, size_t precision, size_t dot_precision)
{
    print(stream, precision, dot_precision);
    stream << std::endl;
}


double piecewise_polynomial::operator()(double x) const
{
    return valueAt(x);
}

double piecewise_polynomial::operator[](size_t idx) const
{
    return _bounds[idx];
}


piecewise_polynomial &piecewise_polynomial::push_front(double bound, polynomial const &poly)
{
    if (_bounds.front() - bound < std::numeric_limits<double>::epsilon())
    {
        throw std::invalid_argument("The bound included into bounds interval");
    }
    
    _bounds.push_front(bound);
    _polynomials.push_front(poly);
    _degree = std::max(_degree, poly.degree());
    
    return *this;
}

piecewise_polynomial &piecewise_polynomial::push_back(double bound, polynomial const &poly)
{
    if (bound - _bounds.back() < std::numeric_limits<double>::epsilon())
    {
        throw std::invalid_argument("Violation of adjacent bound constraint");
    }
    
    _bounds.push_back(bound);
    _polynomials.push_back(poly);
    _degree = std::max(_degree, poly.degree());
    
    return *this;
}

double piecewise_polynomial::get_bound(size_t idx) const
{
    if (idx >= _bounds.size())
    {
        throw std::range_error("Index is out of bounds");
    }
    
    return _bounds[idx];
}

void piecewise_polynomial::set_bound(size_t idx, double x)
{
    if (idx >= _bounds.size())
    {
        throw std::range_error("Index is out of bounds");
    }
    
    double eps = std::numeric_limits<double>::epsilon();
    
    if ((idx > 0 && (_bounds[idx-1] - x >= -eps)) ||
            (idx + 1 < _bounds.size() && (x - _bounds[idx+1] >= -eps)))
    {
        throw std::invalid_argument("Violation of adjacent bound constraint");
    }
    
    _bounds[idx] = x;
}


size_t piecewise_polynomial::bound_count() const
{
    return _bounds.size();
}

size_t piecewise_polynomial::degree() const
{
    return _degree;
}

double piecewise_polynomial::valueAt(double x) const
{
    size_t bound_idx = std::lower_bound(_bounds.begin(), _bounds.end(), x) - _bounds.begin();
    
    if (bound_idx == 0 || bound_idx == _bounds.size())
    {
        throw std::invalid_argument("Point out of the domain");
    }
    
    return _polynomials[bound_idx - 1].valueAt(x);
}