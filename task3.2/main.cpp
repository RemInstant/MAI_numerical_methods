#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <cmath>
#include <numbers>
#include <format>
#include <algorithm>

#include <checks.h>
#include <polynomial.h>
#include <piecewise_polynomial.h>
#include <algorithms.h>

// int get_interpolation_accuracy(
//     polynomial interpolation,
//     std::vector<double> points,
//     std::vector<double> values)
// {
//     if (points.size() != values.size())
//     {
//         throw std::invalid_argument("Point and values count are not equal");
//     }
    
//     double eps = 1.0;
//     size_t accuracy = 0;
    
//     for (; eps > std::numeric_limits<double>::epsilon(); ++accuracy, eps /= 10)
//     {
//         for (size_t i = 0; i < points.size(); ++i)
//         {
//             if (std::abs(interpolation.valueAt(points[i]) - values[i]) > eps)
//             {
//                 return accuracy - 1;
//             }
//         }
//     }
    
//     return accuracy;
// }


int main()
{
    std::freopen("input.txt", "r", stdin);
    
    std::vector<double> points = {0,1,2,3,4};
    std::vector<double> values = {0,1.8415,2.9093,3.1411,3.2432};
    
    piecewise_polynomial spline = algorithms::build_spline(points, values);
    
    spline.println(std::cout, 6);
}