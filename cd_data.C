// C++ Includes
#include <cmath>
#include <complex>

// Mesh library includes
#include "libmesh/libmesh_common.h"
#include "libmesh/equation_systems.h"
#include "../operator/play.h"
#include "cd_data.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

extern int CASE;

Number initial_value_u (const Point & p,
                    const Parameters & parameters,
                    const std::string &,
                    const std::string &)
{
  std::complex<double> z(p(0), p(1));
  std::complex<double> root = std::sqrt(z); // 复数平方根
  return root.real(); // 取实部  
}

Gradient initial_dvalue_u (const Point & p,
                    const Parameters & parameters,
                    const std::string &,
                    const std::string &)
{
    double x = p(0), y = p(1);
    const double eps = 1e-14;
    double r = sqrt(x*x + y*y);
    if (r < eps) {
        // near origin: limit behavior
        // v -> 0, gradient magnitude ~ (1/2) / sqrt(r) -> infinite,
        // but we return finite safe values (could also throw or set nan)
        return {0.0, 0.0, 0.0};
    }

    double theta = atan2(y, x); // (-pi, pi]
    double half = 0.5 * theta;
    double cos_half = cos(half);
    double sin_half = sin(half);

    double v = sqrt(r) * cos_half;
    double coeff = 0.5 / sqrt(r);
    double dvx = coeff * cos_half;
    double dvy = coeff * sin_half;

    return Gradient(dvx,dvy);
}

Number initial_value_u (const Point & p,
                    const Parameters & parameters)
{
  std::complex<double> z(p(0), p(1));
  std::complex<double> root = std::sqrt(z); // 复数平方根
  return root.real(); // 取实部  
}

Gradient initial_dvalue_u (const Point & p,
                    const Parameters & parameters)
{
    double x = p(0), y = p(1);
    const double eps = 1e-14;
    double r = sqrt(x*x + y*y);
    if (r < eps) {
        // near origin: limit behavior
        // v -> 0, gradient magnitude ~ (1/2) / sqrt(r) -> infinite,
        // but we return finite safe values (could also throw or set nan)
        return {0.0, 0.0, 0.0};
    }

    double theta = atan2(y, x); // (-pi, pi]
    double half = 0.5 * theta;
    double cos_half = cos(half);
    double sin_half = sin(half);

    double v = sqrt(r) * cos_half;
    double coeff = 0.5 / sqrt(r);
    double dvx = coeff * cos_half;
    double dvy = coeff * sin_half;

    return Gradient(dvx,dvy);
}

Number boundary_value_u (const Point & p, const Real time = 0)
{
  std::complex<double> z(p(0), p(1));
  std::complex<double> root = std::sqrt(z); // 复数平方根
  return root.real(); // 取实部  
}

Number rhs_value_f (const Point & p, const Real time = 0)
{
    return 0;
}
