// C++ Includes
#include <cmath>

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
  long double pi = acos((long double)-1.0);
  switch (CASE)
  {
    case 1: return p(0)*(1-p(0))*p(1)*(1-p(1))*p(2)*(1-p(2)) + 1;
    default: return p(0)*(1-p(0))*p(1)*(1-p(1))*p(2)*(1-p(2)) + 1;
  }
}

Gradient initial_dvalue_u (const Point & p,
                    const Parameters & parameters,
                    const std::string &,
                    const std::string &)
{
  long double pi = acos((long double)-1.0);
  switch (CASE)
  {
    case 1: return Gradient((1-2*p(0))*p(1)*(1-p(1))*p(2)*(1-p(2)), p(0)*(1-p(0))*(1-2*p(1))*p(2)*(1-p(2)), p(0)*(1-p(0))*p(1)*(1-p(1))*(1-2*p(2)));
    default: return Gradient((1-2*p(0))*p(1)*(1-p(1))*p(2)*(1-p(2)), p(0)*(1-p(0))*(1-2*p(1))*p(2)*(1-p(2)), p(0)*(1-p(0))*p(1)*(1-p(1))*(1-2*p(2)));
  }
}

Number initial_value_u (const Point & p,
                    const Parameters & parameters)
{
  long double pi = acos((long double)-1.0);
  switch (CASE)
  {
    case 1: return p(0)*(1-p(0))*p(1)*(1-p(1))*p(2)*(1-p(2)) + 1;
    default: return p(0)*(1-p(0))*p(1)*(1-p(1))*p(2)*(1-p(2)) + 1;
  }
}

Gradient initial_dvalue_u (const Point & p,
                    const Parameters & parameters)
{
  long double pi = acos((long double)-1.0);
  switch (CASE)
  {
    case 1: return Gradient((1-2*p(0))*p(1)*(1-p(1))*p(2)*(1-p(2)), p(0)*(1-p(0))*(1-2*p(1))*p(2)*(1-p(2)), p(0)*(1-p(0))*p(1)*(1-p(1))*(1-2*p(2)));
    default: return Gradient((1-2*p(0))*p(1)*(1-p(1))*p(2)*(1-p(2)), p(0)*(1-p(0))*(1-2*p(1))*p(2)*(1-p(2)), p(0)*(1-p(0))*p(1)*(1-p(1))*(1-2*p(2)));
  }
}

Number boundary_value_u (const Point & p, const Real time = 0)
{
  switch (CASE)
  {
    case 1: return p(0)*(1-p(0))*p(1)*(1-p(1))*p(2)*(1-p(2)) + 1;
    default: return p(0)*(1-p(0))*p(1)*(1-p(1))*p(2)*(1-p(2)) + 1;
  }
}

Number rhs_value_f (const Point & p, const Real time = 0)
{
    switch (CASE)
    {
      case 1: return 2*p(1)*(1-p(1))*p(2)*(1-p(2)) + p(0)*(1-p(0))*(2)*p(2)*(1-p(2)) + p(0)*(1-p(0))*p(1)*(1-p(1))*(2);
      default: return 2*p(1)*(1-p(1))*p(2)*(1-p(2)) + p(0)*(1-p(0))*(2)*p(2)*(1-p(2)) + p(0)*(1-p(0))*p(1)*(1-p(1))*(2);
    }
}
