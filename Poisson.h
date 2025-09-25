// The nonlinear solver and system we will be using
#include "libmesh/nonlinear_solver.h"
#include "libmesh/nonlinear_implicit_system.h"
// Bring in everything from the libMesh namespace
using namespace libMesh;

class Poisson : public NonlinearImplicitSystem::ComputeResidual,
                                   public NonlinearImplicitSystem::ComputeJacobian
{
private:
  EquationSystems & es;
  std::string sys_name, var_name; 
  Number (* boundary_value) (const Point & p, const Real time);
  Number (* rhs_value) (const Point & p, const Real time);
public:

  Poisson (EquationSystems & es_in, std::string sys_name_in, std::string var_name_in, Number boundary_value_in (const Point & p, const Real time), Number rhs_value_in (const Point & p, const Real time)) :
    es(es_in), sys_name(sys_name_in), var_name(var_name_in), boundary_value(boundary_value_in), rhs_value(rhs_value_in)
  {}

  /**
   * Evaluate the Jacobian of the nonlinear system.
   */
  virtual void jacobian (const NumericVector<Number> & soln,
                         SparseMatrix<Number> & jacobian,
                         NonlinearImplicitSystem & /*sys*/);
  

  /**
   * Evaluate the residual of the nonlinear system.
   */
  virtual void residual (const NumericVector<Number> & soln,
                         NumericVector<Number> & residual,
                         NonlinearImplicitSystem & /*sys*/);
  


};
