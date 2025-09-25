// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



// <h1>Transient Example 1 - Solving a Transient Linear System in Parallel</h1>
// \author Benjamin S. Kirk
// \date 2003
//
// This example shows how a simple, linear transient
// system can be solved in parallel.  The system is simple
// scalar Parabolic with a specified external
// velocity.  The initial condition is given, and the
// solution is advanced in time with a standard implicit Euler
// time-stepping strategy.

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <sstream>
#include <math.h>
#include <fstream>

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/gmv_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/gnuplot_io.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/mesh_generation.h"

// This example will solve a linear transient system,
// so we need to include the TransientLinearImplicitSystem definition.
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/vector_value.h"

// The nonlinear solver and system we will be using
#include "libmesh/nonlinear_solver.h"
#include "libmesh/nonlinear_implicit_system.h"

// The definition of a geometric element
#include "libmesh/elem.h"

// To impose Dirichlet boundary conditions
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/analytic_function.h"
#include "libmesh/zero_function.h"
#include "libmesh/boundary_info.h"

#include "libmesh/getpot.h"
#include "libmesh/exact_solution.h"
#include "libmesh/enum_xdr_mode.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

#include "cd_data.h"
#include "Poisson.h"


// #define REFERENCE_GEN

// #define ERROR_CALC

int CASE = 1;
int mesh_partition = 40; 
         
// We can now begin the main program.  Note that this
// example will fail if you are using complex numbers
// since it was designed to be run only with real numbers.
int main (int argc, char ** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);

  // This example requires a linear solver package.
  libmesh_example_requires(libMesh::default_solver_package() != INVALID_SOLVER_PACKAGE,
                           "--enable-petsc, --enable-trilinos, or --enable-eigen");


  GetPot command_line (argc, argv);

  // Brief message to the user regarding the program name
  // and command line arguments.
  libMesh::out << "Running " << argv[0];

  for (int i=1; i<argc; i++)
    libMesh::out << " " << argv[i];

  libMesh::out << std::endl << std::endl;

  if (command_line.search(1, "-case"))
    CASE = command_line.next(CASE);

  if (command_line.search(1, "-N"))
    mesh_partition = command_line.next(mesh_partition);

  int dim = 1;
  if (command_line.search(1, "-dim"))
    dim = command_line.next(dim);

  int visual = 0;
  if (command_line.search(1, "-visual"))
      visual = command_line.next(visual);

  // This example requires Adaptive Mesh Refinement support - although
  // it only refines uniformly, the refinement code used is the same
  // underneath
#ifndef LIBMESH_ENABLE_AMR
  libmesh_example_requires(false, "--enable-amr");
#else

  // Skip this 2D example if libMesh was compiled as 1D-only.
  libmesh_example_requires(2 <= LIBMESH_DIM, "2D support");

  // Read the mesh from file.  This is the coarse mesh that will be used
  // in example 10 to demonstrate adaptive mesh refinement.  Here we will
  // simply read it in and uniformly refine it 5 times before we compute
  // with it.
  //
  // Create a mesh object, with dimension to be overridden later,
  // distributed across the default MPI communicator.
  Mesh mesh(init.comm());

  // Create a mesh with user-defined dimension.
  // Read number of elements from command line

  Real halfwidth = dim > 1 ? 1. : 0.;
  Real halfheight = dim > 2 ? 1. : 0.;
  MeshTools::Generation::build_cube (mesh,
                                      mesh_partition,
                                      (dim>1) ? mesh_partition : 0,
                                      (dim>2) ? mesh_partition : 0,
                                      0., 1.,
                                      0., halfwidth,
                                      0., halfheight,
                                      (dim==1)    ? EDGE2 :
                                      ((dim == 2) ? TRI3 : TET4));  
  // Print information about the mesh to the screen.
  mesh.print_info();
  // mesh.get_boundary_info().print_info();


  // A pretty update message
  libMesh::out << " System solution step " <<  "..." << std::endl;

  EquationSystems elliptic_systems (mesh);

  NonlinearImplicitSystem & system =
    elliptic_systems.add_system<NonlinearImplicitSystem> ("elliptic");

  unsigned int u_var = system.add_variable ("u", FIRST);


  Poisson lde(elliptic_systems, "elliptic", "u", boundary_value_u, rhs_value_f);
  system.nonlinear_solver->residual_object = &lde;
  system.nonlinear_solver->jacobian_object = &lde;

  elliptic_systems.init();

  system.solve();
  
  
  if (visual) {
  // If Exodus is available, we'll write all timesteps to the same file
  // rather than one file per timestep.
    std::string exodus_filename = "Dim_" + std::to_string(dim) + "_CASE_" + std::to_string(CASE) + "_N_" + std::to_string(mesh_partition) + ".e";
    ExodusII_IO(mesh).write_equation_systems (exodus_filename, elliptic_systems);
  }

  libMesh::out << " System solution finished " <<  "..." << std::endl << std::endl;

  ExactSolution comparison(elliptic_systems);
  comparison.attach_exact_value(initial_value_u);
  comparison.attach_exact_deriv(initial_dvalue_u);
  Real error;

  comparison.compute_error("elliptic", "u");
  
  error = comparison.l2_error("elliptic", "u");
  libMesh::out << "L2 error of u w.r.t reference: " << error << std::endl << std::endl;
 
  error = comparison.h1_error("elliptic", "u");
  libMesh::out << "H1 error of u w.r.t reference: " << error << std::endl << std::endl;


#endif // #ifdef LIBMESH_ENABLE_AMR

  return 0;
}


