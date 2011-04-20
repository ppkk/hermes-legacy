#define HERMES_REPORT_ALL
#include "hermes2d.h"

using namespace RefinementSelectors;

const double THRESHOLD = 0.2;                     // This is a quantitative parameter of the adapt(...) function and
                                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                           // Adaptive strategy:
                                                  // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                                  //   error is processed. If more elements have similar errors, refine
                                                  //   all to keep the mesh symmetric.
                                                  // STRATEGY = 1 ... refine all elements whose error is larger
                                                  //   than THRESHOLD times maximum element error.
                                                  // STRATEGY = 2 ... refine all elements whose error is larger
                                                  //   than THRESHOLD.
                                                  // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO_H;        // Predefined list of element refinement candidates. Possible values are
                                                  // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO, H2D_HP_ANISO_H
                                                  // H2D_HP_ANISO_P, H2D_HP_ANISO. See User Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double ERR_STOP = 1.0;                      // Stopping criterion for adaptivity (rel. error tolerance between the
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
                                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,

const bool HERMES_VISUALIZATION = true;           // Set to "false" to suppress Hermes OpenGL visualization. 
const bool VTK_OUTPUT = true;                     // Set to "true" to enable VTK output.
const int P_INIT = 3;                             // Uniform polynomial degree of mesh elements.

std::complex<double> jj(0,1);
const double FREQ = 1.;

// Boundary markers.
const int BDY_BOTTOM = 1, BDY_OUTER = 2, BDY_LEFT = 3, BDY_INNER = 4;

template<typename Real>
bool civka(Real x, Real y)
{
  const Real C_TOP = 0.125;
  const Real C_BOTTOM = -0.125;
  const Real C_LEFT = 0.125;
  const Real C_RIGHT = 0.25;
  if ((x > C_LEFT) && (x < C_RIGHT) && (y < C_TOP) && (y > C_BOTTOM))
    return true;
  return false;
}

template<typename Real>
bool naboj(Real x, Real y)
{
  const Real C_TOP = -0.25;
  const Real C_BOTTOM = -0.375;
  const Real C_LEFT = 0;
  const Real C_RIGHT = -100; //0.0625;
  if ((x > C_LEFT) && (x < C_RIGHT) && (y < C_TOP) && (y > C_BOTTOM))
    return true;
  return false;
}

//proudova hustota
template<typename Real, typename Scalar>
Scalar F(Real x, Real y)
{
  if (civka(x, y)) 
    return 1.2e7;
  else
    return 0.;
}
    
//vodivost    
template<typename Real, typename Scalar>
Scalar Gamma(Real x, Real y)
{
  if (civka(x, y))
    return 0.; //5.6e7;
  else if (naboj(x,y))
    return 1.5e6;
  else
    return 0.;
}

template<typename Real, typename Scalar>
Scalar RelPermInv(Real x, Real y)
{
  if (naboj(x, y))
    return 0.001;
  else
    return 1.;
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Perform initial mesh refinements (optional).
//    mesh.refine_all_elements();
//    mesh.refine_all_elements();
//    mesh.refine_all_elements();
  mesh.refine_towards_boundary(4,4);

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(Hermes::vector<int>(BDY_BOTTOM, BDY_OUTER, BDY_LEFT, BDY_INNER));

  // Enter zero Dirichlet boundary values.
  BCValues bc_values;

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bc_types, &bc_values, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof = %d", ndof);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(bilinear_form));
  wf.add_vector_form(callback(linear_form));

  
    // Initialize coarse and reference mesh solution.
  Solution sln, ref_sln;
  
  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize views.
  ScalarView view_r("Solution real", new WinGeom(0, 0, 600, 1000));
  ScalarView view_i("Solution imag", new WinGeom(600, 0, 600, 1000));
  view_r.fix_scale_width(50);
  view_r.show_mesh(false);
  view_i.fix_scale_width(50);
  view_i.show_mesh(false);
  OrderView  oview("Polynomial orders", new WinGeom(1200, 0, 600, 1000));
 
    
  // DOF and CPU convergence graphs initialization.
  SimpleGraph graph_dof, graph_cpu;
  
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Adaptivity loop:
  int as = 1; 
  bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh and setup reference space.
    Space* ref_space = construct_refined_space(&space);

    // Initialize matrix solver.
    SparseMatrix* matrix = create_matrix(matrix_solver);
    Vector* rhs = create_vector(matrix_solver);
    Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

    // Assemble reference problem.
    info("Solving on reference mesh.");
    bool is_linear = true;
    DiscreteProblem* dp = new DiscreteProblem(&wf, ref_space, is_linear);
    dp->assemble(matrix, rhs);

    // Time measurement.
    cpu_time.tick();
    
    // Solve the linear system of the reference problem. 
    // If successful, obtain the solution.
    if(solver->solve()) Solution::vector_to_solution(solver->get_solution(), ref_space, &ref_sln);
    else error ("Matrix solver failed.\n");

    // Project the fine mesh solution onto the coarse mesh.
    info("Projecting reference solution on coarse mesh.");
    OGProjection::project_global(&space, &ref_sln, &sln, matrix_solver); 

    // Time measurement.
    cpu_time.tick();
   
    // VTK output.
    if (VTK_OUTPUT) {
      // Output solution in VTK format.
      Linearizer lin;
      char* title = new char[100];
      sprintf(title, "sln-%d.vtk", as);
      lin.save_solution_vtk(&sln, title, "Potential", false);
      info("Solution in VTK format saved to file %s.", title);

      // Output mesh and element orders in VTK format.
      Orderizer ord;
      sprintf(title, "ord-%d.vtk", as);
      ord.save_orders_vtk(&space, title);
      info("Element orders in VTK format saved to file %s.", title);
    }

    // View the coarse mesh solution and polynomial orders.
    if (HERMES_VISUALIZATION) {
      RealFilter real(&sln);
      view_r.show(&real);
      ImagFilter imag(&sln);
      view_i.show(&imag);
      oview.show(&space);
    }

    // Skip visualization time.
    cpu_time.tick(HERMES_SKIP);

    // Calculate element errors and total error estimate.
    info("Calculating error estimate."); 
    Adapt* adaptivity = new Adapt(&space);
    bool solutions_for_adapt = true;
    // In the following function, the Boolean parameter "solutions_for_adapt" determines whether 
    // the calculated errors are intended for use with adaptivity (this may not be the case, for example,
    // when error wrt. an exact solution is calculated). The default value is solutions_for_adapt = true, 
    // The last parameter "error_flags" determine whether the total and element errors are treated as 
    // absolute or relative. Its default value is error_flags = HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL.
    // In subsequent examples and benchmarks, these two parameters will be often used with 
    // their default values, and thus they will not be present in the code explicitly.
    double err_est_rel = adaptivity->calc_err_est(&sln, &ref_sln, solutions_for_adapt, 
                         HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL) * 100;

    // Report results.
    info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%", 
      Space::get_num_dofs(&space), Space::get_num_dofs(ref_space), err_est_rel);

    // Time measurement.
    cpu_time.tick();

    // Add entry to DOF and CPU convergence graphs.
    graph_dof.add_values(Space::get_num_dofs(&space), err_est_rel);
    graph_dof.save("conv_dof_est.dat");
    graph_cpu.add_values(cpu_time.accumulated(), err_est_rel);
    graph_cpu.save("conv_cpu_est.dat");

    // If err_est too large, adapt the mesh.
    if (err_est_rel < ERR_STOP) done = true;
    else 
    {
      info("Adapting coarse mesh.");
      done = adaptivity->adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
      
      // Increase the counter of performed adaptivity steps.
      if (done == false)  as++;
    }
    if (Space::get_num_dofs(&space) >= NDOF_STOP) done = true;

    // Clean up.
    delete solver;
    delete matrix;
    delete rhs;
    delete adaptivity;
    if(done == false) delete ref_space->get_mesh();
    delete ref_space;
    delete dp;
    
  }
  while (done == false);
  
  verbose("Total running time: %g s", cpu_time.accumulated());

  // Show the reference solution - the final result.
  view_r.set_title("Fine mesh solution - real");
  RealFilter real(&ref_sln);
//   SimpleFilter bxReal(Hermes::vector<MeshFunction*>(&real), Hermes::vector<int>(H2D_FN_DY));
  MagFilter bxReal(Hermes::vector<MeshFunction*>(&ref_sln, &ref_sln), Hermes::vector<int>(H2D_FN_DY, H2D_FN_DY));
  view_r.show_mesh(false);
  view_r.show(&bxReal);
  
  
  view_i.set_title("Fine mesh solution - imag");
  ImagFilter imag(&ref_sln);
  MagFilter by(Hermes::vector<MeshFunction*>(&ref_sln, &ref_sln), Hermes::vector<int>(H2D_FN_DX, H2D_FN_DX));
  view_i.show_mesh(false);
  view_i.show(&by);

  // Wait for all views to be closed.
  View::wait();

  return 0;
}

  
  
/*  
  // Initialize the FE problem.
  bool is_linear = true;
  DiscreteProblem dp(&wf, &space, is_linear);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Initialize the solution.
  Solution sln;

  // Assemble the stiffness matrix and right-hand side vector.
  info("Assembling the stiffness matrix and right-hand side vector.");
  dp.assemble(matrix, rhs);

  // Solve the linear system and if successful, obtain the solution.
  info("Solving the matrix problem.");
  if(solver->solve()) Solution::vector_to_solution(solver->get_solution(), &space, &sln);
  else error ("Matrix solver failed.\n");

  // VTK output.
  if (VTK_OUTPUT) {
    // Output solution in VTK format.
    Linearizer lin;
    bool mode_3D = true;
    lin.save_solution_vtk(&sln, "sln.vtk", "Temperature", mode_3D);
    info("Solution in VTK format saved to file %s.", "sln.vtk");

    // Output mesh and element orders in VTK format.
    Orderizer ord;
    ord.save_orders_vtk(&space, "ord.vtk");
    info("Element orders in VTK format saved to file %s.", "ord.vtk");
  }

  // Visualize the solution.
  if (HERMES_VISUALIZATION) {
    ScalarView view_r("Solution real", new WinGeom(0, 0, 800, 1000));
    RealFilter real(&sln);
    view_r.show(&real);

    ScalarView view_i("Solution imag", new WinGeom(800, 0, 800, 1000));
    ImagFilter imag(&sln);
    view_i.show(&imag);

    View::wait();
  }

  // Clean up.
  delete solver;
  delete matrix;
  delete rhs;

  return 0;
}
*/
