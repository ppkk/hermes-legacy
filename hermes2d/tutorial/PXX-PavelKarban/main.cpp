#define HERMES_REPORT_ALL
#include "hermes2d.h"

// This example shows how to solve a simple PDE that describes stationary 
// heat transfer in an object consisting of two materials (aluminum and 
// copper). The object is heated by constant volumetric heat sources
// generated by a DC electric current. The temperature on the boundary 
// is fixed. We will learn how to:
//
//   - load the mesh,
//   - perform initial refinements,
//   - create a H1 space over the mesh,
//   - define weak formulation,
//   - initialize matrix solver,
//   - assemble and solve the matrix system,
//   - output the solution and element orders in VTK format 
//     (to be visualized, e.g., using Paraview),
//   - visualize the solution using Hermes' native OpenGL-based functionality.
//
// PDE: Poisson equation -div(LAMBDA grad u) - VOLUME_HEAT_SRC = 0.
//
// Boundary conditions: Dirichlet u(x, y) = FIXED_BDY_TEMP on the boundary.
//
// The following parameters can be changed:

const bool HERMES_VISUALIZATION = true;           // Set to "false" to suppress Hermes OpenGL visualization. 
const bool VTK_VISUALIZATION = false;              // Set to "true" to enable VTK output.
const int P_MAG_INIT = 3;                             // Uniform polynomial degree of mesh elements.
const int P_TEMP_INIT = 3;
const int INIT_REF_NUM = 1;                       // Number of initial uniform mesh refinements.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem parameters.
const double A_INIT = 0.0;
const double TEMP_INIT = 20.0;
const double DK_INIT = 0.0;

const double TIME_STEP = 0.1;
const double TIME_FINAL = 20.;

#include "tables.cpp"

// Weak forms.
#include "definitions.cpp"

int main(int argc, char* argv[])
{
    // Instantiate a class with global functions.
    Hermes2D hermes2d;

    // Initialize tables from the file tables.cpp
    initTables();

    // Load the mesh.
    Mesh mesh_mag, mesh_temp, mesh_elast_drzak, mesh_elast_stopka;
    H2DReader mloader;
    mloader.load("mesh_mag.mesh", &mesh_mag);
    mloader.load("mesh_temp.mesh", &mesh_temp);
    mloader.load("mesh_elast_drzak.mesh", &mesh_elast_drzak);
    mloader.load("mesh_elast_stopka.mesh", &mesh_elast_stopka);

//    MeshView mv;
//    mv.show(&mesh_mag);
//    MeshView mv_temp;
//    mv_temp.show(&mesh_temp);
//    MeshView mv_elast_stopka;
//    mv_elast_stopka.show(&mesh_elast_stopka);
//    MeshView mv_elast_drzak;
//    mv_elast_drzak.show(&mesh_elast_drzak);

    // Perform initial mesh refinements (optional).
     for (int i=0; i < INIT_REF_NUM; i++){
         mesh_temp.refine_all_elements();
         mesh_mag.refine_all_elements();
         mesh_elast_drzak.refine_all_elements();
         mesh_elast_stopka.refine_all_elements();
    }
    // Initialize the weak formulation.
    WeakFormMagnetic wf(2);
    wf.registerForms();

    // Initialize boundary conditions.
    DefaultEssentialBCConst bc_essential(Hermes::vector<std::string>("16", "17", "18", "39", "40", "41"), 0.0);
    EssentialBCs bcs(&bc_essential);

    // Create an H1 space with default shapeset.
    H1Space space_mag_real(&mesh_mag, &bcs, P_MAG_INIT);
    H1Space space_mag_imag(&mesh_mag, &bcs, P_MAG_INIT);
    // ndof
    int ndof = Space::get_num_dofs(Hermes::vector<Space *>(&space_mag_real, &space_mag_imag));
    std::cout << "ndof: " << ndof << std::endl;

    Solution *sln_mag_real = new Solution();
    sln_mag_real->set_const(&mesh_mag, A_INIT);
    Solution *sln_mag_imag = new Solution();
    sln_mag_imag->set_const(&mesh_mag, A_INIT);

    // Set up the solver, matrix, and rhs according to the solver selection.
    SparseMatrix* matrix = create_matrix(matrix_solver);
    Vector* rhs = create_vector(matrix_solver);
    Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

    // Initialize the FE problem.
    DiscreteProblem dp(&wf, Hermes::vector<Space *>(&space_mag_real, &space_mag_imag));
    dp.assemble(matrix, rhs);

    if (solver->solve())
        Solution::vector_to_solutions(solver->get_solution(),
                                      Hermes::vector<Space *>(&space_mag_real, &space_mag_imag),
                                      Hermes::vector<Solution *>(sln_mag_real, sln_mag_imag));
    else
        error ("Matrix solver failed.\n");

    WjFilter wjfilter(sln_mag_real, sln_mag_imag);
    MagneticVectorPotentialFilter afilter(sln_mag_real);

    // Visualize the solution.
    ScalarView view_a("Ar - real", new WinGeom(0, 0, 440, 750));
    view_a.show(&afilter, HERMES_EPS_NORMAL);
    ScalarView view_wj("wj", new WinGeom(450, 0, 440, 750));
    //view_wj.set_min_max_range(-0.000001, 0.000001);
    view_wj.show(&wjfilter, HERMES_EPS_NORMAL);
 //   View::wait();


//****************** TEMPERATURE **********************************************

    //Solution sln_temp(&mesh_temp, TEMP_INIT);
    Solution* sln_temp = new Solution();
    sln_temp->set_const(&mesh_temp, TEMP_INIT);
    WeakFormTemp wf_temp(TIME_STEP);
    wf_temp.registerForms(sln_temp, &wjfilter);

    double current_time = 0;

    EssentialBCs bcs_temp;

    for(int i = 0; i < NUM_EDGES; i++){
        if(heatEdge[i].type == PhysicFieldBC_Heat_Temperature){
            char label[5];
            sprintf(label, "%d", i);
            DefaultEssentialBCConst *bc = new DefaultEssentialBCConst(label, heatEdge[i].temperature);
            bcs_temp.add_boundary_condition(bc);
        }
    }

//    // Initialize temperature boundary conditions.
//    DefaultEssentialBCConst bc_essential_temp1(Hermes::vector<std::string>("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"), TEMP_INIT);
//    DefaultEssentialBCConst bc_essential_temp2(Hermes::vector<std::string>("11", "12", "13", "14", "15", "35", "36", "37", "38", "40"), TEMP_INIT);
//    EssentialBCs bcs_temp(Hermes::vector<EssentialBoundaryCondition*>(&bc_essential_temp1, &bc_essential_temp2));

    // Create an H1 space with default shapeset.
    H1Space space_temp(&mesh_temp, &bcs_temp, P_TEMP_INIT);
    int ndof_temp = space_temp.get_num_dofs();
    info("temperature ndof = %d", ndof_temp);

    // Initialize the FE problem.
    bool is_linear = true;
    DiscreteProblem dp_temp(&wf_temp, &space_temp, is_linear);

    // Set up the solver, matrix, and rhs according to the solver selection.
    SparseMatrix* matrix_temp = create_matrix(matrix_solver);
    Vector* rhs_temp = create_vector(matrix_solver);
    Solver* solver_temp = create_linear_solver(matrix_solver, matrix_temp, rhs_temp);
    solver_temp->set_factorization_scheme(HERMES_REUSE_FACTORIZATION_COMPLETELY);

    // Initialize views.
    ScalarView Tview("Temperature", new WinGeom(0, 0, 450, 600));
    //Tview.set_min_max_range(0,30);
    //Tview.fix_scale_width(30);
    Tview.show(sln_temp);
    //Tview.wait();

    // Time stepping:
    int ts = 1;
    do
    {
      info("---- Time step %d, time %3.5f s", ts, current_time);

      // First time assemble both the stiffness matrix and right-hand side vector,
      // then just the right-hand side vector.
      //wf.set_current_time(current_time);
      info("Assembling the stiffness matrix and right-hand side vector.");
      dp_temp.assemble(matrix_temp, rhs_temp);
      FILE *matfile;
      matfile = fopen("matice.txt", "w");
      matrix_temp->dump(matfile, "matrix");
      rhs_temp->dump(matfile, "rhs");

      // Solve the linear system and if successful, obtain the solution.
      info("Solving the temperature matrix problem.");
      if(solver_temp->solve())
          Solution::vector_to_solution(solver_temp->get_solution(), &space_temp, sln_temp);
      else error ("Matrix solver failed.\n");

      // Visualize the solution.
      char title[100];
      sprintf(title, "Time %3.2f s", current_time);
      Tview.set_title(title);
      Tview.show(sln_temp);
    //  Tview.wait();

      // Increase current time and time step counter.
      current_time += TIME_STEP;
      ts++;
    }
    while (current_time < TIME_FINAL);

    Tview.wait();

    // Clean up.
    delete solver;
    delete matrix;
    delete rhs;

    return 0;
}
