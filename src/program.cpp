#include "program.h"
#include "quad_bcs.h"
#include "global_variables.h"
#include "preprocessor.h"
#include <ctime>
#include "quad_bcs_plus.h"
#include "domain_geometry.h"
#include "Uniform_Mesh.h"
#include "Boundary_Conditions.h"
#include "external_forces.h"
#include "Solution.h"
#include "Solver.h"
#include <fstream>

program::program()
{
    //ctor
}

program::~program()
{
    //dtor
}


void program::run(char* xml_input){


    quad_bcs_plus bcs;
    global_variables globals;
    domain_geometry domain;
    initial_conditions initial_conds;
    int mg =0; // first multigrid cycle
    std::clock_t start;
    double duration;

    preprocessor pre_processor;

    pre_processor.initialise_program_variables(xml_input, globals, domain,initial_conds,bcs);

    copyfile(xml_input,globals.output_file);
     // create Mesh
    Uniform_Mesh mesh(domain);

    // create boundary conditions
    Boundary_Conditions bc(mesh.get_num_x(), mesh.get_num_y());
    bc.assign_boundary_conditions(mesh.get_num_x(), mesh.get_num_y(),bcs);

    // assign external force terms
    external_forces source_term(mesh.get_total_nodes());
    source_term.set_uniform_force(initial_conds.pressure_gradient); //pressure force applied through external force

    //create solution
    Solution soln(mesh.get_total_nodes());
    Solution residual(mesh.get_total_nodes());
    soln.assign_pressure_gradient(initial_conds.rho_gradient,initial_conds.origin_loc,
                            initial_conds.rho_origin_mag,mesh);
    soln.assign_velocity_gradient(initial_conds.vel_gradient,initial_conds.origin_loc,
                            initial_conds.vel_origin_mag,mesh);
    soln.set_average_rho(initial_conds.average_rho);
    // Solvec

    Solver solve;

    solve.Uniform_Mesh_Solver(mesh,soln,bc,source_term,globals,domain,initial_conds,bcs,
                              mg,residual);

    soln.post_process(globals.pre_conditioned_gamma);
    soln.output(globals.output_file);



}

    // copy in binary mode
void program::copyfile( char* SRC,  std::string  DEST)
{
    DEST.append("/input.xml");

    std::ifstream src(SRC, std::ios::binary);
    std::ofstream dest(DEST, std::ios::binary);
    dest << src.rdbuf();

}

