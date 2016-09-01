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
    int fmg =0; // first full multigrid cycle
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
    if ( globals.fmg_levels > 0 ){

        fmg_cycle(fmg,soln,soln,mesh,bcs,initial_conds,globals,domain);



    }else{

        soln.assign_pressure_gradient(initial_conds.rho_gradient,initial_conds.origin_loc,
                                initial_conds.rho_origin_mag,mesh);
        soln.assign_velocity_gradient(initial_conds.vel_gradient,initial_conds.origin_loc,
                                initial_conds.vel_origin_mag,mesh);
        soln.set_average_rho(initial_conds.average_rho);


    }
    // Solvec

    Solver solve;

    solve.Uniform_Mesh_Solver(mesh,soln,bc,source_term,globals,domain,initial_conds,bcs,
                              mg,residual,fmg);

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

void program::fmg_cycle(int &fmg,Solution &residual , Solution &soln,
                                      Uniform_Mesh &fine_mesh, quad_bcs_plus &bcs,
                                      initial_conditions &initial_conds,
                                      global_variables globals, domain_geometry &fine_domain ){

    fmg = fmg +1;

    int mg = 0;


    // create new coarse Mesh with double up dimensions
    domain_geometry coarse_domain = fine_mesh.create_coarse_mesh_domain();

    Uniform_Mesh coarse_mesh (coarse_domain);

    globals.update_tau(coarse_domain);

    Boundary_Conditions bc(coarse_mesh.get_num_x(), coarse_mesh.get_num_y());
    bc.assign_boundary_conditions(coarse_mesh.get_num_x(), coarse_mesh.get_num_y(),bcs);

    Solution coarse_soln( coarse_mesh.get_total_nodes());
    Solution coarse_residual( coarse_mesh.get_total_nodes());
    // goto coarsest level
    if( fmg < globals.fmg_levels){
            fmg_cycle(fmg,coarse_residual,coarse_soln,coarse_mesh,bcs,initial_conds,globals,
                      coarse_domain);

    }else{

        //apply initial conditions at coarsest level
        coarse_soln.assign_pressure_gradient(initial_conds.rho_gradient,initial_conds.origin_loc,
                                initial_conds.rho_origin_mag,coarse_mesh);
        coarse_soln.assign_velocity_gradient(initial_conds.vel_gradient,initial_conds.origin_loc,
                                initial_conds.vel_origin_mag,coarse_mesh);
        coarse_soln.set_average_rho(initial_conds.average_rho);

    }

    external_forces source_term(coarse_mesh.get_total_nodes());
    source_term.set_uniform_force(initial_conds.pressure_gradient); //pressure force applied through external force


    Solver solve_coarse;

    solve_coarse.Uniform_Mesh_Solver(coarse_mesh,coarse_soln,bc,source_term,globals,
                                     coarse_domain,initial_conds,bcs,mg,coarse_residual,fmg);

    Solution temp_soln(coarse_mesh.get_total_nodes());
    // prolongation occurs here
    soln.prolongation( coarse_soln, temp_soln, soln,coarse_mesh, fine_mesh,bc,true);
    soln.set_average_rho(initial_conds.average_rho);
    globals.update_tau(fine_domain);
    fmg = fmg -1;

}

