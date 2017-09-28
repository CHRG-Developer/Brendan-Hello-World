#include "program.h"
#include "quad_bcs.h"
#include "global_variables.h"
#include "preprocessor.h"
#include <ctime>
#include "quad_bcs_plus.h"
#include "domain_geometry.h"
#include "Mesh.h"
#include "Boundary_Conditions.h"
#include "external_forces.h"
#include "Solution.h"
#include "Solver.h"
#include <fstream>
#include <iostream>
#include "postprocessor.h"
#include <stdlib.h>
#include <stdio.h>
#include "TECIO.h"
#include <sstream>
#include <tecplot_output.h>

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
    std::clock_t start = clock();
    double duration;


    preprocessor pre_processor;
    //postprocessor post_processor;

    pre_processor.initialise_program_variables(xml_input, globals, domain,initial_conds,bcs);

    copyfile(xml_input,globals.output_file);
     // create Mesh
    Mesh mesh(domain,globals);


    // create boundary conditions
    Boundary_Conditions bc(mesh.get_num_x(), mesh.get_num_y());
    bc.assign_boundary_conditions(mesh.get_num_x(), mesh.get_num_y(),bcs);

    // assign external force terms
    external_forces source_term(mesh.get_total_cells());
    source_term.set_uniform_force(initial_conds.pressure_gradient); //pressure force applied through external force

    //create solution
    Solution soln(mesh.get_total_cells());
    Solution residual(mesh.get_total_cells());
    if ( globals.fmg_levels > 0 ){

        fmg_cycle(fmg,soln,soln,mesh,bcs,initial_conds,globals,domain, bc);



    }else{

        soln.assign_pressure_gradient(initial_conds.rho_gradient,initial_conds.origin_loc,
                                initial_conds.rho_origin_mag,mesh,globals);
        soln.assign_velocity_gradient(initial_conds.vel_gradient,initial_conds.origin_loc,
                                initial_conds.vel_origin_mag,mesh,globals);
        soln.set_average_rho(initial_conds.average_rho);


    }
    // Solvec

    Solver solve;
    tecplot_output grid(globals,mesh,soln,bc,1,0.0);

    solve.Uniform_Mesh_Solver_Clean_MK2(mesh,soln,bc,source_term,globals,domain,initial_conds,bcs,
                              mg,residual,fmg);

    soln.post_process(globals.pre_conditioned_gamma,mesh, globals,initial_conds);
    soln.output(globals.output_file, globals, domain);

    //post_processor.output_vtk_mesh(globals.output_file,globals,domain);

    std::clock_t end = clock();

    duration = double(end-start)/ CLOCKS_PER_SEC;

    output_globals(globals,duration);

    std::cout << duration << std::endl;
    std::cout << "CPU Cycles:" << double(end-start) << std::endl;


}
//void program::output_tecplot(global_variables &globals, Mesh &Mesh, Solution &Soln,
//                             Boundary_Conditions &bcs){
//
//    std::string output_location;
//    output_location = globals.output_file + "/output.plt";
//    std::string reynolds_text;
//    reynolds_text = globals.reynolds_number;
//    std::string zone_name;
//    std::stringstream ss;
//
//   enum fileType_e { FULL = 0, GRID = 1, SOLUTION = 2 };
//
//     float *nodx, *nody, *z, *p,*u, *v,*y ,*x, *u_err, *u_exact;
//    int *connectivity;
//    double solTime;
//    INTEGER4 debug, i, j, k, dIsDouble, vIsDouble, zoneType, strandID, parentZn, isBlock;
//    INTEGER4 iCellMax, jCellMax, kCellMax, nFConns, fNMode, shrConn, fileType;
//
//    int valueLocation[] = {1,1,1,0,0,0,0,0,0,0};
//    INTEGER4 fileFormat; // 0 == PLT, 1 == SZPLT
//    fileFormat = 0;
//
//    INTEGER4 nNodes, nCells, nFaces, connectivityCount, index;
//    int XDIM, YDIM, ZDIM; // nodes
//    int t;
//    XDIM = Mesh.get_num_x() +1 -2;
//    YDIM = Mesh.get_num_y() +1 -2;
//    ZDIM = 2;
//
//    debug     = 1;
//    vIsDouble = 0;
//    dIsDouble = 0;
//    nNodes = XDIM * YDIM * ZDIM;
//    nCells = (XDIM - 1) * (YDIM - 1) * (ZDIM - 1);
//    nFaces = 6; /* Not used */
//    zoneType  = 5;      /* Brick */
//    solTime   = 360.0;
//    strandID  = 0;     /* StaticZone */
//    parentZn  = 0;      /* No Parent */
//    isBlock   = 1;      /* Block */
//    iCellMax  = 0;
//    jCellMax  = 0;
//    kCellMax  = 0;
//    nFConns   = 0;
//    fNMode    = 0;
//    shrConn   = 0;
//    fileType  = FULL;
//    ss << nCells;
//    zone_name = ss.str();
//    /*
//     * Open the file and write the tecplot datafile
//     * header information
//     */
//    i = TECINI142((char*) "Couette Flow" ,
//                  (char*)"nodx nody z p u v x y u_error u_exact",
//                  (char*) output_location.c_str(),
//                  (char*) ".",
//                  &fileFormat,
//                  &fileType,
//                  &debug,
//                  &vIsDouble);
//
//    i = TECAUXSTR142("Re" ,  reynolds_text.c_str());
//
//    nodx  = (float*)malloc(nNodes * sizeof(float));
//    nody  = (float*)malloc(nNodes * sizeof(float));
//    z  = (float*)malloc(nNodes * sizeof(float));
//    p  = (float*)malloc(nNodes * sizeof(float));
//    u = (float*)malloc(nCells * sizeof(float));
//    v = (float*)malloc(nCells * sizeof(float));
//    y = (float*)malloc(nCells * sizeof(float));
//    x = (float*)malloc(nCells * sizeof(float));
//    u_err = (float*)malloc(nCells * sizeof(float));
//    u_exact = (float*)malloc(nCells * sizeof(float));
//    for (k = 0; k < ZDIM; k++)
//        for (j = 0; j < YDIM; j++)
//            for (i = 0; i < XDIM; i++)
//            {
//                index = (k * YDIM + j) * XDIM + i;
//                nodx[index] = (float)(i + 1);
//                nody[index] = (float)(j + 1);
//                z[index] = (float)(k + 1);
//
//            }
//
//    t = 0;
//    for (i = 0; i < Mesh.get_total_cells(); ++i){
//        if( bcs.get_bc(i) == false){
//            p[t] = (float) Soln.get_rho(i);
//            u[t] = (float) Soln.get_u(i)/globals.mach_number *sqrt(3);
//            v[t] = (float) Soln.get_v(i)/globals.mach_number *sqrt(3);
//            x[t] = (float) Mesh.get_centroid_x(i)/Mesh.get_X();
//            y[t] = (float) Mesh.get_centroid_y(i)/Mesh.get_Y();
//            u_err[t] = (float) Soln.get_u_error(i);
//            u_exact[t] = (float) Soln.get_u_exact(i)/globals.mach_number *sqrt(3);
//            t++;
//        }
//
//    }
//
//    connectivityCount = 8 * nCells;
//    connectivity = (INTEGER4*)malloc(connectivityCount * sizeof(INTEGER4));
//    for (k = 0; k < ZDIM - 1; k++)
//        for (j = 0; j < YDIM - 1; j++)
//            for (i = 0; i < XDIM - 1; i++)
//            {
//                index = ((k * (YDIM - 1) + j) * (XDIM - 1) + i) * 8;
//                connectivity[index] = (k * YDIM + j) * XDIM + i + 1;
//                connectivity[index + 1] = connectivity[index] + 1;
//                connectivity[index + 2] = connectivity[index] + XDIM + 1;
//                connectivity[index + 3] = connectivity[index] + XDIM;
//                connectivity[index + 4] = connectivity[index] + XDIM * YDIM;
//                connectivity[index + 5] = connectivity[index + 1] + XDIM * YDIM;
//                connectivity[index + 6] = connectivity[index + 2] + XDIM * YDIM;
//                connectivity[index + 7] = connectivity[index + 3] + XDIM * YDIM;
//            }
//
//
//       /*
//     * Write the zone header information.
//     */
//     std::cout << zone_name.c_str() << endl;
//   i = TECZNE142(
//                 //(char*) zone_name.c_str(),
//                 zone_name.c_str(),
//                  &zoneType,
//                  &nNodes,
//                  &nCells,
//                  &nFaces,
//                  &iCellMax,
//                  &jCellMax,
//                  &kCellMax,
//                  &solTime,
//                  &strandID,
//                  &parentZn,
//                  &isBlock,
//                  &nFConns,
//                  &fNMode,
//                  0,              /* TotalNumFaceNodes */
//                  0,              /* NumConnectedBoundaryFaces */
//                  0,              /* TotalNumBoundaryConnections */
//                  NULL,           /* PassiveVarList */
//                  valueLocation,  /* ValueLocation = Nodal */
//                  NULL,           /* SharVarFromZone */
//                  &shrConn);
///*
//     * Write out the field data.
//     */
//    i = TECDAT142(&nNodes, nodx, &dIsDouble);
//    i = TECDAT142(&nNodes, nody, &dIsDouble);
//    i = TECDAT142(&nNodes, z, &dIsDouble);
//    i = TECDAT142(&nCells, p, &dIsDouble);
//    i = TECDAT142(&nCells, u, &dIsDouble);
//    i = TECDAT142(&nCells, v, &dIsDouble);
//    i = TECDAT142(&nCells, x, &dIsDouble);
//    i = TECDAT142(&nCells, y, &dIsDouble);
//    i = TECDAT142(&nCells, u_err, &dIsDouble);
//    i = TECDAT142(&nCells, u_exact, &dIsDouble);
//    free(nodx);
//    free(nody);
//    free(z);
//    free(p);
//    free(u);
//    free(v);
//    free(x);
//    free(y);
//    free(u_err);
//    free(u_exact);
//
//    i = TECNODE142(&connectivityCount, connectivity);
//    free(connectivity);
//
//    i = TECEND142();
//}
void program::output_globals (global_variables globals,double duration){

    std::string output_location;
    std::string filename;
    std::ofstream globals_txt ;
    std::string globals_file;
    output_location = globals.output_file;
    filename = globals.simulation_name;
    globals_file = output_location + "/globals.txt";
    double cycles;

    cycles = duration * CLOCKS_PER_SEC;
/// Generic Load Case Input

    globals_txt.open(globals_file.c_str(), ios::out);

    globals_txt << "Runtime:" << duration << "s" << endl;
    globals_txt << "CPU Cycles:" << cycles << endl;
    globals_txt << "File:" << filename.c_str()  << endl;

    globals_txt << "Tau:"  << globals.tau << endl;
    globals_txt << "Mach" << globals.max_velocity *sqrt(3)<< endl;
    globals_txt << "Reynolds" << globals.reynolds_number << endl;

    std::string reynolds_text;
    reynolds_text = globals.reynolds_number;

    globals_txt.close();


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
                                      Mesh &fine_mesh, quad_bcs_plus &bcs,
                                      initial_conditions &initial_conds,
                                      global_variables globals, domain_geometry &fine_domain,
                                        Boundary_Conditions &fine_bcs){

    fmg = fmg +1;

    int mg = 0;


    // create new coarse Mesh with double up dimensions
    domain_geometry coarse_domain = fine_mesh.create_coarse_mesh_domain();

    Mesh coarse_mesh (coarse_domain,globals);

    globals.update_tau(coarse_domain);
    globals.magnify_time_step();
    Boundary_Conditions bc(coarse_mesh.get_num_x(), coarse_mesh.get_num_y());
    bc.assign_boundary_conditions(coarse_mesh.get_num_x(), coarse_mesh.get_num_y(),bcs);

    Solution coarse_soln( coarse_mesh.get_total_cells());
    Solution coarse_residual( coarse_mesh.get_total_cells());
    // goto coarsest level
    if( fmg < globals.fmg_levels){
            fmg_cycle(fmg,coarse_residual,coarse_soln,coarse_mesh,bcs,initial_conds,globals,
                      coarse_domain,bc);

    }else{

        //apply initial conditions at coarsest level
        coarse_soln.assign_pressure_gradient(initial_conds.rho_gradient,initial_conds.origin_loc,
                                initial_conds.rho_origin_mag,coarse_mesh,globals);
        coarse_soln.assign_velocity_gradient(initial_conds.vel_gradient,initial_conds.origin_loc,
                                initial_conds.vel_origin_mag,coarse_mesh,globals);
        coarse_soln.set_average_rho(initial_conds.average_rho);

    }

    external_forces source_term(coarse_mesh.get_total_cells());
    source_term.set_uniform_force(initial_conds.pressure_gradient); //pressure force applied through external force


    Solver solve_coarse;

    solve_coarse.Uniform_Mesh_Solver_Clean(coarse_mesh,coarse_soln,bc,source_term,globals,
                                     coarse_domain,initial_conds,bcs,mg,coarse_residual,fmg);

    Solution temp_soln(coarse_mesh.get_total_cells());
    // prolongation occurs here

    coarse_soln.update_bcs(bc,coarse_mesh,coarse_domain);

    soln.prolongation( coarse_soln, temp_soln, soln,coarse_mesh, fine_mesh,fine_bcs,true);
    soln.set_average_rho(initial_conds.average_rho);
    soln.update_bcs(fine_bcs,fine_mesh,fine_domain);
    globals.update_tau(fine_domain);
     globals.reduce_time_step();
    fmg = fmg -1;

}

