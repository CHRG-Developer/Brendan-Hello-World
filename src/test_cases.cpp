#include "test_cases.h"
#include <iostream>
#include "Uniform_Mesh.h"
#include "vector_var.h"
#include <stdio.h>      /* printf */
#include <iostream>
#include <math.h>
#include "Boundary_Conditions.h"
#include "Solution.h"
#include "Solver.h"
#include "quad_bcs.h"

using namespace std;

test_cases::test_cases()
{
    //ctor
}

test_cases::~test_cases()
{
    //dtor
}
void test_cases::lid_driven_cavity_N(){

    double X,Y,dx,dy,dt; // dt is streaming time step
    double reynolds,kine_viscosity,tau;
    double U;
    double simulation_length;
    double delta_t; // time stepping step
    quad_bcs bcs;
    double cs;

    /// Parameters unique to test case
    reynolds = 1000;
    X= 100;
    Y=100;
    dx=1; // grid spacing
    dy = 1;  // grid spacing
    dt = 0.05;  // streaming time step let dt =dx = dy i.e. lattice spacing
    U = 1;
    simulation_length = 200;
    kine_viscosity = U * X/ reynolds;
    //kine_viscosity = 20;
    delta_t = 0.1;
    cs = 1/sqrt(3);
    tau = kine_viscosity + 0.5* pow(cs,2) *dt;
    //tau =0.75;

    // set boundary conditions for this test case
    bcs.w_rho = 1;
    bcs.w_u = 0;
    bcs.w_v = 0;
    bcs.w_w = 0;

    bcs.s_rho = 1;
    bcs.s_u = 0;
    bcs.s_v = 0;
    bcs.s_w = 0;

    bcs.e_rho = 1;
    bcs.e_u = 0;
    bcs.e_v = 0;
    bcs.e_w = 0;

    bcs.n_rho = 1;
    bcs.n_u = 1;
    bcs.n_v = 0;
    bcs.n_w = 0;

    /// Methods to run the test case
    //vector_var_tests();

    // create Mesh
    Uniform_Mesh mesh(X,Y,dx,dy);

    // create boundary conditions
    Boundary_Conditions bc(mesh.get_num_x(), mesh.get_num_y());
    bc.assign_boundary_conditions(mesh.get_num_x(), mesh.get_num_y(),bcs);


    //create solution
    Solution soln(mesh.get_total_nodes());

    // Solve

    Solver solve;
    solve.Uniform_Mesh_Solver(dt,tau,mesh,soln,bc,simulation_length, delta_t);

    tau = 1;

}
void test_cases::west_to_east_1d(){

    double X,Y,dx,dy,dt; // dt is streaming time step
    double reynolds,kine_viscosity,tau;
    double U;
    double simulation_length;
    double delta_t; // time stepping step
    quad_bcs bcs;
    double cs;

    /// Parameters unique to test case
    reynolds = 1000;
    X= 100;
    Y=1;
    dx=1; // grid spacing
    dy = 1;  // grid spacing
    dt = 0.05;  // streaming time step let dt =dx = dy i.e. lattice spacing
    U = 1;
    simulation_length = 200;
    kine_viscosity = U * X/ reynolds;
    //kine_viscosity = 20;
    delta_t = 0.1;
    cs = 1/sqrt(3);
    tau = kine_viscosity + 0.5* pow(cs,2) *dt;
    //tau =0.75;

    // set boundary conditions for this test case
    bcs.w_rho = 0;
    bcs.w_u = 1;
    bcs.w_v = 0;
    bcs.w_w = 0;

    bcs.s_rho = 0;
    bcs.s_u = 0;
    bcs.s_v = 0;
    bcs.s_w = 0;

    bcs.e_rho = 0;
    bcs.e_u = 0;
    bcs.e_v = 0;
    bcs.e_w = 0;

    bcs.n_rho = 0;
    bcs.n_u = 0;
    bcs.n_v = 0;
    bcs.n_w = 0;

    /// Methods to run the test case
    //vector_var_tests();

    // create Mesh
    Uniform_Mesh mesh(X,Y,dx,dy);

    // create boundary conditions
    Boundary_Conditions bc(mesh.get_num_x(), mesh.get_num_y());
    bc.assign_boundary_conditions(mesh.get_num_x(), mesh.get_num_y(),bcs);


    //create solution
    Solution soln(mesh.get_total_nodes());

    // Solve

    Solver solve;
    solve.Uniform_Mesh_Solver(dt,tau,mesh,soln,bc,simulation_length, delta_t);

    tau = 1;

}

void test_cases::vector_var_tests(){
     vector_var a,b,c,d,e,f;

    a.x = 3;
    a.y = 4;
    a.z = 0;

    b.x = 2;
    b.y = 2;
    b.z = 2;

    c.x = 4;
    c.y = -3;
    c.z = 0;

    d.x= 0;
    d.y = 0;
    d.z = 0;

    f.x = 5;
    f.y = 5;
    f.z = 0;




    double mag , dp, ang;
    mag = a.Magnitude();

    cout <<  "Magnitude of a:" << mag << endl ;
    dp = a.Dot_Product(b);
    cout << "Dot product of a and b: " << dp << endl ;
    ang = a.Angle_Between_Vectors(c);
    ang = ang *360/2/M_PI;
    cout << "angle between a and c: " << ang << endl ;
    e.Get_Gradient(10,20,d,b);
    cout << "gradient vector between d and b =" << e.x << "," << e.y << "," << e.z << endl ;



}
