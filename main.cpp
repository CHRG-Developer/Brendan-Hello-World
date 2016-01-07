#include <iostream>
#include "Uniform_Mesh.h"
#include "vector_var.h"
#include <stdio.h>      /* printf */
#include <iostream>
#include <math.h>
#include "Boundary_Conditions.h"
#include "Solution.h"
#include "Solver.h"


using namespace std;

void vector_var_tests();

int main()
{
    double X,Y,dx,dy,dt;
    double reynolds,kine_viscosity;
    double U;


    reynolds = 1000;
    X= 4;
    Y=4;
    dx=1;
    dy = 1;
    dt = 1;
    U = 10;

    kine_viscosity = U * X/ reynolds;

    //vector_var_tests();

    // create Mesh
    Uniform_Mesh mesh(X,Y,dx,dy);

    // create boundary conditions
    Boundary_Conditions bc(mesh.get_num_x(), mesh.get_num_y());
    bc.assign_boundary_conditions(mesh.get_num_x(), mesh.get_num_y());


    //create solution
    Solution soln(mesh.get_total_nodes());

    // Solve

    Solver solve;
    solve.Uniform_Mesh_Solver(dt,kine_viscosity,mesh,soln,bc);









    return 0;
}


void vector_var_tests(){
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
