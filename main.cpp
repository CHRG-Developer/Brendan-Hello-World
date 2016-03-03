#include <iostream>
#include "Uniform_Mesh.h"
#include "vector_var.h"
#include <stdio.h>      /* printf */
#include <iostream>
#include <math.h>
#include "Boundary_Conditions.h"
#include "Solution.h"
#include "Solver.h"
#include "test_cases.h"

using namespace std;


int main()
{

    test_cases couette_flow;

    test_cases poiseuille_flow;
   poiseuille_flow.west_to_east_poiseuille_flow();



    //test_cases lid_driven_cav;
    //lid_driven_cav.lid_driven_cavity_N();
     // couette_flow.west_to_east_couette_flow();
 //   couette_flow.east_to_west_couette_flow();
    // couette_flow.north_to_south_couette_flow();
   // couette_flow.south_to_north_couette_flow();

    return 0;
}


