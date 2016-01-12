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

    test_cases w_to_e;
    w_to_e.west_to_east_1d();
    test_cases lid_driven_cav;
    //lid_driven_cav.lid_driven_cavity_N();


    return 0;
}


