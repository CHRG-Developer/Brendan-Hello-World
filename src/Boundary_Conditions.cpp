#include "Boundary_Conditions.h"
#include <stdlib.h>

Boundary_Conditions::Boundary_Conditions(int num_x_nodes, int num_y_nodes)
{
    //ctor

    int total_nodes = num_x_nodes * num_y_nodes;



      n_bc = (bool*) malloc (total_nodes+1);
        if (n_bc==NULL) exit (1);
      s_bc = (bool*) malloc (total_nodes+1);
        if (s_bc==NULL) exit (1);
      w_bc = (bool*) malloc (total_nodes+1);
        if (w_bc==NULL) exit (1);
      e_bc = (bool*) malloc (total_nodes+1);
        if (e_bc==NULL) exit (1);

     n_rho = (double*) malloc (total_nodes+1);
        if (n_rho==NULL) exit (1);
      s_rho = (double*) malloc (total_nodes+1);
        if (s_rho==NULL) exit (1);
      w_rho = (double*) malloc (total_nodes+1);
        if (w_rho==NULL) exit (1);
      e_rho = (double*) malloc (total_nodes+1);
        if (e_rho==NULL) exit (1);

      n_u = (double*) malloc (total_nodes+1);
        if (n_u==NULL) exit (1);
      s_u = (double*) malloc (total_nodes+1);
        if (s_u==NULL) exit (1);
      w_u = (double*) malloc (total_nodes+1);
        if (w_u==NULL) exit (1);
      e_u = (double*) malloc (total_nodes+1);
        if (e_u==NULL) exit (1);

      n_v = (double*) malloc (total_nodes+1);
        if (n_v==NULL) exit (1);
      s_v = (double*) malloc (total_nodes+1);
        if (s_v==NULL) exit (1);
      w_v = (double*) malloc (total_nodes+1);
        if (w_v==NULL) exit (1);
      e_v = (double*) malloc (total_nodes+1);
        if (e_v==NULL) exit (1);





}

Boundary_Conditions::~Boundary_Conditions()
{
    //dtor
}


void Boundary_Conditions::assign_boundary_conditions(){


}
