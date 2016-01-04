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


void Boundary_Conditions::assign_boundary_conditions(int num_x, int num_y){


    /// this method should be user defined to reflect the geometry of the problem
    int t = 0;

    //lid driven cavity conditions
    for (int i =0; i < num_x; i++){
        for( int j=0; j < num_y; j++){

            // West boundary
            if( i ==0){
                w_bc[t] = true;
                w_rho[t] = 1;
                w_u[t] = 0;
                w_v[t] = 0;


            }

            // east boundary
            if ( i == num_x){
                e_bc[t] = true;
                e_rho[t] = 1;
                e_u[t] = 0;
                e_v[t] = 0;


            }

            // south boundary
            if(j == 0){
                s_bc[t] = true;
                s_rho[t] = 1;
                s_u[t] = 0;
                s_v[t] = 0;

            }

            // north boundary
            if( j == num_y){
                n_bc[t] = true;
                n_rho[t] = 1;
                n_u[t] = 10; //key value of U in lid driven cavity
                n_v[t] = 0;

            }

            t++;
        }


    }

}
