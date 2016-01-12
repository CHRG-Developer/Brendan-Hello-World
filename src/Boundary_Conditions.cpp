#include "Boundary_Conditions.h"
#include <stdlib.h>


Boundary_Conditions::Boundary_Conditions(int num_x_nodes, int num_y_nodes)
{
    //ctor

    int total_nodes = num_x_nodes * num_y_nodes;



      n_bc = (bool*) malloc (sizeof(bool)*(total_nodes));
        if (n_bc==NULL) exit (1);
      s_bc = (bool*) malloc (sizeof(bool)*(total_nodes));
        if (s_bc==NULL) exit (1);
      w_bc = (bool*) malloc (sizeof(bool)*(total_nodes));
        if (w_bc==NULL) exit (1);
      e_bc = (bool*) malloc (sizeof(bool)*(total_nodes));
        if (e_bc==NULL) exit (1);

     n_rho = (double*) malloc (sizeof(double)*(total_nodes));
        if (n_rho==NULL) exit (1);
      s_rho = (double*) malloc (sizeof(double)*(total_nodes));
        if (s_rho==NULL) exit (1);
      w_rho = (double*) malloc (sizeof(double)*(total_nodes));
        if (w_rho==NULL) exit (1);
      e_rho = (double*) malloc (sizeof(double)*(total_nodes));
        if (e_rho==NULL) exit (1);

      n_u = (double*) malloc (sizeof(double)*(total_nodes));
        if (n_u==NULL) exit (1);
      s_u = (double*) malloc (sizeof(double)*(total_nodes));
        if (s_u==NULL) exit (1);
      w_u = (double*) malloc (sizeof(double)*(total_nodes));
        if (w_u==NULL) exit (1);
      e_u = (double*) malloc (sizeof(double)*(total_nodes));
        if (e_u==NULL) exit (1);

      n_v = (double*) malloc (sizeof(double)*(total_nodes));
        if (n_v==NULL) exit (1);
      s_v = (double*) malloc (sizeof(double)*(total_nodes));
        if (s_v==NULL) exit (1);
      w_v = (double*) malloc (sizeof(double)*(total_nodes));
        if (w_v==NULL) exit (1);
      e_v = (double*) malloc (sizeof(double)*(total_nodes));
        if (e_v==NULL) exit (1);





}

Boundary_Conditions::~Boundary_Conditions()
{
    //dtor
}


void Boundary_Conditions::assign_boundary_conditions(int num_x, int num_y, quad_bcs _bc){


    /// this method should be user defined to reflect the geometry of the problem
    /// currently set up for quad problem domain-> potential overload of this operator.
    int t = 0;

    //lid driven cavity conditions
    for (int i =0; i < num_x; i++){
        for( int j=0; j < num_y; j++){

            // West boundary
            if( i ==0){
                w_bc[t] = true;
                w_rho[t] = _bc.w_rho;
                w_u[t] = _bc.w_u;
                w_v[t] = _bc.w_v;


            }else{
                w_bc[t] =false;
            }

            // east boundary
            if ( i == (num_x -1)){
                e_bc[t] = true;
                e_rho[t] = _bc.e_rho;
                e_u[t] = _bc.e_u;
                e_v[t] = _bc.e_v;


            }else {
                    e_bc[t] = false;
            }

            // south boundary
            if(j == 0){
                s_bc[t] = true;
                s_rho[t] = _bc.s_rho;
                s_u[t] = _bc.s_u;
                s_v[t] = _bc.s_v;

            }else{
                s_bc[t] = false;
            }

            // north boundary
            if( j == (num_y-1)){
                n_bc[t] = true;
                n_rho[t] = _bc.n_rho;
                n_u[t] = _bc.n_u;
                n_v[t] = _bc.n_v;

            }else {
                n_bc[t] = false;
            }

            t++;
        }


    }

}
