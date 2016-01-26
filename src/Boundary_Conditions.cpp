#include "Boundary_Conditions.h"
#include <stdlib.h>


Boundary_Conditions::Boundary_Conditions(int num_x_nodes, int num_y_nodes)
{
    //ctor

    int total_nodes = num_x_nodes * num_y_nodes;


      // if boundary condition is present for this cell face
      n_bc = new bool [total_nodes +1];
        if (n_bc==NULL) exit (1);
      s_bc = new bool [total_nodes +1];
        if (s_bc==NULL) exit (1);
      w_bc = new bool [total_nodes +1];
        if (w_bc==NULL) exit (1);
      e_bc = new bool [total_nodes +1];
        if (e_bc==NULL) exit (1);

     n_rho = new double [total_nodes +1];
        if (n_rho==NULL) exit (1);
      s_rho = new double [total_nodes +1];
        if (s_rho==NULL) exit (1);
      w_rho = new double [total_nodes +1];
        if (w_rho==NULL) exit (1);
      e_rho = new double [total_nodes +1];
        if (e_rho==NULL) exit (1);

      n_u = new double [total_nodes +1];
        if (n_u==NULL) exit (1);
      s_u = new double [total_nodes +1];
        if (s_u==NULL) exit (1);
      w_u = new double [total_nodes +1];
        if (w_u==NULL) exit (1);
      e_u = new double [total_nodes +1];
        if (e_u==NULL) exit (1);

      n_v = new double [total_nodes +1];
        if (n_v==NULL) exit (1);
      s_v = new double [total_nodes +1];
        if (s_v==NULL) exit (1);
      w_v = new double [total_nodes +1];
        if (w_v==NULL) exit (1);
      e_v = new double [total_nodes +1];
        if (e_v==NULL) exit (1);

        /// integer describing the boundary condition type

        // 1: Dirichlet Boundary Condition i.e. constant value at boundary
        // 2: Neumann boundary condition i.e. gradient = fixed value.
        // 3: Periodic Boundary Condition

     n_type = new int [total_nodes +1];
        if (n_v==NULL) exit (1);
      s_type = new int [total_nodes +1];
        if (s_v==NULL) exit (1);
      w_type = new int [total_nodes +1];
        if (w_v==NULL) exit (1);
      e_type = new int [total_nodes +1];
        if (e_v==NULL) exit (1);

     // node u[n] which will source u[0]
     periodic_node = new int [total_nodes +1];
        if (e_v==NULL) exit (1);
}

Boundary_Conditions::~Boundary_Conditions()
{
    //dtor
    delete [] (n_bc);
    n_bc = NULL;
    delete []  (s_bc);
    s_bc = NULL;
    delete []  (w_bc);
    w_bc = NULL;
    delete []  (e_bc);
    e_bc = NULL;
    delete []  (n_u);
    n_u = NULL;
    delete []  (n_v);
    n_v = NULL;
    delete []  (n_rho);
    n_rho = NULL;
    delete []  (e_u);
    e_u = NULL;
    delete []  (e_v);
    e_v = NULL;
    delete []  (e_rho);
    e_rho = NULL;
    delete []  (w_u);
    w_u = NULL;
    delete []  (w_v);
    w_v = NULL;
    delete []  (w_rho);
    w_rho = NULL;
    delete []  (s_u);
    s_u = NULL;
    delete []  (s_v);
    s_v = NULL;
    delete []  (s_rho);
    s_rho = NULL;

    delete [] n_type;
    n_type = NULL;
    delete [] e_type;
    e_type = NULL;
    delete [] w_type;
    w_type = NULL;
    delete [] s_type;
    s_type = NULL;

    delete [] periodic_node;
    periodic_node = NULL;
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
                w_type[t] = _bc.w_type;

                if (_bc.w_type == 3){
                    periodic_node[t] = (num_x-1) * (num_y ) + t;
                }

            }else{
                w_bc[t] =false;
            }

            // east boundary
            if ( i == (num_x -1)){
                e_bc[t] = true;
                e_rho[t] = _bc.e_rho;
                e_u[t] = _bc.e_u;
                e_v[t] = _bc.e_v;
                e_type[t] = _bc.e_type;

                if (_bc.e_type == 3){
                    periodic_node[t] = t - (num_x-1) * (num_y );
                }


            }else {
                    e_bc[t] = false;
            }

            // south boundary
            if(j == 0){
                s_bc[t] = true;
                s_rho[t] = _bc.s_rho;
                s_u[t] = _bc.s_u;
                s_v[t] = _bc.s_v;
                s_type[t] = _bc.s_type;

                if (_bc.s_type == 3){
                    periodic_node[t] = t + (num_y-1);
                }


            }else{
                s_bc[t] = false;
            }

            // north boundary
            if( j == (num_y-1)){
                n_bc[t] = true;
                n_rho[t] = _bc.n_rho;
                n_u[t] = _bc.n_u;
                n_v[t] = _bc.n_v;
                n_type[t] = _bc.n_type;

                 if (_bc.n_type == 3){
                    periodic_node[t] = t - (num_y-1);
                }


            }else {
                n_bc[t] = false;
            }

            t++;
        }


    }

}
