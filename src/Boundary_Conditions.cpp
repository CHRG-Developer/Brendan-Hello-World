#include "Boundary_Conditions.h"
#include <stdlib.h>


Boundary_Conditions::Boundary_Conditions(int num_x_nodes, int num_y_nodes)
{
    //ctor

    int total_nodes = num_x_nodes * num_y_nodes;


      // if boundary condition is present for this cell face
      bc = new bool [total_nodes +1];
        if (n_bc==NULL) exit (1);

      rho = new double [total_nodes +1];
        if (n_rho==NULL) exit (1);


      u = new double [total_nodes +1];

      v = new double [total_nodes +1];
        if (n_v==NULL) exit (1);


        /// integer describing the boundary condition type

        // 1: Dirichlet Boundary Condition i.e. constant value at boundary
        // 2: Neumann boundary condition i.e. gradient = fixed value.
        // 3: Periodic Boundary Condition

     type_rho = new int [total_nodes +1];
        if (n_type_rho==NULL) exit (1);

     type_vel = new int [total_nodes +1];
        if (n_type_vel==NULL) exit (1);


     // node u[n] which will source u[0]
     periodic_node = new int [total_nodes +1];
        if (e_v==NULL) exit (1);
    neighbour = new int[total_nodes +1]
        if(neighbour == NULL) exit (1);
}

Boundary_Conditions::~Boundary_Conditions()
{
    //dtor
    delete [] (bc);
    bc = NULL;

    delete []  (u);
    u = NULL;
    delete []  (v);
    v = NULL;
    delete []  (rho);
    rho = NULL;


    delete [] type_rho;
    type_rho = NULL;


    delete [] type_vel;
    type_vel = NULL;


    delete [] periodic_node;
    periodic_node = NULL;

    delete [] neighbour;
    neighbour = NULL;
}

void Boundary_Conditions::assign_boundary_conditions(int num_x, int num_y, quad_bcs_plus _bc){

    /// this method should be user defined to reflect the geometry of the problem
    /// currently set up for quad problem domain-> potential overload of this operator.
    int t = 0;

    //lid driven cavity conditions
    for (int i =0; i < num_x; i++){
        for( int j=0; j < num_y; j++){

            // West boundary
            if( i ==0){
                bc[t] = true;
                rho[t] = _bc.w_rho;
                u[t] = _bc.w_u;
                v[t] = _bc.w_v;
                type_vel[t] = _bc.w_type_vel;
                type_rho[t] = _bc.w_type_rho;

                periodic_node[t] = (num_x-1) * (num_y ) + t;
                neighbour[t] = t + num_y;

            }else{
                bc[t] =false;
            }

            // east boundary
            if ( i == (num_x -1)){
                bc[t] = true;
                rho[t] = _bc.e_rho;
                u[t] = _bc.e_u;
                v[t] = _bc.e_v;
                type_vel[t] = _bc.e_type_vel;
                type_rho[t] = _bc.e_type_rho;


                periodic_node[t] = t - (num_x-1) * (num_y );
                neighbour[t] = t - num_y;


            }else {
                    bc[t] = false;
            }

            // south boundary
            if(j == 0){
                bc[t] = true;
                rho[t] = _bc.s_rho;
                u[t] = _bc.s_u;
                v[t] = _bc.s_v;
                type_vel[t] = _bc.s_type_vel;
                type_rho[t] = _bc.s_type_rho;

                periodic_node[t] = t + (num_y-1);
                neighbour[t] = t + 1;




            }else{
                bc[t] = false;
            }

            // north boundary
            if( j == (num_y-1)){
                bc[t] = true;
                rho[t] = _bc.n_rho;
                u[t] = _bc.n_u;
                v[t] = _bc.n_v;
                type_vel[t] = _bc.n_type_vel;
                type_rho[t] = _bc.n_type_rho;

                periodic_node[t] = t - (num_y-1);
                neighbour[t] = t - 1;


            }else {
                bc[t] = false;
            }

            t++;
        }


    }
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
                bc[t] = true;
                rho[t] = _bc.w_rho;
                u[t] = _bc.w_u;
                v[t] = _bc.w_v;
                type_vel[t] = _bc.w_type;
                type_rho[t] = _bc.w_type;

                periodic_node[t] = (num_x-2) * (num_y ) + t;
                neighbour[t] = t + num_y;

            }else{
                bc[t] =false;
            }

            // east boundary
            if ( i == (num_x -1)){
                bc[t] = true;
                rho[t] = _bc.e_rho;
                u[t] = _bc.e_u;
                v[t] = _bc.e_v;
                type_vel[t] = _bc.e_type;
                type_rho[t] = _bc.e_type;


                periodic_node[t] = t - (num_x-2) * (num_y );
                neighbour[t] = t - num_y;


            }else {
                    bc[t] = false;
            }

            // south boundary
            if(j == 0){
                bc[t] = true;
                rho[t] = _bc.s_rho;
                u[t] = _bc.s_u;
                v[t] = _bc.s_v;
                type_vel[t] = _bc.s_type;
                type_rho[t] = _bc.s_type;

                periodic_node[t] = t + (num_y-2);
                neighbour[t] = t + 1;




            }else{
                bc[t] = false;
            }

            // north boundary
            if( j == (num_y-1)){
                bc[t] = true;
                rho[t] = _bc.n_rho;
                u[t] = _bc.n_u;
                v[t] = _bc.n_v;
                type_vel[t] = _bc.n_type;
                type_rho[t] = _bc.n_type;

                periodic_node[t] = t - (num_y-2);
                neighbour[t] = t - 1;


            }else {
                bc[t] = false;
            }


            t++;
        }


    }

}
