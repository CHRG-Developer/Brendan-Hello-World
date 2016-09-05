#include "artificial_dissipation.h"
#include <stdlib.h>
#include <math.h>

artificial_dissipation::artificial_dissipation()
{
    //ctor
}
artificial_dissipation::artificial_dissipation(int total_nodes)
{
    //ctor
     global_JST_switch_x = new double [total_nodes +1];
        if (global_JST_switch_x==NULL) exit (1);
    global_JST_switch_y = new double [total_nodes +1];
    if (global_JST_switch_y==NULL) exit (1);
}

artificial_dissipation::~artificial_dissipation()
{
    //dtor
    delete [] (global_JST_switch_x);
    global_JST_switch_x = NULL;
    delete [] (global_JST_switch_y);
    global_JST_switch_y = NULL;
}

void artificial_dissipation::get_global_jst(Solution &soln, Boundary_Conditions &bcs,
                                            Uniform_Mesh &Mesh, domain_geometry &domain)
{
    int neighbour;

    for(int i = 0; i < Mesh.get_total_nodes(); i++){


        //first go x-direction
        // m1 = minus 1
        //p1 = plus1
        // zero = current node
        double m1,p1,zero;
         neighbour = Mesh.get_w_node(i);
       zero = soln.get_rho(i);
        if ( bcs.get_w_bc(i)){
            if(bcs.get_w_type_rho(i) == 1){

                m1 =(bcs.get_w_rho(i)-zero) * 2 + zero;

            // Neumann BC -> set constant gradient
            }else if(bcs.get_w_type_rho(i) == 2){
                /// u_bc - u_node = gradient * dx/2

                m1= bcs.get_w_rho(i) * domain.dx/2 + zero;


            }else if (bcs.get_w_type_rho(i) == 3){

                m1 = soln.get_rho( bcs.get_periodic_node(i));
            }

        }else{
            m1 = soln.get_rho(neighbour);

        }
        neighbour = Mesh.get_e_node(i);
         if ( bcs.get_e_bc(i)){
            if(bcs.get_e_type_rho(i) == 1){

                p1 =(bcs.get_e_rho(i)-zero) * 2 + zero;

            // Neumann BC -> set constant gradient
            }else if(bcs.get_e_type_rho(i) == 2){
                /// u_bc - u_node = gradient * dx/2

                p1= bcs.get_e_rho(i) * domain.dx/2 + zero;


            }else if (bcs.get_e_type_rho(i) == 3){

                p1 = soln.get_rho( bcs.get_periodic_node(i));
            }

        }else{
            p1 = soln.get_rho(neighbour);

        }

        global_JST_switch_x[i] = fabs( (m1 - 2*zero + p1)/ (m1 + 2*zero + p1));
        neighbour = Mesh.get_s_node(i);
         if ( bcs.get_s_bc(i)){
            if(bcs.get_s_type_rho(i) == 1){

                m1 =(bcs.get_s_rho(i)-zero) * 2 + zero;

            // Neumann BC -> set constant gradient
            }else if(bcs.get_s_type_rho(i) == 2){
                /// u_bc - u_node = gradient * dx/2

                m1= bcs.get_s_rho(i) * domain.dx/2 + zero;


            }else if (bcs.get_s_type_rho(i) == 3){

                m1 = soln.get_rho( bcs.get_periodic_node(i));
            }

        }else{
            m1 = soln.get_rho(neighbour);

        }
        neighbour = Mesh.get_n_node(i);
         if ( bcs.get_n_bc(i)){
            if(bcs.get_n_type_rho(i) == 1){

                p1 =(bcs.get_n_rho(i)-zero) * 2 + zero;

            // Neumann BC -> set constant gradient
            }else if(bcs.get_n_type_rho(i) == 2){
                /// u_bc - u_node = gradient * dx/2

                p1= bcs.get_n_rho(i) * domain.dx/2 + zero;


            }else if (bcs.get_n_type_rho(i) == 3){

                p1 = soln.get_rho( bcs.get_periodic_node(i));
            }

        }else{
            p1 = soln.get_rho(neighbour);

        }

        global_JST_switch_y[i] = fabs( (m1 - 2*zero + p1)/ (m1 + 2*zero + p1));

    }



}

void artificial_dissipation::set_local_jst(){

}

void artificial_dissipation::reset_local_jst_switch(){
    local_jst_switch_x = 0.0;
    local_jst_switch_y = 0.0;

}
