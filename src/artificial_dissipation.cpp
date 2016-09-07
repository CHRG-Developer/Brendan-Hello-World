#include "artificial_dissipation.h"
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <algorithm>

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

        if (bcs.get_bc(i)){
            global_JST_switch_x[i] = 0.0;
             global_JST_switch_y[i] = 0.0;

        }else{


        //first go x-direction
        // m1 = minus 1
        //p1 = plus1
        // zero = current node
        double m1,p1,zero;
        neighbour = Mesh.get_w_node(i);
        zero = soln.get_rho(i);
        m1 = soln.get_rho(neighbour);
        neighbour = Mesh.get_e_node(i);
        p1 = soln.get_rho(neighbour);

        global_JST_switch_x[i] = fabs( (m1 - 2*zero + p1)/ (m1 + 2*zero + p1));
        neighbour = Mesh.get_s_node(i);
        m1 = soln.get_rho(neighbour);

        neighbour = Mesh.get_n_node(i);
        p1 = soln.get_rho(neighbour);



        global_JST_switch_y[i] = fabs( (m1 - 2*zero + p1)/ (m1 + 2*zero + p1));
        }
    }



}
void artificial_dissipation::get_local_coeffs(Solution &soln, Boundary_Conditions &bcs,
                                            Uniform_Mesh &Mesh, Solution &local_soln, domain_geometry &domain)
{
    int neighbour;

    for(int i = 0; i < Mesh.get_total_nodes(); i++){

        if (bcs.get_bc(i)){
            global_2nd_order_x = 0.0;
             global_2nd_order_y = 0.0;

        }else{

        /// Get dissipation coefficients
        double m1,p1,zero,p2;

        zero = local_jst_switch_x;
        neighbour = Mesh.get_w_node(i);
        m1 = global_JST_switch_x[neighbour];
        neighbour = Mesh.get_e_node(i);
        p1 = global_JST_switch_x[neighbour];
        neighbour = Mesh.get_e_node(neighbour);
        p2 = global_JST_switch_x[neighbour];



        global_2nd_order_x = std::max(m1,std::max(zero,std::max(p1,p2))) * kappa_2;
        global_4th_order_x = std::max(0.0,(kappa_4 - global_2nd_order_y));

        zero = local_jst_switch_y;
        neighbour = Mesh.get_s_node(i);
        m1 = global_JST_switch_y[neighbour];
        neighbour = Mesh.get_n_node(i);
        p1 = global_JST_switch_y[neighbour];
        neighbour = Mesh.get_n_node(neighbour);
        p2 = global_JST_switch_y[neighbour];


        global_2nd_order_y = std::max(m1,std::max(zero,std::max(p1,p2))) * kappa_2;
        global_4th_order_y = std::max(0.0,(kappa_4 - global_2nd_order_y));


        /// Get spectral radii for euler equations (inviscid)

        /// needs allowances for pre conditioning


        /// artificial dissipation calcs to get lambda x
        neighbour = Mesh.get_e_node(i);

        delta_x = Mesh.centroid_x(neighbour) - Mesh.centroid_x(i);
        delta_y = Mesh.centroid_y(neighbour) - Mesh.centroid_y(i);


        // x direction radii
        spectral_radii[0] = fabs( local_soln.get_u(i)* delta_y - local_soln.get_v(i) *delta_x)
                + domain.cs*sqrt(pow(delta_y,2)+ pow(delta_x,2));
        spectral_radii[1] = fabs( soln.get_u(neighbour)* delta_y - soln.get_v(neighbour) *delta_x)
                + domain.cs*sqrt(pow(delta_y,2)+ pow(delta_x,2));

        neighbour = Mesh.get_n_node(i);

        delta_x = Mesh.centroid_x(neighbour) - Mesh.centroid_x(i);
        delta_y = Mesh.centroid_y(neighbour) - Mesh.centroid_y(i);

        neighbour = Mesh.get_e_node(i);

        // y direction radii
        spectral_radii[2] = fabs( local_soln.get_v(i)* delta_x - local_soln.get_u(i) *delta_y)
                + domain.cs*sqrt(pow(delta_y,2)+ pow(delta_x,2));
        spectral_radii[3] = fabs( soln.get_v(neighbour)* delta_x - soln.get_u(neighbour) *delta_y)
                + domain.cs*sqrt(pow(delta_y,2)+ pow(delta_x,2));


        // Martinelli scaling

        phi_i = 1 +pow( spectral_radii[2]/spectral_radii[0], martinelli_exponent) ;

        phi_i_p1 = 1 + pow(spectral_radii[3]/spectral_radii[1],martinelli_exponent);


        lambda_x = 0.5 * ( phi_i * spectral_radii[0] + phi_i_p1* spectral_radii[1]);


        /// artificial dissipation calcs to get lambda y



        neighbour = Mesh.get_e_node(i);

        delta_x = Mesh.centroid_x(neighbour) - Mesh.centroid_x(i);
        delta_y = Mesh.centroid_y(neighbour) - Mesh.centroid_y(i);

        neighbour = Mesh.get_n_node(i);

        // y direction radii
        spectral_radii[5] = fabs( local_soln.get_v(i)* delta_x - local_soln.get_u(i) *delta_y)
                + domain.cs*sqrt(pow(delta_y,2)+ pow(delta_x,2));
        spectral_radii[6] = fabs( soln.get_v(neighbour)* delta_x - soln.get_u(neighbour) *delta_y)
                + domain.cs*sqrt(pow(delta_y,2)+ pow(delta_x,2


        delta_x = Mesh.centroid_x(neighbour) - Mesh.centroid_x(i);
        delta_y = Mesh.centroid_y(neighbour) - Mesh.centroid_y(i);

         // x direction radii
        spectral_radii[7] = fabs( local_soln.get_u(i)* delta_y - local_soln.get_v(i) *delta_x)
                + domain.cs*sqrt(pow(delta_y,2)+ pow(delta_x,2));
        spectral_radii[8] = fabs( soln.get_u(neighbour)* delta_y - soln.get_v(neighbour) *delta_x)
                + domain.cs*sqrt(pow(delta_y,2)+ pow(delta_x,2));


        // Martinelli scaling

        phi_j = 1 +pow( spectral_radii[7]/spectral_radii[5], martinelli_exponent) ;
        phi_j_p1 = 1 + pow(spectral_radii[8]/spectral_radii[6],martinelli_exponent);

        lambda_y = 0.5 * (phi_j * spectral_radii[5] + phi_j_p1* spectral_radii[6] );

        disp_2_x = ;
        disp_2_y =
        disp_4_x = ;
        disp_4_y = ;





        }
    }



}

void artificial_dissipation::add_local_jst(int j, double rho_local, double rho_neighbour){

    if( fmod(j,2) < 1){
            jst_x_num = jst_x_num + rho_neighbour - rho_local;
            jst_x_den = jst_x_den + rho_neighbour + rho_local;
    }else{
            jst_y_num = jst_y_num + rho_neighbour - rho_local;
            jst_y_den = jst_y_den + rho_neighbour + rho_local;

    }

}

void artificial_dissipation::reset_local_jst_switch(){
    local_jst_switch_x = 0.0;
    local_jst_switch_y = 0.0;
    jst_x_num = 0.0;
    jst_x_den = 0.0;
    jst_y_num = 0.0;
    jst_y_den = 0.0;
}

