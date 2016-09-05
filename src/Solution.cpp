#include <stdlib.h>
#include "Solution.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <math.h>
using namespace std;

Solution::Solution(){

}
Solution::Solution(int _total_nodes)
{
    //ctor
    total_nodes = _total_nodes;
     //rho = (double*) malloc (sizeof(double)*(total_nodes));
     rho = new double [total_nodes +1];
        if (rho==NULL) exit (1);
     u = new double [total_nodes+1];
        if (u==NULL) exit (1);
     v = new double [total_nodes +1];
        if (v==NULL) exit (1);
     w = new double [total_nodes +1];
        if (w==NULL) exit (1);
    Initialise();

}

Solution::~Solution()
{
    //dtor
    delete [] rho;
    rho = NULL;
    delete [] u;
    u = NULL;
    delete [] v;
    v= NULL;
    delete [] w;
    w= NULL;

}

void Solution::Initialise() {

    std::fill_n(rho, total_nodes , 0.00);
    std::fill_n(u, total_nodes, 0.0);
    std::fill_n(v, total_nodes , 0.0);
    std::fill_n(w, total_nodes , 0.0);

    average_rho = 0.0; //default value

}

void Solution::assign_pressure_gradient( vector_var _gradient, vector_var gradient_origin,
    vector_var origin_magnitude, Uniform_Mesh &Mesh){

   vector_var displacement;
   vector_var rho_temp;

   for( int t =0 ; t< Mesh.get_total_nodes(); t++){


            displacement.x = Mesh.get_centroid_x(t)-gradient_origin.x;
            displacement.y = Mesh.get_centroid_y(t)- gradient_origin.y;
            displacement.z = Mesh.get_centroid_z(t) - gradient_origin.z;

            rho_temp = rho_temp.line_magnitude(origin_magnitude,_gradient,displacement);
            rho[t] = rho_temp.Magnitude();

        }
    displacement.add(rho_temp) ;

   }

void Solution::assign_velocity_gradient( vector_var _gradient, vector_var gradient_origin,
    vector_var origin_magnitude, Uniform_Mesh &Mesh){

   vector_var displacement;
   vector_var vel_temp;

   for( int t =0 ; t< Mesh.get_total_nodes(); t++){


            displacement.x = Mesh.get_centroid_x(t)-gradient_origin.x;
            displacement.y = Mesh.get_centroid_y(t)- gradient_origin.y;
            displacement.z = Mesh.get_centroid_z(t) - gradient_origin.z;

            vel_temp = vel_temp.line_magnitude(origin_magnitude,_gradient,displacement);
            u[t] = vel_temp.Magnitude();

        }
    displacement.add(vel_temp) ;

   }

void Solution::update ( double _rho, double _u, double _v, double _w , int i){

    rho[i] =_rho;
    u[i] = _u;
    v[i] = _v;
    w[i] = _w;
}

void Solution::output (std::string output_location){

    std::ofstream rho_txt,u_txt,v_txt ;
    std::string rho_file, u_file, v_file;
    rho_file = output_location + "/rho.txt";
    u_file = output_location + "/u.txt";
    v_file = output_location + "/v.txt";

    rho_txt.open(rho_file.c_str(), ios::out);
    u_txt.open(u_file.c_str(), ios::out);
    v_txt.open(v_file.c_str(), ios::out);

    for( int i = 0; i < total_nodes; i++){

        rho_txt << i << " ,"  << rho[i] << endl;
        u_txt << i << " ,"  << u[i] << endl;
        v_txt << i << " ,"  << v[i] << endl;


    }

    rho_txt.close();
    u_txt.close();
    v_txt.close();

}

void Solution::clone( Solution &soln_a){

        for (int i =0; i< total_nodes; i++){
            rho[i] = soln_a.get_rho(i);
            u[i] = soln_a.get_u(i);
            v[i] = soln_a.get_v(i);
            w[i] = soln_a.get_w(i);

        }
        average_rho = soln_a.get_average_rho();
}


void Solution::post_process(double gamma){

 for (int i =0; i< total_nodes; i++){
            rho[i] = rho[i] /gamma;

        }


}


// update bc nodes to allow for changes in solution
void Solution::update_bcs(Boundary_Conditions &bcs,Uniform_Mesh &mesh,domain_geometry &domain){

    for(int i =0; i< mesh.get_total_nodes();i++){

        // if bc present
        if (bcs.get_bc(i)){

            ///NEEDS to be modified for non-uniform solver
            // 1 = dirichlet, 2 = neumann, 3 = periodic
            if(bcs.get_rho_type(i) == 1){
                soln.rho[i] = bcs.rho[i] - (soln.rho[bcs.neighbour[i] -bcs.rho[i]);

            }else if(bcs.get_rho_type(i) == 2){
                soln.rho[i] = soln.rho[bcs.neighbour[i]] + domain.dx*bcs.rho[i]

            }else if(bs.get_rho_type(i) == 3){
                soln.rho[i] = soln.rho[bcs.periodic_node[i]];

            }


            if(bcs.get_vel_type(i) == 1){
                soln.u[i] = bcs.u[i] - (soln.u[bcs.neighbour[i] -bcs.u[i]);
                soln.v[i] = bcs.v[i] - (soln.v[bcs.neighbour[i] -bcs.v[i]);
            }else if(bcs.get_vel_type(i) == 2){
                soln.u[i] = soln.u[bcs.neighbour[i]] + domain.dx*bcs.u[i]
                soln.v[i] = soln.v[bcs.neighbour[i]] + domain.dx*bcs.v[i]
            }else if(bs.get_vel_type(i) == 3){
                soln.u[i] = soln.u[bcs.periodic_node[i]];
                 soln.v[i] = soln.v[bcs.periodic_node[i]];
            }




        }

    }


}
void Solution::restriction(Solution &coarse_soln,Uniform_Mesh &coarse_mesh,
                           Uniform_Mesh &fine_mesh, Boundary_Conditions &bc){

    int coarse_x, coarse_y;
    int coarse_i;

    double coarse_rho,coarse_u,coarse_v,coarse_w;

    coarse_soln.set_average_rho( average_rho);

    //may need to swap the approach here around for parrelisation
    // i.e. loop through coarse mesh
    for (int i =0; i< total_nodes; i++){

            if(!bc.get_bc(i)){
                // get index in terms of x and y

                // 0.5 allows for ghost cells
                coarse_x = floor( (i/ fine_mesh.get_num_y()+0.5)/2.0);
                coarse_y = floor((fmod(i, fine_mesh.get_num_y()) +0.5)/2.0);

                coarse_i = coarse_mesh.get_num_y()* coarse_x + coarse_y;

                // Uniform Mesh -> get area is not needed-> just divide by 4


                // add area_fine/area_coarse * var to coarse_i


                coarse_soln.add_rho(coarse_i, rho[i]/4.0);
                coarse_soln.add_u(coarse_i, u[i]/4.0);
                coarse_soln.add_v(coarse_i, v[i]/4.0);
                coarse_soln.add_w(coarse_i,w[i]/4.0);
            }


    }


}



void Solution::prolongation(Solution &coarse_soln, Solution &temp_soln, Solution &soln,
                            Uniform_Mesh &coarse_mesh, Uniform_Mesh &fine_mesh,
                    Boundary_Conditions &bc ,bool fmg){
        double mg_delta_rho, mg_delta_u, mg_delta_v, mg_delta_w;
            //loop through the finer mesh as this will enable parrelisation later

        int edge_cell_x,edge_cell_y,vertex_cell,coarse_i,coarse_x,coarse_y;
        int fine_x, fine_y;

        double mg_factor[4] = {9.0/16.0 ,3.0/16.0, 3.0/16.0, 1./16.0 };
        double rho_edge_factor, vel_edge_factor;
        double rho_bc_contribution,u_bc_contribution,v_bc_contribution,w_bc_contribution;
        bool calculate;
        Solution debug_correction(total_nodes);

       for(int i =0; i< total_nodes; i++){



            // get index in terms of x and y
            coarse_x = floor(i/ fine_mesh.get_num_y()/2.0);
            coarse_y = floor(fmod(i, fine_mesh.get_num_y())/2.0);

            // for finer cells within a coarse cell

            for(int j = 0; j <4; j++){
                    calculate = true;
                    rho_bc_contribution = 0.0;
                    u_bc_contribution = 0.0;
                    v_bc_contribution = 0.0;
                    w_bc_contribution = 0.0;
                    rho_edge_factor = 1.0;
                    vel_edge_factor = 1.0;
                 switch(j) {

                    case 0: // nearest coarse cell
                        coarse_i = coarse_mesh.get_num_y()* coarse_x + coarse_y;
                        break;

                    case 1:
                        //North/South edge cell contribution
                        edge_cell_y = coarse_y + pow(-1.0,1+ floor(fmod(fmod(i, fine_mesh.get_num_y()),2.0)));
                        coarse_i = coarse_mesh.get_num_y() * coarse_x + edge_cell_y;
                        if( edge_cell_y < 0 || edge_cell_y > (coarse_mesh.get_num_y()-1)){
                            calculate = false;
                            coarse_i = coarse_mesh.get_num_y()* coarse_x + coarse_y;
                        }
                        break;
                    case 2:
                        //East/West edge cell contribution
                        edge_cell_x = coarse_x + pow(-1.0 ,1 + floor(fmod(floor(i/ fine_mesh.get_num_y()),2.0)));
                        coarse_i = coarse_mesh.get_num_y()* edge_cell_x + coarse_y;
                        if( edge_cell_x< 0 || edge_cell_x > (coarse_mesh.get_num_x()-1)){
                            calculate = false;
                            coarse_i = coarse_mesh.get_num_y()* coarse_x + coarse_y;
                        }
                        break;
                    case 3:
                        // Vertex coarse cell Contribution
                        coarse_i = coarse_mesh.get_num_y()* edge_cell_x + edge_cell_y;

                        if( edge_cell_x< 0 || edge_cell_x > (coarse_mesh.get_num_x()-1) ||
                            edge_cell_y< 0 || edge_cell_y > (coarse_mesh.get_num_y()-1)){
                            calculate = false;
                            coarse_i = coarse_mesh.get_num_y()* coarse_x + coarse_y;
                        }

                        break;
                }

                // check if fine cell is on corner or edge

                // get index in terms of x and y
                fine_x = floor(i/ fine_mesh.get_num_y());
                fine_y = floor(fmod(i, fine_mesh.get_num_y()));

                // West Edge




                if (fine_x == 0) {
                    //both corners on west edge

                        ///south west corner
                    if (fine_y == 0){
                                                //using barycentric interpolation
                        // coarse node provides 0.75 of contribution
                        //rho_edge_factor = 4.0/3.0;
                        //vel_edge_factor = 4.0/3.0;

                        corner_interpolation(bc.get_w_type_rho(coarse_i),bc.get_s_type_rho(coarse_i),rho_edge_factor,
                                             bc.get_w_rho(coarse_i),bc.get_s_rho(coarse_i),rho_bc_contribution,
                                             coarse_soln.get_rho(coarse_i),coarse_mesh.get_dx(),
                                             coarse_mesh.get_dy(),coarse_i,bc,coarse_soln,1,fmg);

                        corner_interpolation(bc.get_w_type_vel(coarse_i),bc.get_s_type_vel(coarse_i),vel_edge_factor,
                                             bc.get_w_u(coarse_i),bc.get_s_u(coarse_i),u_bc_contribution,
                                             coarse_soln.get_u(coarse_i),coarse_mesh.get_dx(),
                                             coarse_mesh.get_dy(),coarse_i,bc,coarse_soln,2,fmg);

                        corner_interpolation(bc.get_w_type_vel(coarse_i),bc.get_s_type_vel(coarse_i),vel_edge_factor,
                                             bc.get_w_v(coarse_i),bc.get_s_v(coarse_i),v_bc_contribution,
                                             coarse_soln.get_rho(coarse_i),coarse_mesh.get_dx(),
                                             coarse_mesh.get_dy(),coarse_i,bc,coarse_soln,3,fmg );


                        ///north west corner
                    }else if( fine_y == (fine_mesh.get_num_y()-1)){

                        //using barycentric interpolation
                        // coarse node provides 0.75 of contribution

                        corner_interpolation(bc.get_w_type_rho(coarse_i),bc.get_n_type_rho(coarse_i),rho_edge_factor,
                                             bc.get_w_rho(coarse_i),bc.get_n_rho(coarse_i),rho_bc_contribution,
                                             coarse_soln.get_rho(coarse_i),coarse_mesh.get_dx(),
                                             coarse_mesh.get_dy(),coarse_i,bc,coarse_soln,1,fmg);

                        corner_interpolation(bc.get_w_type_vel(coarse_i),bc.get_n_type_vel(coarse_i),vel_edge_factor,
                                             bc.get_w_u(coarse_i),bc.get_n_u(coarse_i),u_bc_contribution,
                                             coarse_soln.get_u(coarse_i),coarse_mesh.get_dx(),
                                             coarse_mesh.get_dy(),coarse_i,bc,coarse_soln,2,fmg);

                        corner_interpolation(bc.get_w_type_vel(coarse_i),bc.get_n_type_vel(coarse_i),vel_edge_factor,
                                             bc.get_w_v(coarse_i),bc.get_n_v(coarse_i),v_bc_contribution,
                                             coarse_soln.get_v(coarse_i),coarse_mesh.get_dx(),
                                             coarse_mesh.get_dy(),coarse_i,bc,coarse_soln,3,fmg);




                    }else {

                        ///West Edge
                        //dirichlet condition
                        if( bc.get_w_type_rho(coarse_i)  == 1){

                            rho_edge_factor = 2.0/3.0;
                            if(fmg){
                                rho_bc_contribution = bc.get_w_rho(coarse_i)* mg_factor[j] *2.0/3.0;
                            }else{
                                rho_bc_contribution =0.0;
                            }
                        }else if( bc.get_w_type_rho(coarse_i) ==3){
                            rho_bc_contribution = coarse_soln.get_rho(bc.get_periodic_node(coarse_i))
                            * mg_factor[j]/3.0;

                        }else{
                             //Neumann give constant value in x- direction
                            rho_edge_factor = 4.0/3.0;
                        }


                        if( bc.get_w_type_vel(coarse_i)  == 1){
                            vel_edge_factor = 2.0/3.0;
                            if(fmg){
                                u_bc_contribution = bc.get_w_u(coarse_i)* mg_factor[j] *2.0/3.0;
                                v_bc_contribution = bc.get_w_v(coarse_i)* mg_factor[j] *2.0/3.0;
                            }else{
                                u_bc_contribution = 0.0;
                                v_bc_contribution = 0.0;

                            }

                            //periodic node
                        }else if( bc.get_w_type_vel(coarse_i) ==3){
                            u_bc_contribution = coarse_soln.get_u(bc.get_periodic_node(coarse_i))
                            * mg_factor[j] /3.0;
                            v_bc_contribution = coarse_soln.get_v(bc.get_periodic_node(coarse_i))
                            * mg_factor[j] /3.0;

                        }else{  //Neumann and periodic give constant value in x- direction
                            vel_edge_factor = 4.0/3.0;
                        }

                    }
                }else if(fine_x == (fine_mesh.get_num_x()-1)){
                    //both corners on east edge

                        ///south east corner
                        if (fine_y == 0){
                            corner_interpolation(bc.get_e_type_rho(coarse_i),bc.get_s_type_rho(coarse_i),rho_edge_factor,
                                             bc.get_e_rho(coarse_i),bc.get_s_rho(coarse_i),rho_bc_contribution,
                                            coarse_soln.get_rho(coarse_i),coarse_mesh.get_dx(),
                                             coarse_mesh.get_dy(),coarse_i,bc,coarse_soln,1,fmg);
                            corner_interpolation(bc.get_e_type_vel(coarse_i),bc.get_s_type_vel(coarse_i),vel_edge_factor,
                                                 bc.get_e_u(coarse_i),bc.get_s_u(coarse_i),u_bc_contribution,
                                                 coarse_soln.get_u(coarse_i),coarse_mesh.get_dx(),
                                             coarse_mesh.get_dy(),coarse_i,bc,coarse_soln,2,fmg);
                            corner_interpolation(bc.get_e_type_vel(coarse_i),bc.get_s_type_vel(coarse_i),vel_edge_factor,
                                                 bc.get_e_v(coarse_i),bc.get_s_v(coarse_i),v_bc_contribution,
                                                 coarse_soln.get_v(coarse_i),coarse_mesh.get_dx(),
                                             coarse_mesh.get_dy(),coarse_i,bc,coarse_soln,3,fmg);

                            ///north east corner
                        }else if( fine_y == (fine_mesh.get_num_y()-1)){
                            corner_interpolation(bc.get_e_type_rho(coarse_i),bc.get_n_type_rho(coarse_i),rho_edge_factor,
                                             bc.get_e_rho(coarse_i),bc.get_n_rho(coarse_i),rho_bc_contribution,
                                                 coarse_soln.get_rho(coarse_i),coarse_mesh.get_dx(),
                                             coarse_mesh.get_dy(),coarse_i,bc,coarse_soln,1,fmg);
                            corner_interpolation(bc.get_e_type_vel(coarse_i),bc.get_n_type_vel(coarse_i),vel_edge_factor,
                                                 bc.get_e_u(coarse_i),bc.get_n_u(coarse_i),u_bc_contribution,
                                                 coarse_soln.get_u(coarse_i),coarse_mesh.get_dx(),
                                             coarse_mesh.get_dy(),coarse_i,bc,coarse_soln,2,fmg);
                            corner_interpolation(bc.get_e_type_vel(coarse_i),bc.get_n_type_vel(coarse_i),vel_edge_factor,
                                                 bc.get_e_v(coarse_i),bc.get_n_v(coarse_i),v_bc_contribution,
                                                 coarse_soln.get_v(coarse_i),coarse_mesh.get_dx(),
                                             coarse_mesh.get_dy(),coarse_i,bc,coarse_soln,3,fmg);

                             // east edge
                        }else{
                                ///East Edge
                            //dirichlet condition
                            if( bc.get_e_type_rho(coarse_i)  == 1){
                                rho_edge_factor = 2.0/3.0;
                                if(fmg){
                                    rho_bc_contribution = bc.get_e_rho(coarse_i)* mg_factor[j] *2.0/3.0;
                                }else{
                                    rho_bc_contribution =0.0;
                                }

                            }else if( bc.get_e_type_rho(coarse_i) ==3){
                                rho_bc_contribution = coarse_soln.get_rho(bc.get_periodic_node(coarse_i))
                                * mg_factor[j]/3.0;

                            }else{
                                 //Neumann give constant value in x- direction
                                rho_edge_factor = 4.0/3.0;
                            }


                            if( bc.get_e_type_vel(coarse_i)  == 1){
                                vel_edge_factor = 2.0/3.0;
                                if(fmg){
                                    u_bc_contribution = bc.get_e_u(coarse_i)* mg_factor[j] *2.0/3.0;
                                    v_bc_contribution = bc.get_e_v(coarse_i)* mg_factor[j] *2.0/3.0;
                                }else{
                                    u_bc_contribution = 0.0;
                                    v_bc_contribution = 0.0;
                                }

                                //periodic node
                            }else if( bc.get_e_type_vel(coarse_i) ==3){
                                u_bc_contribution = coarse_soln.get_u(bc.get_periodic_node(coarse_i))
                                * mg_factor[j]/3.0;
                                v_bc_contribution = coarse_soln.get_v(bc.get_periodic_node(coarse_i))
                                * mg_factor[j]/3.0;

                            }else{  //Neumann and periodic give constant value in x- direction
                                vel_edge_factor = 4.0/3.0;
                            }


                        }





                    ///South Edge
                }else if ( fine_y == 0){
                        ///South Edge
                        //dirichlet condition
                        if( bc.get_s_type_rho(coarse_i)  == 1){
                            rho_edge_factor = 2.0/3.0;
                            if(fmg){
                                rho_bc_contribution = bc.get_s_rho(coarse_i)* mg_factor[j] *2.0/3.0;
                            }else{
                                rho_bc_contribution =0.0;
                            }

                        }else if( bc.get_s_type_rho(coarse_i) ==3){
                            rho_bc_contribution = coarse_soln.get_rho(bc.get_periodic_node(coarse_i))
                            * mg_factor[j]/3.0;

                        }else{
                             //Neumann give constant value in x- direction
                            rho_edge_factor = 4.0/3.0;
                        }


                        if( bc.get_s_type_vel(coarse_i)  == 1){
                            vel_edge_factor = 2.0/3.0;
                            if(fmg){
                                u_bc_contribution = bc.get_s_u(coarse_i)* mg_factor[j] *2.0/3.0;
                                v_bc_contribution = bc.get_s_v(coarse_i)* mg_factor[j] *2.0/3.0;
                            }else{
                                u_bc_contribution = 0.0;
                                v_bc_contribution = 0.0;
                            }

                            //periodic node
                        }else if( bc.get_s_type_vel(coarse_i) ==3){
                            u_bc_contribution = coarse_soln.get_u(bc.get_periodic_node(coarse_i))
                            * mg_factor[j]/3.0;
                            v_bc_contribution = coarse_soln.get_v(bc.get_periodic_node(coarse_i))
                            * mg_factor[j]/3.0;

                        }else{  //Neumann and periodic give constant value in x- direction
                            vel_edge_factor = 4.0/3.0;
                        }

                    /// North Edge
                }else if ( fine_y == (fine_mesh.get_num_y()-1)){
                         ///North Edge
                        //dirichlet condition
                        if( bc.get_n_type_rho(coarse_i)  == 1){
                            rho_edge_factor = 2.0/3.0;
                            if(fmg){
                                rho_bc_contribution = bc.get_n_rho(coarse_i)* mg_factor[j] *2.0/3.0;
                            }else{
                                rho_bc_contribution = 0.0;
                            }

                        }else if( bc.get_n_type_rho(coarse_i) ==3){
                            rho_bc_contribution = coarse_soln.get_rho(bc.get_periodic_node(coarse_i))
                            * mg_factor[j]/3.0;

                        }else{
                             //Neumann give constant value in x- direction
                            rho_edge_factor = 4.0/3.0;
                        }


                        if( bc.get_n_type_vel(coarse_i)  == 1){
                            vel_edge_factor = 2.0/3.0;
                            if(fmg ){
                                u_bc_contribution = bc.get_n_u(coarse_i)* mg_factor[j] *2.0/3.0;
                                v_bc_contribution = bc.get_n_v(coarse_i)* mg_factor[j] *2.0/3.0;
                            }else{
                                u_bc_contribution = 0.0;
                                v_bc_contribution = 0.0;
                            }

                            //periodic node
                        }else if( bc.get_n_type_vel(coarse_i) ==3){
                            u_bc_contribution = coarse_soln.get_u(bc.get_periodic_node(coarse_i))
                            * mg_factor[j]/3.0;
                            v_bc_contribution = coarse_soln.get_v(bc.get_periodic_node(coarse_i))
                            * mg_factor[j]/3.0;

                        }else{  //Neumann and periodic give constant value in x- direction
                            vel_edge_factor = 4.0/3.0;
                        }
                }

                /// coarse_soln = Q2h
                /// temp_soln = Q2h_(0)
                /// soln = Qh
                if (calculate == true){
                    //_delta_rho = 1*mg_factor[j]* edge_factor;
                    mg_delta_rho = (coarse_soln.get_rho(coarse_i) - temp_soln.get_rho(coarse_i)) *mg_factor[j]* rho_edge_factor
                            + rho_bc_contribution;
                    mg_delta_u = (coarse_soln.get_u(coarse_i) - temp_soln.get_u(coarse_i))*mg_factor[j]* vel_edge_factor
                            + u_bc_contribution;
                    mg_delta_v = (coarse_soln.get_v(coarse_i) - temp_soln.get_v(coarse_i)) *mg_factor[j]* vel_edge_factor
                            + v_bc_contribution;
                    mg_delta_w = (coarse_soln.get_w(coarse_i) - temp_soln.get_w(coarse_i))*mg_factor[j]* vel_edge_factor;

                    soln.add_rho(i, mg_delta_rho ) ;
                    soln.add_u(i, mg_delta_u);
                    soln.add_v(i,mg_delta_v);
                    soln.add_w(i,mg_delta_v);
                    debug_correction.set_rho(i,mg_delta_rho);
                    debug_correction.set_u(i,mg_delta_u);
                    debug_correction.set_v(i,mg_delta_v);
                }
    }

       }
}



/// 1 indicates first boundary edge, 2 indicate second boundary edge
void Solution::corner_interpolation(int type_1, int type_2, double &factor, double value1, double value2,
                          double &edge_contribution , double coarse_value,double d1,double d2, int coarse_i,
                          Boundary_Conditions &bc, Solution &coarse_soln, int variable, bool fmg){

         // no alteration for dirichlet

         //periodic node should be passed in as value1/value2

         if(type_1 == 3){
            switch( variable){
            case 1:
                value1 =coarse_soln.get_rho(bc.get_periodic_node(coarse_i));
                break;
            case 2:
                value1 =coarse_soln.get_u(bc.get_periodic_node(coarse_i));
                break;
            case 3:
                value1 =coarse_soln.get_v(bc.get_periodic_node(coarse_i));
                break;

            }

         }
        if(type_2 == 3){
            switch( variable){
            case 1:
                value2 =coarse_soln.get_rho(bc.get_periodic_node(coarse_i));
                break;
            case 2:
                value2 =coarse_soln.get_u(bc.get_periodic_node(coarse_i));
                break;
            case 3:
                value2 =coarse_soln.get_v(bc.get_periodic_node(coarse_i));
                break;

            }

         }
        // alter values for v cycle prolongation
        //dirichlet conditions have zero error correction
        if( !fmg){
            if(type_1 ==1){
                value1 = 0;
            }
            if(type_2 == 1){
                value2 = 0;

            }

        }




         // neumann needs alteration for gradient
         //dx is passed

         if(type_1 == 2){
            value1 = value1* d1 + coarse_value;

         }
         if(type_2 == 2){

            value2 = value2*d2 + coarse_value;
         }

        //using barycentric interpolation
        //two dirichlet conditions
        if(type_1 ==1 && type_2 == 1){
            factor = 0.0;
            edge_contribution = 0.5*value1 + 0.5*value2;

        // 1 dirichlet condition -> barycentric interpolation
        }else if (type_1 ==1){
            factor = 4.0/9.0;
            edge_contribution = 4.0/8.0 *value1 + 2.0/8.0 *value2;
        // 1 dirichlet condition
        }else if (type_2 ==1){
            factor = 4.0/9.0;
            edge_contribution = 4.0/8.0 *value2 + 2.0/8.0 *value1;

        // No dirichlet condition

        }else{
            factor = 8.0/9.0;
            edge_contribution = 1.0/4.0*value1 + 1.0/4.0* value2;
        }

  }

