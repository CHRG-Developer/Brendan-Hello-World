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
            u[t] = vel_temp.x + vel_temp.y +vel_temp.z;

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
                rho[i] = bcs.get_rho(i) - (rho[bcs.get_neighbour(i)] -bcs.get_rho(i));

            }else if(bcs.get_rho_type(i) == 2){
                rho[i] = rho[bcs.get_neighbour(i)] + domain.dx*bcs.get_rho(i);

            }else if(bcs.get_rho_type(i) == 3){
                rho[i] = rho[bcs.get_periodic_node(i)];

            }


            if(bcs.get_vel_type(i) == 1){
                u[i] = bcs.get_u(i) - (u[bcs.get_neighbour(i)] -bcs.get_u(i));
                v[i] = bcs.get_v(i) - (v[bcs.get_neighbour(i)] -bcs.get_v(i));
            }else if(bcs.get_vel_type(i) == 2){
                u[i] = u[bcs.get_neighbour(i)] + domain.dx*bcs.get_u(i);
                v[i] = v[bcs.get_neighbour(i)] + domain.dx*bcs.get_v(i);
            }else if(bcs.get_vel_type(i) == 3){
                u[i] = u[bcs.get_periodic_node(i)];
                 v[i] = v[bcs.get_periodic_node(i)];
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
                coarse_x = floor( (i/ fine_mesh.get_num_y())/2.0 + 0.5);
                coarse_y = floor((fmod(i, fine_mesh.get_num_y()) )/2.0 +0.5);

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

        int edge_cell_x,edge_cell_y,coarse_i,coarse_x,coarse_y;
   

        double mg_factor[4] = {9.0/16.0 ,3.0/16.0, 3.0/16.0, 1./16.0 };
        
        bool calculate;
        Solution debug(fine_mesh.get_total_nodes());
        
        debug.Initialise();
        
       for(int i =0; i< total_nodes; i++){

            if(! bc.get_bc(i)){
          
            // get index in terms of x and y
            coarse_x = floor(i/ fine_mesh.get_num_y()/2.0 +0.5);
            coarse_y = floor(fmod(i, fine_mesh.get_num_y())/2.0+ 0.5);

            // for finer cells within a coarse cell

            for(int j = 0; j <4; j++){
                    calculate = true;

                 switch(j) {

                    case 0: // nearest coarse cell
                        coarse_i = coarse_mesh.get_num_y()* coarse_x + coarse_y;
                        break;

                    case 1:
                        //North/South edge cell contribution
                        edge_cell_y = coarse_y + pow(-1.0,floor(fmod(fmod(i, fine_mesh.get_num_y()),2.0)));
                        coarse_i = coarse_mesh.get_num_y() * coarse_x + edge_cell_y;

                        break;
                    case 2:
                        //East/West edge cell contribution
                        edge_cell_x = coarse_x + pow(-1.0 , floor(fmod(floor(i/ fine_mesh.get_num_y()),2.0)));
                        coarse_i = coarse_mesh.get_num_y()* edge_cell_x + coarse_y;

                        break;
                    case 3:
                        // Vertex coarse cell Contribution
                        coarse_i = coarse_mesh.get_num_y()* edge_cell_x + edge_cell_y;

                        break;
                }

                // check if fine cell is on corner or edge

                // West Edge


                /// coarse_soln = Q2h
                /// temp_soln = Q2h_(0)
                /// soln = Qh
                if (calculate == true){
                    //_delta_rho = 1*mg_factor[j]* edge_factor;
                    mg_delta_rho = (coarse_soln.get_rho(coarse_i) - temp_soln.get_rho(coarse_i)) *mg_factor[j];
                    mg_delta_u = (coarse_soln.get_u(coarse_i) - temp_soln.get_u(coarse_i))*mg_factor[j];
                    mg_delta_v = (coarse_soln.get_v(coarse_i) - temp_soln.get_v(coarse_i)) *mg_factor[j];
                    mg_delta_w = (coarse_soln.get_w(coarse_i) - temp_soln.get_w(coarse_i))*mg_factor[j];

                    soln.add_rho(i, mg_delta_rho ) ;
                    soln.add_u(i, mg_delta_u);
                    soln.add_v(i,mg_delta_v);
                    soln.add_w(i,mg_delta_v);
                    debug.add_rho(i, mg_delta_rho);
                    debug.add_u(i,mg_delta_u);
                    debug.add_v(i,mg_delta_v);
                }
            }
        }
    }
      calculate = true;  
}


