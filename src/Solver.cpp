#include <math.h>
#include <cmath>
#include "Solver.h"

#include "vector_var.h"
#include <iostream>
#include "Solution.h"
#include <fstream>
#include <global_variables.h>
#include "residuals.h"
#include <cstdio>
#include <ctime>
#include "artificial_dissipation.h"
#include <boost/math/special_functions/sign.hpp>
#include <limits>

using namespace std;
Solver::Solver()
{
    //ctor
}

Solver::~Solver()
{
    //dtor
}

void Solver::Uniform_Mesh_Solver( Uniform_Mesh &Mesh , Solution &soln, Boundary_Conditions &bcs,
                                   external_forces &source,global_variables &globals, domain_geometry &domain,
                                   initial_conditions &init_conds, quad_bcs_plus &quad_bcs_orig, int mg,
                                   Solution &residual, int fmg)
{
    //dtor
    dt = 1.0; // timestepping for streaming // non-dim equals 1

    c = 1; // assume lattice spacing is equal to streaming timestep
    cs = domain.cs;

    tau = globals.tau;
    Solution temp_soln(Mesh.get_total_nodes());
    Solution soln_t0(Mesh.get_total_nodes());
    
    
    artificial_dissipation arti_dis(Mesh.get_total_nodes(),globals);

    temp_soln.clone(soln);
    soln_t0.clone(soln);
    
//    Solution RK1(Mesh.get_total_nodes()), RK2(Mesh.get_total_nodes()), RK3(Mesh.get_total_nodes())
//        RK4(Mesh.get_total_nodes());
    
    
    flux_var RK;
    double RK_delta_t,RK_weight;
    double delta_t = globals.time_marching_step;
    std::clock_t start;
    double duration;

    // separate toelrances for convergence for FMG and normal solver
    double local_tolerance;
    if(fmg > 0){
        local_tolerance = globals.fmg_tolerance;


    }else{
        local_tolerance = globals.tolerance;
    }

    //time increment on runge kutta must satisfy CFL condition
    //U delta_t/delta_x  < 1  ->> delta_t < 1/ U
    // U is equal to max theoretical velocity in solution

    std::ofstream error_output ;
    std::string output_dir;
    if (mg == 0){

        if (fmg < 0){
                output_dir = globals.output_file +"/fmg_error.txt";
        }else{
            output_dir = globals.output_file +"/error.txt";
        }

        // error_output.open("/home/brendan/Dropbox/PhD/Test Cases/Couette Flow/error.txt", ios::out);
        error_output.open(output_dir.c_str(), ios::out);
    }

    
     std::ofstream mesh_output ;

       
    output_dir = globals.output_file +"/mesh.txt";
     

        // error_output.open("/home/brendan/Dropbox/PhD/Test Cases/Couette Flow/error.txt", ios::out);
        mesh_output.open(output_dir.c_str(), ios::out);
        for (int mesh_t= 0; mesh_t < Mesh.get_total_nodes(); mesh_t++){
             mesh_output << mesh_t << "," << Mesh.get_centroid_x(mesh_t) << "," <<
                    Mesh.get_centroid_y(mesh_t) 
                     << "," << Mesh.get_centroid_z(mesh_t) << endl;
            
            
        }
    mesh_output.close();
    
    vector_var cell_1, cell_2, interface_node, lattice_node, delta_u, delta_v ,delta_w,delta_rho;
    vector_var relative_interface;
    vector_var e_alpha, u_lattice,  rho_u_interface , u_interface;
    vector_var flux_e_alpha;
    vector_var cell_normal;
    double rho_interface,feq_interface,fneq_interface;
    residuals convergence_residual;
    flux_var x_flux , y_flux;
    flux_var cell_flux ;
    flux_var *mg_forcing_term;
    
    mg_forcing_term= new flux_var [Mesh.get_total_nodes() +1 ];
    if (mg_forcing_term==NULL) exit (1);
    flux_var debug [4] ,debug_flux[4],arti_debug [4];
    flux_var dbug [4];
    flux_var int_debug[4];
    double interface_area;

    bc_var bc;

    int neighbour;

    double rho_lattice ;
    double feq_lattice [9], feq_int_debug[9], fneq_int_debug[9];
    double u_lattice_deb[9], v_lattice[9], rho_lattice_deb[9];
    double lattice_weight;
    double u_magnitude;
    double u_bc, rho_bc,  v_bc;

    int bc_node;
    //calculate timesteps

    int timesteps;
    double time;
    if( mg ==0){
        timesteps = ceil( globals.simulation_length/delta_t);
    }else{
        timesteps = 1;
    }


    // loop through each cell
    for (int t= 0; t < timesteps; t++){
       
                                   // soln is the solution at the start of every 
                                // RK step.(rk = n) Temp_soln holds the values at end of
                                // step.(rk = n+1)
        soln_t0.clone(soln);    // soln_t0 holds solution at start of time step
                                // t= 0, rk = 0

        for( int rk=0; rk< 4; rk++){
            
            //update temp_soln boundary conditions
             temp_soln.update_bcs(bcs,Mesh,domain);
             soln.clone(temp_soln);
             
            //temp_soln.update_bcs(bcs,Mesh,domain);

            convergence_residual.reset();
            //arti_dis.get_global_jst(soln,bcs, Mesh,domain);
            for (int i=0 ; i < Mesh.get_total_nodes() ; i ++) {


                // skip if a boundary node
                if(! bcs.get_bc(i)){


                    //arti_dis.reset_local_jst_switch();

                    interface_area = 0.0;


                    cell_1.x = Mesh.get_centroid_x(i);
                    cell_1.y = Mesh.get_centroid_y(i);
                    cell_1.z = Mesh.get_centroid_z(i);
                    // add in reset function
                    cell_flux.P =0.0;
                    cell_flux.momentum_x =0.0;
                    cell_flux.momentum_y = 0.0;
                    cell_flux.momentum_z = 0.0;
                    // loop through cell interfaces
                    for (int j= 0; j <4; j++ ){
                        bc.present = false;
                        cell_interface_variables( j, i,interface_node, neighbour,
                                                 interface_area,cell_normal, bcs, bc, Mesh);

                        // initialise variables
                         // add in reset function
                        rho_interface = 0;

                        rho_u_interface.x =0;
                        rho_u_interface.y = 0;
                        rho_u_interface.z = 0;

                        x_flux.P = 0;
                        x_flux.momentum_x =0;
                        x_flux.momentum_y =0;
                        x_flux.momentum_z =0;

                        y_flux.P = 0;
                        y_flux.momentum_x =0;
                        y_flux.momentum_y =0;
                        y_flux.momentum_z =0;

                        // include w in 3d

                         // add in reset function
                        delta_w.x = 0;
                        delta_w.y = 0;
                        delta_w.z = 0; //update for 3d

                        //calculate slope of macro variables
                        cell_2.x = Mesh.get_centroid_x(neighbour);
                        cell_2.y = Mesh.get_centroid_y((neighbour));
                        cell_2.z = Mesh.get_centroid_z(neighbour);

                        // use temp soln for neighbour values as these refelct real boundary conditions
                        // temp_soln should update continuously through RK stepping
                        delta_rho.Get_Gradient(temp_soln.get_rho(i), soln.get_rho(neighbour),cell_1,cell_2 );
                        delta_u.Get_Gradient(temp_soln.get_u(i), soln.get_u(neighbour),cell_1,cell_2 );
                        delta_v.Get_Gradient(temp_soln.get_v(i), soln.get_v(neighbour),cell_1,cell_2 );


                        dbug[j].P = delta_rho.x;
                        dbug[j].momentum_x = delta_u.x;
                        dbug[j].momentum_y = delta_u.y;

                 

                        // using D2Q9 , loop through each lattice node
                        for (int k =0 ; k<9; k++){



                             /// GET change in magnitude across the lattice
                            e_alpha = get_e_alpha(k,lattice_weight,c,globals.PI);


                             //f( r- e*c*dt) relative to cell_centroid
                            lattice_node.x = interface_node.x -cell_1.x - e_alpha.x * dt;
                            lattice_node.y = interface_node.y -cell_1.y - e_alpha.y * dt;
                            lattice_node.z = 0; // update in 3d
                            // y = mx + c
                            relative_interface.x = interface_node.x -cell_1.x;
                            relative_interface.y = interface_node.y -cell_1.y;

                            // lattice densities are scalars but need to have gradients relative to cell normal

                            rho_lattice = temp_soln.get_rho(i) + delta_rho.Dot_Product(relative_interface) // rho at interface
                                +   pow(-1.0,signbit(cell_normal.x))*-1 *delta_rho.Dot_Product(e_alpha)*dt;
                                            

                            // lattice velocities are relative to the cell_normal
                            // fluxes are then calculated in the global reference from density distribution
                            //functions

                            u_lattice.x = pow(-1.0,signbit(cell_normal.x)) *(delta_u.Dot_Product(lattice_node)
                                            + temp_soln.get_u(i) );
                            u_lattice.y = pow(-1.0,signbit(cell_normal.y)) *(delta_v.Dot_Product(lattice_node)
                                            + temp_soln.get_v(i)) ;
                            u_lattice.z = 0;

                            rho_lattice_deb[k] =rho_lattice;
                            u_lattice_deb[k] = u_lattice.x;
                            v_lattice[k] = u_lattice.y;

//                            u_magnitude = u_lattice.Magnitude();
//                            feq_lattice[k] = 1.0 * rho_lattice ;
//                            feq_lattice[k] = feq_lattice[k]
//                                + e_alpha.Dot_Product(u_lattice) / pow(cs,2) * temp_soln.get_average_rho();
//                            feq_lattice[k] = feq_lattice[k]
//                                + ( pow(e_alpha.Dot_Product(u_lattice),2)  - pow((u_magnitude* cs),2) )
//                            / (2.0 * pow(cs,4)* globals.pre_conditioned_gamma) * temp_soln.get_average_rho();
//                            feq_lattice[k] = feq_lattice[k] *lattice_weight ;
                            
                            u_magnitude = u_lattice.Magnitude();
                            feq_lattice[k] = 1.0  ;
                            feq_lattice[k] = feq_lattice[k]
                                + e_alpha.Dot_Product(u_lattice) / pow(cs,2) ;
                            feq_lattice[k] = feq_lattice[k]
                                + ( pow(e_alpha.Dot_Product(u_lattice),2)  - pow((u_magnitude* cs),2) )
                            / (2.0 * pow(cs,4)* globals.pre_conditioned_gamma) ;
                            feq_lattice[k] = feq_lattice[k] *lattice_weight *rho_lattice ;
                            
                            
                            
                            
                            rho_interface = rho_interface + feq_lattice[k];
                            rho_u_interface.x = rho_u_interface.x + feq_lattice[k] * e_alpha.x;
                            rho_u_interface.y = rho_u_interface.y + feq_lattice[k] * e_alpha.y;
                            rho_u_interface.z = rho_u_interface.z + feq_lattice[k] * e_alpha.z;
                        }
                        
                        //rho_interface = rho_interface + temp_soln.get_rho(i);
                        
                        // divide rho * u to get u but only after complete summation

                            u_interface.x = rho_u_interface.x /rho_interface;
                            u_interface.y = rho_u_interface.y /rho_interface;
                            u_interface.z = rho_u_interface.z /rho_interface;
//                        u_interface.x = rho_u_interface.x /temp_soln.get_average_rho();
//                        u_interface.y = rho_u_interface.y /temp_soln.get_average_rho();
//                        u_interface.z = rho_u_interface.z /temp_soln.get_average_rho();
                        u_magnitude = u_interface.Magnitude();
                        
                   

                        int_debug[j].momentum_x =u_interface.x;
                        int_debug[j].momentum_y = u_interface.y;
                        int_debug[j].P = rho_interface;
                           
                        for (int k =0 ; k<9; k++){


                            e_alpha = get_e_alpha(k,lattice_weight,c,globals.PI);

            //
                            // get feq at cell interface
                            
                            // different order of adding may aid rounding errors
                            
//                            feq_interface = e_alpha.Dot_Product(u_interface) / pow(cs,2) *temp_soln.get_average_rho();
//                            feq_interface = feq_interface
//                                + ( pow(e_alpha.Dot_Product(u_interface),2)  - pow((u_magnitude* cs),2) )
//                                    / (2 * pow(cs,4) * globals.pre_conditioned_gamma)
//                                    *temp_soln.get_average_rho();
//                            feq_interface = feq_interface + 1 * rho_interface;
                            
                            feq_interface = e_alpha.Dot_Product(u_interface) / pow(cs,2) ;
                            feq_interface = feq_interface
                                + ( pow(e_alpha.Dot_Product(u_interface),2)  - pow((u_magnitude* cs),2) )
                                    / (2 * pow(cs,4) * globals.pre_conditioned_gamma)
                                    ;
                            feq_interface = (feq_interface + 1 )* rho_interface;
                           
                            
                            
                            feq_interface = feq_interface  *lattice_weight ;
                            feq_int_debug[k] = feq_interface;

                            //get fneq at cell interface
                            fneq_interface = -tau * ( feq_interface -feq_lattice[k]);
                            fneq_int_debug[k] = fneq_interface;

                            //calculate fluxes from feq and fneq

                            //as DEnsity flux is a scaler, we need to reverse local coordinates
                            // into global coordinates

                            // where as momentume flux is the dot product of e and returns a scaler

                            x_flux.P = x_flux.P + e_alpha.x* feq_interface*cell_normal.x;
                            y_flux.P = y_flux.P +  e_alpha.y * feq_interface*cell_normal.y;
                            //x_flux.momentum_x = x_flux.momentum_x + pow(e_alpha.x,2) *( feq_interface +
                                                                             //(1-1/(2*tau))*fneq_interface);

                            x_flux.momentum_x = x_flux.momentum_x + e_alpha.x * (e_alpha.x) *( feq_interface
                                    + (1-1/(2*tau))*fneq_interface);
                            x_flux.momentum_y = x_flux.momentum_y + e_alpha.x*(e_alpha.y) *( feq_interface
                                                    + (1-1/(2*tau))*fneq_interface);

                            y_flux.momentum_x = y_flux.momentum_x + e_alpha.y*(e_alpha.x) *( feq_interface
                                                    + (1-1/(2*tau))*fneq_interface);
                            //y_flux.momentum_y = y_flux.momentum_y + pow(e_alpha.y,2) *( feq_interface
                                                                             //+ (1-1/(2*tau))*fneq_interface);
                            y_flux.momentum_y = y_flux.momentum_y + e_alpha.y * (e_alpha.y) *( feq_interface
                                                    + (1-1/(2*tau))*fneq_interface);




                        }
                            //truncate_flux(x_flux);
                            //truncate_flux(y_flux);

                            debug_flux[j].P = x_flux.P  ;
                            debug_flux[j].momentum_x = x_flux.momentum_x* cell_normal.x ;
                            debug_flux[j].momentum_y = x_flux.momentum_y* cell_normal.x ;
                            debug_flux[j].momentum_z = x_flux.P* cell_normal.x ;

                       // account for rounding errors

                        //artificial dissipation calcs

                        //arti_dis.get_local_coeffs( soln,bcs,Mesh,temp_soln,domain,j,i);

                        
                            //cell_flux.P = 0;
                            cell_flux.P = cell_flux.P
                                + (-1)*interface_area*
                                ( (x_flux.P +arti_dis.local_flux.P) * cell_normal.x
                                 + (y_flux.P +arti_dis.local_flux.P) *cell_normal.y );
                            cell_flux.momentum_x = cell_flux.momentum_x
                                + (-1)*interface_area*
                                    ( (x_flux.momentum_x + arti_dis.local_flux.momentum_x)* cell_normal.x
                                     + (y_flux.momentum_x + arti_dis.local_flux.momentum_x) *cell_normal.y ) ;
                            cell_flux.momentum_y = cell_flux.momentum_y
                                + (-1)*interface_area*
                                    ( (x_flux.momentum_y +arti_dis.local_flux.momentum_y ) * cell_normal.x
                                       + ( y_flux.momentum_y +arti_dis.local_flux.momentum_y) *cell_normal.y );

                            /// debug

                                debug[j].P = (-1)*interface_area/ 1*
                                         ( x_flux.P * cell_normal.x + y_flux.P *cell_normal.y );
                                debug[j].momentum_x = (-1)*interface_area/ 1*
                                        ( x_flux.momentum_x * cell_normal.x + y_flux.momentum_x *cell_normal.y ) ;
                                debug[j].momentum_y =  (-1)*interface_area/ 1*
                                        ( x_flux.momentum_y * cell_normal.x + y_flux.momentum_y *cell_normal.y );

                                arti_debug[j].P =  arti_dis.local_flux.P * cell_normal.x
                                 + arti_dis.local_flux.P *cell_normal.y ;
                                 arti_debug[j].momentum_x =  arti_dis.local_flux.momentum_x * cell_normal.x
                                 + arti_dis.local_flux.momentum_x *cell_normal.y ;
                                 arti_debug[j].momentum_x =  arti_dis.local_flux.momentum_x * cell_normal.x
                                 + arti_dis.local_flux.momentum_x *cell_normal.y ;


                        

                    }
                    // divide sum of fluxes by zero    
                    //cell_flux.div_volume(Mesh.get_cell_volume(i),globals.small_number);
                    
                    //

                    // account for rounding errors
                    //truncate_flux(cell_flux);
                    //(cell_flux);

                    // store RK fluxes
                    if (rk == 0){
                        RK = cell_flux;
                        if( mg > 0){
                            if(t == 0){
                                mg_forcing_term[i].P =  -2 *residual.get_rho(i) + cell_flux.P;
                                mg_forcing_term[i].momentum_x = -2* residual.get_u(i) + cell_flux.momentum_x;
                                mg_forcing_term[i].momentum_y =  -2 *residual.get_v(i) + cell_flux.momentum_y;
                            }
                            //add momentum z later


                        }else{
                            mg_forcing_term[i].P = 0;
                            mg_forcing_term[i].momentum_x = 0;
                            mg_forcing_term[i].momentum_y = 0;

                        }
                        // timestep for calculating second step
                        RK_delta_t = delta_t/2;
                        RK_weight = 1.0/6.0;
                        
                    }else if (rk == 1){
                        RK = cell_flux;
                        RK_delta_t = delta_t/2;
                        RK_weight = 2.0/6.0;
                    
                    
                    }else if (rk == 2){
                        RK = cell_flux;
                        RK_delta_t = delta_t;
                        RK_weight = 2.0/6.0;
                    }else{
                        RK = cell_flux;
                    }
                    double f1,f2,f3,temp_force, R1,R2,R3;
                   
                    // old forward euler method adapted to RK4

                    temp_force = source.get_force(i)* Mesh.get_cell_volume(i) * temp_soln.get_rho(i);
                    // try removing
                    f1= soln_t0.get_rho(i) + RK_delta_t * (RK.P/Mesh.get_cell_volume(i) + mg_forcing_term[i].P);
//                    f2 =  soln_t0.get_average_rho()* soln_t0.get_u(i) + (RK_delta_t *
//                            (RK.momentum_x + mg_forcing_term[i].momentum_x +
//                            source.get_force(i)* Mesh.get_cell_volume(i) * temp_soln.get_average_rho()));
//                    f3 =  soln_t0.get_average_rho() * soln_t0.get_v(i) + (RK_delta_t *
//                                                (RK.momentum_y + mg_forcing_term[i].momentum_y)) ;
                    
                    f2 =  soln_t0.get_rho(i)* soln_t0.get_u(i) + (RK_delta_t *
                            (RK.momentum_x/Mesh.get_cell_volume(i) + mg_forcing_term[i].momentum_x +
                            source.get_force(i)* Mesh.get_cell_volume(i) * temp_soln.get_rho(i)));
                    f3 =  soln_t0.get_rho(i) * soln_t0.get_v(i) + (RK_delta_t *
                                                (RK.momentum_y/Mesh.get_cell_volume(i) + mg_forcing_term[i].momentum_y)) ;
                    
                    //f3 = 0.0;
                    switch (rk){

                        case 0:
                            residual.set_rho(i,((RK.P + mg_forcing_term[i].P)*(RK_weight)));
                            residual.set_u(i,((RK.momentum_x + mg_forcing_term[i].momentum_x)*RK_weight));
                            residual.set_v(i,((RK.momentum_y + mg_forcing_term[i].momentum_z)*RK_weight));

                            break;

                        case 1:
                        case 2:  
                        case 3:
                            residual.add_rho(i,((RK.P + mg_forcing_term[i].P)*(RK_weight)));
                            residual.add_u(i,((RK.momentum_x + mg_forcing_term[i].momentum_x)*RK_weight));
                            residual.add_v(i,((RK.momentum_y + mg_forcing_term[i].momentum_z)*RK_weight));


                            break;

                    }
                        
                    // on final loop get Runge Kutta integration
                    if( rk == 3){
                       

                        f1 = soln_t0.get_rho(i) + residual.get_rho(i)/Mesh.get_cell_volume(i)* delta_t;
//                        f2 = soln_t0.get_average_rho()* soln_t0.get_u(i) 
//                                + residual.get_u(i) *  delta_t;
//                        f3 = soln_t0.get_average_rho() * soln_t0.get_v(i) 
//                                + residual.get_v(i) * delta_t;
                        
                         f2 = soln_t0.get_rho(i)* soln_t0.get_u(i) 
                                + residual.get_u(i)/Mesh.get_cell_volume(i) *  delta_t;
                        f3 = soln_t0.get_rho(i) * soln_t0.get_v(i) 
                                + residual.get_v(i)/Mesh.get_cell_volume(i) * delta_t;
                        
                        //f3 = 0.0;

                    }




//                    f2 = f2 /soln_t0.get_average_rho();
//                    f3 = f3 /soln_t0.get_average_rho();
//                    
                       f2 = f2 /soln_t0.get_rho(i);
                    f3 = f3 /soln_t0.get_rho(i);

                    temp_soln.update(f1,f2,f3,0.0, i);

                    // error calculations
                    if( rk ==3){
    //                        residual.add_l2_norm_residuals(f1,soln.get_rho(i),f2,soln.get_u(i),
    //                            f3,soln.get_v(i));
                            convergence_residual.add_ansys_l2_norm_residuals(f1,soln.get_rho(i),f2,soln.get_u(i),
                                f3,soln.get_v(i));

                          if (std::isnan(temp_soln.get_rho(i)) || std::isnan(temp_soln.get_u(i))) {
                                if( mg == 0){
                                    error_output.close();
                                }
                                return;
                        }
                    }

                }
            }
        }

        int cycle_no;




        // get root of squared approximation error

//        residual.l2_norm_rms();
        convergence_residual.ansys_5_iter_rms(t);

        if( mg == 0){
            error_output << t << ", "  << convergence_residual.max_error()   << ", " <<
            convergence_residual.rho_rms << ", " << convergence_residual.u_rms << ", " <<
            convergence_residual.v_rms << " , FMG cycle: " << fmg << endl;
        }
        soln.clone(temp_soln);

        time = t*delta_t;



        cout << "time t=" << time << std::endl;


        if ( convergence_residual.max_error() < local_tolerance){
            if( mg == 0){
                error_output.close();
            }
            return ;

        }



        //Give the program 5 iterations to find residuals
        //without MG and also perform MG cycle every 3 iterations


        ///Multigrid Agglomeration/Prolongation
        if( mg < globals.max_mg_levels && fmg == 0 && t > 3){
                multi_grid_agglomoration(residual,soln,cycle_no,Mesh,quad_bcs_orig,
                                         init_conds,mg,globals,domain,bcs);
        }else if( mg ==0){

        }

    }

    error_output.close();
    delete []((mg_forcing_term));
    mg_forcing_term = NULL;
} 


void Solver::multi_grid_agglomoration( Solution &residual , Solution &soln, int cycle_no,
                                      Uniform_Mesh &fine_mesh, quad_bcs_plus &bcs,
                                      initial_conditions &initial_conds, int &mg,
                                      global_variables globals, domain_geometry &fine_domain,
                                      Boundary_Conditions &fine_bc){


    mg = mg +1; // increase multigrid levels by 1
    int fmg = 0;
    // create new coarse Mesh with double up dimensions
    domain_geometry coarse_domain = fine_mesh.create_coarse_mesh_domain();

    Uniform_Mesh coarse_mesh (coarse_domain);

    //globals.update_tau(coarse_domain);
    //globals.update_coarse_tau();
    //Apply Boundary Conditions to new Mesh

    Boundary_Conditions bc(coarse_mesh.get_num_x(), coarse_mesh.get_num_y());
    bc.assign_boundary_conditions(coarse_mesh.get_num_x(), coarse_mesh.get_num_y(),bcs);

    external_forces source_term(coarse_mesh.get_total_nodes());
    source_term.set_uniform_force(initial_conds.pressure_gradient); //pressure force applied through external force

    Solution coarse_soln( coarse_mesh.get_total_nodes());
    Solution coarse_residual_0( coarse_mesh.get_total_nodes()) ;

    soln.restriction(coarse_soln,coarse_mesh,fine_mesh,fine_bc);

    residual.restriction(coarse_residual_0,coarse_mesh,fine_mesh,fine_bc);

    // call solver for coarser mesh

    Solver solve_coarse;
    Solution temp_soln(coarse_mesh.get_total_nodes());
    coarse_soln.update_bcs(bc,coarse_mesh,coarse_domain);
    temp_soln.clone(coarse_soln);

    solve_coarse.Uniform_Mesh_Solver(coarse_mesh,coarse_soln,bc,source_term,globals,
                                     coarse_domain,initial_conds,bcs,mg,coarse_residual_0,fmg);

    coarse_soln.update_bcs(bc,coarse_mesh,coarse_domain);
    temp_soln.update_bcs(bc,coarse_mesh,coarse_domain);
    // prolongation occurs here
    soln.prolongation(coarse_soln, temp_soln, soln,coarse_mesh, fine_mesh,fine_bc,false);
    soln.set_average_rho(initial_conds.average_rho);
    //globals.update_fine_tau();
    mg = mg -1;

}

void Solver::truncate_flux(flux_var &flux){
    
    double temp = std::numeric_limits<double>::epsilon();
    int digits;
    digits = 10;

    if(fabs(flux.P -0.0) < temp){
        flux.P  = 0.0;
    }
    
    

    double f = pow(10,digits);
   

    temp = floor(fabs(flux.P *f))/f;
    if (flux.P < 0.0){
        flux.P = -temp;
    }else{
        flux.P = temp;
    }
    temp = floor(fabs(flux.momentum_x *f))/f;
    if (flux.momentum_x < 0.0){
        flux.momentum_x = -temp;
    }else{
        flux.momentum_x = temp;
    }
    temp = floor(fabs(flux.momentum_y *f))/f;
    if (flux.momentum_y < 0.0){
        flux.momentum_y = -temp;
    }else{
        flux.momentum_y = temp;
    }
    temp = floor(fabs(flux.momentum_z *f))/f;
    if (flux.momentum_z < 0.0){
        flux.momentum_z = -temp;
    }else{
        flux.momentum_z = temp;
    }

}

void Solver::truncate_flux(double &val){
    
    double temp = std::numeric_limits<double>::epsilon();
    int digits;
    digits = 10;

    if(fabs(val-0.0) < temp){
        val = 0.0;
    }
    
    
    double f = pow(10,digits);
   
    
    temp = floor(fabs(val*f))/f;
    if (val < 0.0){
        val = -temp;
    }else{
        val= temp;
    }


}


vector_var Solver::get_e_alpha(int k, double &lattice_weight, double c, double PI ){

        vector_var temp;
        int x ,y,z;
        //get e_alpha again
        if (k >0 && k< 5){ //

            x = round(cos((k-1)*PI/2 ) * c );
            y = round(sin((k-1)*PI/2 )* c);
            z = 0; //update in 3D
            lattice_weight = 1.0/9.0;
        }else if( k >4){

            x = round(sqrt(2) * cos((k-5)*PI/2 + PI/4 ) * c );
            y = round(sqrt(2) * sin((k-5)*PI/2 + PI/4 ) * c);
            z = 0; //update in 3D
            lattice_weight = 1.0/36.0;

        }else{
            x = 0 ;
            y = 0;
            z = 0;
            lattice_weight = 4.0/9.0;
        }
        temp.x = x;
        temp.y = y;
        temp.z = z;
        

    return temp;
}


void Solver::cell_interface_variables( int j, int i, vector_var &interface_node, int &neighbour, double &interface_area,
                              vector_var &cell_normal, Boundary_Conditions &boundary_conditions,  bc_var &bc,
                              Uniform_Mesh &Mesh) {

          switch(j) {

            case 0: // West
                interface_node.x = Mesh.get_west_x(i);
                interface_node.y = Mesh.get_west_y(i);
                interface_node.z= Mesh.get_west_z(i);
                neighbour = Mesh.get_w_node(i);
                interface_area = Mesh.get_w_area(i);
                cell_normal.x = Mesh.get_w_i(i);
                cell_normal.y = Mesh.get_w_j(i);
                cell_normal.z = Mesh.get_w_k(i);
                break;

            case 1: // South
                interface_node.x = Mesh.get_south_x(i);
                interface_node.y = Mesh.get_south_y(i);
                interface_node.z= Mesh.get_south_z(i);
                neighbour =Mesh.get_s_node(i);
                interface_area = Mesh.get_s_area(i);
                cell_normal.x = Mesh.get_s_i(i);
                cell_normal.y = Mesh.get_s_j(i);
                cell_normal.z = Mesh.get_s_k(i);

                break;
            case 2: // East
                interface_node.x = Mesh.get_east_x(i);
                interface_node.y = Mesh.get_east_y(i);
                interface_node.z= Mesh.get_east_z(i);
                interface_area = Mesh.get_e_area(i);
                neighbour =Mesh.get_e_node(i);
                cell_normal.x = Mesh.get_e_i(i);
                cell_normal.y = Mesh.get_e_j(i);
                cell_normal.z = Mesh.get_e_k(i);

                break;
            case 3: // North
                interface_node.x = Mesh.get_north_x(i);
                interface_node.y = Mesh.get_north_y(i);
                interface_node.z= Mesh.get_north_z(i);
                neighbour =Mesh.get_n_node(i);
                interface_area = Mesh.get_n_area(i);
                cell_normal.x = Mesh.get_n_i(i);
                cell_normal.y = Mesh.get_n_j(i);
                cell_normal.z = Mesh.get_n_k(i);

                break;
            case 4: // Front

                break;
            case 5: // Back
                break;

            }


      }



//                //get position of lattice node relative to cell_centroid
//                if (k >0 && k< 5){ //
//
//                    e_alpha.x = cos((k-1)*M_PI/2 );
//                    e_alpha.y = sin((k-1)*M_PI/2 );
//                    e_alpha.z = 0; //update in 3D
//                    lattice_weight = 1/36;
//                }else if( k >4){
//
//                    e_alpha.x = sqrt(2) * cos((k-5)*M_PI/2 + M_PI/4 ) ;
//                    e_alpha.y = sqrt(2) * sin((k-5)*M_PI/2 + M_PI/4 );
//                    e_alpha.z = 0; //update in 3D
//                    lattice_weight = 1/9;
//
//                }else{
//                    e_alpha.x = 0 ;
//                    e_alpha.y = 0;
//                    e_alpha.z = 0;
//                    lattice_weight = 4/9;
//                }
