#include "Solver.h"
#include <math.h>
#include "vector_var.h"
#include <iostream>
#include "Solution.h"
#include <cmath>
#include <fstream>
using namespace std;
Solver::Solver()
{
    //ctor
}

Solver::~Solver()
{
    //dtor
}

void Solver::Uniform_Mesh_Solver(double _dt, double _dtau, Uniform_Mesh &Mesh, Solution &soln, Boundary_Conditions &boundary_conditions,
                                  double simulation_length, double delta_t, double _dx, double tolerance, std::string output_file)
{
    //dtor
    dt = _dt; // timestepping for streaming

    c = 1; // assume lattice spacing is equal to streaming timestep
    cs = c/sqrt(3.0);

    tau = _dtau;
    Solution temp_soln(Mesh.get_total_nodes());


    //time increment on runge kutta must satisfy CFL condition
    //U delta_t/delta_x  < 1  ->> delta_t < 1/ U
    // U is equal to max theoretical velocity in solution

    std::ofstream error_output ;
    std::string output_dir;
    output_dir = output_file +"/error.txt";

    // error_output.open("/home/brendan/Dropbox/PhD/Test Cases/Couette Flow/error.txt", ios::out);
    error_output.open(output_dir.c_str(), ios::out);


    vector_var cell_1, cell_2, interface_node, lattice_node, delta_u, delta_v ,delta_w,delta_rho;
    vector_var e_alpha, u_lattice,  rho_u_interface , u_interface;
    vector_var cell_normal;
    double rho_interface,feq_interface,fneq_interface;
    double rho_rms_num, u_rms_num, v_rms_num, w_rms_num, error;
    double rho_rms_den, u_rms_den, v_rms_den, w_rms_den;
    double rho_rms, u_rms, v_rms, w_rms;

    flux_var x_flux , y_flux;
    flux_var cell_flux;

    double interface_area;

    bc_var bc;

    int neighbour;

    double rho_lattice ;
    double feq_lattice [9], feq_int_debug[9], fneq_int_debug[9];
    double lattice_weight;
    double u_magnitude;
    double u_bc, rho_bc,  v_bc;
    int bc_node;
    //calculate timesteps

    int timesteps;
    double time;
    timesteps = ceil( simulation_length/delta_t);

    // loop through each cell
    for (int t= 0; t < timesteps; t++){


//        //loop through periodic boundary conditions
//        for ( int m =0 ; m< Mesh.get_total_nodes(); m++){
//
//            if ( boundary_conditions.get_e_type(m) == 3 || boundary_conditions.get_w_type(m) == 3
//                || boundary_conditions.get_n_type(m) == 3 || boundary_conditions.get_s_type(m) == 3){
//
//
//                    u_bc = soln.get_u( boundary_conditions.get_periodic_node(m));
//                    v_bc = soln.get_v( boundary_conditions.get_periodic_node(m));
//                    rho_bc = soln.get_rho( boundary_conditions.get_periodic_node(m));
//
//                    soln.update(rho_bc,u_bc,v_bc,0.0,m);
//
//                }
//
//        }

        rho_rms = 0;
        u_rms =0;
        v_rms = 0;
        w_rms = 0;
        u_rms_num =0;
        v_rms_num = 0;
        w_rms_num = 0;
        u_rms_den = 0;
        v_rms_den = 0;
        w_rms_den = 0;


        for (int i=0 ; i < Mesh.get_total_nodes() ; i ++) {


            interface_area = 0.0;

            cell_1.x = Mesh.get_centroid_x(i);
            cell_1.y = Mesh.get_centroid_y(i);
            cell_1.z = Mesh.get_centroid_z(i);

            cell_flux.P =0.0;
            cell_flux.Momentum_x =0.0;
            cell_flux.Momentum_y = 0.0;
            cell_flux.Momentum_z = 0.0;
            // loop through cell interfaces
            for (int j= 0; j <4; j++ ){
                bc.present = false;
                cell_interface_variables( j, i,interface_node, neighbour, interface_area,cell_normal, boundary_conditions, bc, Mesh);

                // initialise variables

                rho_interface = 0;

                rho_u_interface.x =0;
                rho_u_interface.y = 0;
                rho_u_interface.z = 0;

                x_flux.P = 0;
                x_flux.Momentum_x =0;
                x_flux.Momentum_y =0;
                x_flux.Momentum_z =0;

                y_flux.P = 0;
                y_flux.Momentum_x =0;
                y_flux.Momentum_y =0;
                y_flux.Momentum_z =0;



                // include w in 3d
                delta_w.x = 0;
                delta_w.y = 0;
                delta_w.z = 0; //update for 3d



                // using D2Q9 , loop through each lattice node
                for (int k =0 ; k<9; k++){
                    e_alpha = get_e_alpha(k,lattice_weight,c);


                     //f( r- e*c*dt) relative to cell_centroid
                    lattice_node.x = interface_node.x -cell_1.x - e_alpha.x * dt;
                    lattice_node.y = interface_node.y -cell_1.y - e_alpha.y * dt;
                    lattice_node.z = 0; // update in 3d
                    // y = mx + c

                    // bc present - only dirichlet and Neumann
                    if ( bc.present ){
                        cell_2.x = interface_node.x;
                        cell_2.y = interface_node.y;
                        cell_2.z = interface_node.z;
                        // dirichlet BC -> set correct gradient to generate correct shear stresses
                        if(bc.vel_type == 1){

                            delta_u.Get_Gradient(soln.get_u(i), bc.u,cell_1,cell_2 );
                            delta_v.Get_Gradient(soln.get_v(i), bc.v,cell_1,cell_2 );

                        // Neumann BC -> set constant gradient
                        }else if(bc.vel_type == 2){
                            /// u_bc - u_node = gradient * dx/2

                            // temp variables need to be updated for unstructure formulation
                            u_bc = bc.u * _dx/2 + soln.get_u(i);
                            v_bc = bc.v * _dx/2 + soln.get_v(i);

                            delta_u.Get_Gradient(soln.get_u(i), u_bc,cell_1,cell_2 );
                            delta_v.Get_Gradient(soln.get_v(i), v_bc,cell_1,cell_2 );

                        }else if (bc.vel_type == 3){

                            delta_u.Get_Gradient(soln.get_u(i),
                                soln.get_u( boundary_conditions.get_periodic_node(i)),cell_1,cell_2 );
                            delta_v.Get_Gradient(soln.get_v(i),
                                soln.get_v( boundary_conditions.get_periodic_node(i)),cell_1,cell_2 );
                        }

                        // dirichlet BC -> set correct gradient to generate correct shear stresses
                        if(bc.rho_type == 1){

                            delta_rho.Get_Gradient(soln.get_rho(i), bc.rho,cell_1,cell_2 );

                        // Neumann BC -> set constant gradient
                        }else if(bc.rho_type == 2){
                            /// u_bc - u_node = gradient * dx/2

                            // temp variables need to be updated for unstructure formulation
                            rho_bc = bc.rho * _dx/2 + soln.get_rho(i);
                            delta_rho.Get_Gradient(soln.get_rho(i), rho_bc,cell_1,cell_2 );

                        }else if (bc.rho_type == 3){

                            delta_rho.Get_Gradient(soln.get_rho(i),
                            soln.get_rho( boundary_conditions.get_periodic_node(i)),cell_1,cell_2 );

                        }

                    }else{

                        //calculate slope of macro variables
                        cell_2.x = Mesh.get_centroid_x(neighbour);
                        cell_2.y = Mesh.get_centroid_y((neighbour));
                        cell_2.z = Mesh.get_centroid_z(neighbour);

                        delta_rho.Get_Gradient(soln.get_rho(i), soln.get_rho(neighbour),cell_1,cell_2 );
                        delta_u.Get_Gradient(soln.get_u(i), soln.get_u(neighbour),cell_1,cell_2 );
                        delta_v.Get_Gradient(soln.get_v(i), soln.get_v(neighbour),cell_1,cell_2 );




                    }


                    rho_lattice = delta_rho.Dot_Product(lattice_node)
                                    + soln.get_rho(i) ;

                    u_lattice.x = delta_u.Dot_Product(lattice_node)
                                    + soln.get_u(i) ;
                    u_lattice.y = delta_v.Dot_Product(lattice_node)
                                    + soln.get_v(i) ;
                    u_lattice.z = 0;


                    u_magnitude = u_lattice.Magnitude();
                    feq_lattice[k] = 1 ;
                    feq_lattice[k] = feq_lattice[k] + e_alpha.Dot_Product(u_lattice) / pow(cs,2);
                    feq_lattice[k] = feq_lattice[k] + ( pow(e_alpha.Dot_Product(u_lattice),2)  - pow((u_magnitude* cs),2) )/ (2 * pow(cs,4));
                    feq_lattice[k] = feq_lattice[k] *lattice_weight * rho_lattice;

                    rho_interface = rho_interface + feq_lattice[k];
                    rho_u_interface.x = rho_u_interface.x + feq_lattice[k] * e_alpha.x;
                    rho_u_interface.y = rho_u_interface.y + feq_lattice[k] * e_alpha.y;
                    rho_u_interface.z = rho_u_interface.z + feq_lattice[k] * e_alpha.z;
                }


                // divide rho * u to get u but only after complete summation
                if( rho_interface < pow(10,-5)){
                    u_magnitude = 0;
                    u_interface.x = 0;
                    u_interface.y = 0;
                    u_interface.z = 0;
                }else{
                    u_interface.x = rho_u_interface.x /rho_interface;
                    u_interface.y = rho_u_interface.y /rho_interface;
                    u_interface.z = rho_u_interface.z /rho_interface;
                    u_magnitude = u_interface.Magnitude();
                }

                for (int k =0 ; k<9; k++){


                    e_alpha = get_e_alpha(k,lattice_weight,c);

    //
                    // get feq at cell interface
                    feq_interface = 1 ;
                    feq_interface = feq_interface  + e_alpha.Dot_Product(u_interface) / pow(cs,2);
                    feq_interface = feq_interface  + ( pow(e_alpha.Dot_Product(u_interface),2)  - pow((u_magnitude* cs),2) )/ (2 * pow(cs,4));
                    feq_interface = feq_interface  *lattice_weight * rho_interface;
                    feq_int_debug[k] = feq_interface;

                    //get fneq at cell interface
                    fneq_interface = -tau * ( feq_interface -feq_lattice[k]);
                    fneq_int_debug[k] = fneq_interface;

                    //calculate fluxes from feq and fneq
                    x_flux.P = x_flux.P + e_alpha.x * feq_interface;
                    y_flux.P = y_flux.P + e_alpha.y * feq_interface;
                    //x_flux.Momentum_x = x_flux.Momentum_x + pow(e_alpha.x,2) *( feq_interface + (1-1/(2*tau))*fneq_interface);

                    x_flux.Momentum_x = x_flux.Momentum_x + e_alpha.x * (e_alpha.x) *( feq_interface + (1-1/(2*tau))*fneq_interface);
                    x_flux.Momentum_y = x_flux.Momentum_y + e_alpha.x*(e_alpha.y) *( feq_interface + (1-1/(2*tau))*fneq_interface);


                    y_flux.Momentum_x = y_flux.Momentum_x + e_alpha.y*(e_alpha.x) *( feq_interface + (1-1/(2*tau))*fneq_interface);
                    //y_flux.Momentum_y = y_flux.Momentum_y + pow(e_alpha.y,2) *( feq_interface + (1-1/(2*tau))*fneq_interface);
                    y_flux.Momentum_y = y_flux.Momentum_y + e_alpha.y * (e_alpha.y) *( feq_interface + (1-1/(2*tau))*fneq_interface);



                }

                if( Mesh.get_cell_volume(i) < pow(10,-5)){
                    //do nothing for now

                }else{
                    cell_flux.P = 0;
                    //cell_flux.P = cell_flux.P + (-1)*interface_area/ Mesh.get_cell_volume(i)* ( x_flux.P * cell_normal.x + y_flux.P *cell_normal.y );
                    cell_flux.Momentum_x = cell_flux.Momentum_x + (-1)*interface_area/ Mesh.get_cell_volume(i)*
                                    ( x_flux.Momentum_x * cell_normal.x + y_flux.Momentum_x *cell_normal.y );
                    cell_flux.Momentum_y = cell_flux.Momentum_y + (-1)*interface_area/ Mesh.get_cell_volume(i)*
                                    ( x_flux.Momentum_y * cell_normal.x + y_flux.Momentum_y *cell_normal.y );
                }



            }

            // forward euler method
            double f1,f2,f3;
            // try removing
            f1= soln.get_rho(i) + delta_t * cell_flux.P;
            f2 =  soln.get_rho(i) * soln.get_u(i) + (delta_t * cell_flux.Momentum_x) ;
            f3 =  soln.get_rho(i) * soln.get_v(i) + (delta_t * cell_flux.Momentum_y) ;

            if (f1 < pow(10,-5)){

                 f2 = soln.get_rho(i) * soln.get_u(i);
                 f3 =soln.get_rho(i) * soln.get_v(i) ;

                }else {

                f2 = f2 /f1;
                f3 = f3 /f1;
            }

            temp_soln.update(f1,f2,f3,0.0, i);

            // error calculations

            rho_rms_num = rho_rms_num + pow( f1 -soln.get_rho(i) ,2.0 );
            if( f1 < pow(1,-8)){
                rho_rms_den = rho_rms_den + pow(1,5) ;

            }else{
                rho_rms_den = rho_rms_den + pow(f1,2) ;

            }
            u_rms_num = u_rms_num + pow( f2 -soln.get_u(i) ,2.0 );
            u_rms_den = u_rms_den + pow(f2,2) ;

            v_rms_num = v_rms_num + pow( f3 -soln.get_v(i) ,2.0 ) ;
            v_rms_den = v_rms_den +pow(f3,2);
            v_rms_den = v_rms_den + 1;


        }


        // get root of squared approximation error

        if (rho_rms_den < pow(1, -8)){
            rho_rms = 0 ;
        }else{
            rho_rms = sqrt(rho_rms_num / rho_rms_den);
        }

        if (u_rms_den < pow(1, -8)){
            u_rms = u_rms ;
        }else{
            u_rms =sqrt(u_rms_num / u_rms_den);
        }
        if (v_rms_den < pow(1, -8)){
            v_rms = v_rms;
        }else{
            v_rms = sqrt(v_rms_num / v_rms_den) ;
        }



        error = fmax(rho_rms, u_rms);
        error = fmax(error, v_rms);

        error_output << t << ", "  << error << endl;



        for (int i =0; i< Mesh.get_total_nodes(); i++){
            soln.set_rho(i,temp_soln.get_rho(i));
            soln.set_u(i,temp_soln.get_u(i));
            soln.set_v(i,temp_soln.get_v(i));
            soln.set_w(i,temp_soln.get_w(i));

        }


        time = t*delta_t;
        cout << "time t=" << time << std::endl;

        if ( error < tolerance ){

            error_output.close();
            return ;

        }
    }

    error_output.close();
}

vector_var Solver::get_e_alpha(int k, double &lattice_weight, double c ){

        vector_var temp;
        //get e_alpha again
        if (k >0 && k< 5){ //

            temp.x = cos((k-1)*M_PI/2 ) * c;
            temp.y = sin((k-1)*M_PI/2 )* c;
            temp.z = 0; //update in 3D
            lattice_weight = 1.0/9.0;
        }else if( k >4){

            temp.x = sqrt(2) * cos((k-5)*M_PI/2 + M_PI/4 ) * c ;
            temp.y = sqrt(2) * sin((k-5)*M_PI/2 + M_PI/4 ) * c;
            temp.z = 0; //update in 3D
            lattice_weight = 1.0/36.0;

        }else{
            temp.x = 0 ;
            temp.y = 0;
            temp.z = 0;
            lattice_weight = 4.0/9.0;
        }

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
                if ( boundary_conditions.get_w_bc(i)){
                    bc.present = true;
                    bc.rho = boundary_conditions.get_w_rho(i);
                    bc.u = boundary_conditions.get_w_u(i);
                    bc.v = boundary_conditions.get_w_v(i);
                    bc.vel_type = boundary_conditions.get_w_type_vel(i);
                    bc.rho_type = boundary_conditions.get_w_type_rho(i);
                    bc.periodic_node = boundary_conditions.get_periodic_node(i);

                }
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
                if ( boundary_conditions.get_s_bc(i)){
                    bc.present = true;
                    bc.rho = boundary_conditions.get_s_rho(i);
                    bc.u = boundary_conditions.get_s_u(i);
                    bc.v = boundary_conditions.get_s_v(i);
                    bc.vel_type = boundary_conditions.get_s_type_vel(i);
                    bc.rho_type = boundary_conditions.get_s_type_rho(i);
                    bc.periodic_node = boundary_conditions.get_periodic_node(i);


                }
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
                if ( boundary_conditions.get_e_bc(i)){
                    bc.present = true;
                    bc.rho = boundary_conditions.get_e_rho(i);
                    bc.u = boundary_conditions.get_e_u(i);
                    bc.v = boundary_conditions.get_e_v(i);
                    bc.vel_type = boundary_conditions.get_e_type_vel(i);
                    bc.rho_type = boundary_conditions.get_e_type_rho(i);
                    bc.periodic_node = boundary_conditions.get_periodic_node(i);


                }
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
                if ( boundary_conditions.get_n_bc(i)){
                    bc.present = true;
                    bc.rho = boundary_conditions.get_n_rho(i);
                    bc.u = boundary_conditions.get_n_u(i);
                    bc.v = boundary_conditions.get_n_v(i);
                    bc.vel_type = boundary_conditions.get_n_type_vel(i);
                    bc.rho_type = boundary_conditions.get_n_type_rho(i);
                    bc.periodic_node = boundary_conditions.get_periodic_node(i);


                }
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
