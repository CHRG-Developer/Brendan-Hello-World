#include "Solver.h"
#include <math.h>
#include "vector_var.h"
Solver::Solver()
{
    //ctor
}

Solver::~Solver()
{
    //dtor
}

void Solver::Uniform_Mesh_Solver(double _dt, double _dvis, Uniform_Mesh Mesh, Solution soln, Boundary_Conditions boundary_conditions)
{
    //dtor
    dt = _dt; // timestepping for streaming
    kine_viscosity = _dvis;
    c = Mesh.get_dx() / dt;
    cs = c/sqrt(3.0);
    tau = kine_viscosity + 0.5* pow(cs,2) *dt;

    //time increment on runge kutta must satisfy CFL condition
    //U delta_t/delta_x  < 1  ->> delta_t < 1/ U
    // U is equal to max theoretical velocity in solution

    double delta_t;
    delta_t = 1;
    // loop through each cell

    for (int i=0 ; i < Mesh.get_total_nodes() ; i ++) {

        vector_var cell_1, cell_2, interface_node, lattice_node, delta_u, delta_v ,delta_w,delta_rho;
        vector_var e_alpha, u_lattice,  rho_u_interface ;
        vector_var cell_normal;
        double rho_interface,feq_interface,fneq_interface;
        flux_var x_flux , y_flux;
        flux_var cell_flux;

        double interface_area;

        bc_var bc;

        int neighbour;
        cell_1.x = Mesh.get_centroid_x(i);
        cell_1.y = Mesh.get_centroid_y(i);
        cell_1.z =Mesh.get_centroid_z(i);

        double rho_lattice ;
        double feq_lattice [8];
        double lattice_weight;
        double u_magnitude;
        // loop through cell interfaces
        for (int j= 0; j <6; j++ ){
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

            //calculate slope of macro variables
            cell_2.x = Mesh.get_centroid_x(neighbour);
            cell_2.y = Mesh.get_centroid_y((neighbour));
            cell_2.z = Mesh.get_centroid_z(neighbour);

            delta_rho.Get_Gradient(soln.get_rho(i), soln.get_rho(neighbour),cell_1,cell_2 );
            delta_u.Get_Gradient(soln.get_u(i), soln.get_u(neighbour),cell_1,cell_2 );
            delta_v.Get_Gradient(soln.get_v(i), soln.get_v(neighbour),cell_1,cell_2 );

            // include w in 3d
            delta_w.x = 0;
            delta_w.y = 0;
            delta_w.z = 0; //update for 3d

            // using D2Q9 , loop through each lattice node
            for (int k =0 ; k<9; k++){
                e_alpha = get_e_alpha(k,lattice_weight);


                 //f( r- e*c*dt) relative to cell_centroid
                lattice_node.x = interface_node.x -cell_1.x - e_alpha.x * c * dt;
                lattice_node.y = interface_node.y -cell_1.y - e_alpha.y * c * dt;
                lattice_node.z = 0; // update in 3d
                // y = mx + c

                // bc present
                if ( bc.present){

                    // greater than 90 -> lattice node within domain
                    if( cell_normal.Angle_Between_Vectors(e_alpha) > M_PI/2 ){

                            // gradient is calculated between centroid and cell interface
                            // this is because of the boundary condition
                        delta_rho.Get_Gradient(soln.get_rho(i), bc.rho,cell_1,interface_node );
                        delta_u.Get_Gradient(soln.get_u(i), bc.u,cell_1,interface_node );
                        delta_v.Get_Gradient(soln.get_v(i), bc.v,cell_1,interface_node );

                        rho_lattice = delta_rho.Dot_Product(lattice_node)
                                    + soln.get_rho(i) ;

                        u_lattice.x = delta_u.Dot_Product(lattice_node)
                                        + soln.get_u(i) ;
                        u_lattice.y = delta_v.Dot_Product(lattice_node)
                                        + soln.get_v(i) ;
                        u_lattice.z = 0;

                    // less than 90 -> lattice node outside domain
                    }else{
                        rho_lattice = bc.rho;
                        u_lattice.x = bc.u;
                        u_lattice.y = bc.v;
                        u_lattice.z = 0; //update in 3d

                    }


                }else{
                     //no bc present

                    rho_lattice = delta_rho.Dot_Product(lattice_node)
                                    + soln.get_rho(i) ;

                    u_lattice.x = delta_u.Dot_Product(lattice_node)
                                    + soln.get_u(i) ;
                    u_lattice.y = delta_v.Dot_Product(lattice_node)
                                    + soln.get_v(i) ;
                    u_lattice.z = 0;

                }



                                //3d version
    //            w_lattice = (delta_w.x * lattice_node.x
    //                           + delta_w.y * lattice_node.y
    //                           + delta_w.z * lattice_node.z)
    //                            + soln.get_w(i) ;
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
            rho_u_interface.x = rho_u_interface.x /rho_interface;
            rho_u_interface.y = rho_u_interface.y /rho_interface;
            rho_u_interface.z = rho_u_interface.z /rho_interface;
            u_magnitude = rho_u_interface.Magnitude();
            for (int k =0 ; k<9; k++){


                e_alpha = get_e_alpha(k,lattice_weight);
                e_alpha = get_e_alpha(k, lattice_weight);
//
                // get feq at cell interface
                feq_interface = 1 ;
                feq_interface = feq_interface  + e_alpha.Dot_Product(u_lattice) / pow(cs,2);
                feq_interface = feq_interface  + ( pow(e_alpha.Dot_Product(u_lattice),2)  - pow((u_magnitude* cs),2) )/ (2 * pow(cs,4));
                feq_interface = feq_interface  *lattice_weight * rho_lattice;

                //get fneq at cell interface
                fneq_interface = -tau * ( feq_interface -feq_lattice[k]);


                //calculate fluxes from feq and fneq
                x_flux.P = x_flux.P + e_alpha.x * feq_interface;
                y_flux.P = y_flux.P + e_alpha.y * feq_interface;
                x_flux.Momentum_x = x_flux.Momentum_x + pow(e_alpha.x,2) *( feq_interface + (1-1/(2*tau))*fneq_interface);
                x_flux.Momentum_y = x_flux.Momentum_y + e_alpha.x*e_alpha.y *( feq_interface + (1-1/(2*tau))*fneq_interface);
                y_flux.Momentum_x = y_flux.Momentum_x + e_alpha.x*e_alpha.y *( feq_interface + (1-1/(2*tau))*fneq_interface);
                y_flux.Momentum_y = y_flux.Momentum_y + pow(e_alpha.y,2) *( feq_interface + (1-1/(2*tau))*fneq_interface);

                cell_flux.P = cell_flux.P + interface_area/ Mesh.get_cell_volume(i)* ( x_flux.P * cell_normal.x + y_flux.P *cell_normal.y );
                cell_flux.Momentum_x = cell_flux.Momentum_x + interface_area/ Mesh.get_cell_volume(i)*
                                ( x_flux.Momentum_x * cell_normal.x + y_flux.Momentum_x *cell_normal.y );
                cell_flux.Momentum_y = cell_flux.Momentum_y + interface_area/ Mesh.get_cell_volume(i)*
                                ( x_flux.Momentum_y * cell_normal.x + y_flux.Momentum_y *cell_normal.y );


            }




        }

        // forward euler method
        double f0,f1,f2,f3,f4;

        f1= soln.get_rho(i) + delta_t * cell_flux.P;
        f2 =  soln.get_rho(i) * soln.get_u(i) + (delta_t * cell_flux.Momentum_x) * f1;
        f2 = f2 /f1;
        f3 =  soln.get_rho(i) * soln.get_v(i) + (delta_t * cell_flux.Momentum_y) * f1;
        f3 = f3 /f1;
        soln.update(f1,f2,f3,0.0, 1);




    }
}

vector_var Solver::get_e_alpha(int k, double &lattice_weight ){

        vector_var temp;
        //get e_alpha again
        if (k >0 && k< 5){ //

            temp.x = cos((k-1)*M_PI/2 );
            temp.y = sin((k-1)*M_PI/2 );
            temp.z = 0; //update in 3D
            lattice_weight = 1/36;
        }else if( k >4){

            temp.x = sqrt(2) * cos((k-5)*M_PI/2 + M_PI/4 ) ;
            temp.y = sqrt(2) * sin((k-5)*M_PI/2 + M_PI/4 );
            temp.z = 0; //update in 3D
            lattice_weight = 1/9;

        }else{
            temp.x = 0 ;
            temp.y = 0;
            temp.z = 0;
            lattice_weight = 4/9;
        }

    return temp;
}


void Solver::cell_interface_variables( int j, int i, vector_var interface_node, int neighbour, double interface_area,
                              vector_var cell_normal, Boundary_Conditions boundary_conditions,  bc_var bc,
                              Uniform_Mesh Mesh) {

          switch(j) {

            case 0: // West
                interface_node.x = Mesh.get_west_x(i);
                interface_node.y = Mesh.get_west_y(i);
                interface_node.z= Mesh.get_west_z(i);
                neighbour =Mesh.get_w_node(i);
                interface_area = Mesh.get_w_area(i);
                cell_normal.x = Mesh.get_w_i(i);
                cell_normal.y = Mesh.get_w_j(i);
                cell_normal.z = Mesh.get_w_k(i);
                if ( boundary_conditions.get_w_bc(i)){
                    bc.present = true;
                    bc.rho = boundary_conditions.get_w_rho(i);
                    bc.u = boundary_conditions.get_w_u(i);
                    bc.v = boundary_conditions.get_w_v(i);


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
