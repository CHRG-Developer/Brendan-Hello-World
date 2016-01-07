#include "Uniform_Mesh.h"
#include "Boundary_Conditions.h"
#include "Solution.h"
#include "vector_var.h"

#ifndef SOLVER_H
#define SOLVER_H


class Solver
{
    public:
        Solver();
        virtual ~Solver();
         void Uniform_Mesh_Solver(double _dt,double _vis, Uniform_Mesh Mesh , Solution soln, Boundary_Conditions bc);
    protected:
    private:

        double dt;
        double tau;
        double kine_viscosity;
        double c,cs;

//        struct vector_var {
//            double x;
//            double y;
//            double z;
//        };
        vector_var get_e_alpha(int k, double &lattice_weight );

        struct flux_var {
            double P;
            double Momentum_x;
            double Momentum_y;
            double Momentum_z;
        };

        struct bc_var{

            bool present;
            double rho;
            double u;
            double v;

        };
        void cell_interface_variables( int j,int i ,vector_var &interface_node, int &neighbour, double &interface_area,
                              vector_var &cell_normal, Boundary_Conditions boundary_conditions,  bc_var &bc ,Uniform_Mesh Mesh) ;
};

#endif // SOLVER_H
