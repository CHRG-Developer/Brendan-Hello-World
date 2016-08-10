#include "Uniform_Mesh.h"
#include "Boundary_Conditions.h"
#include "Solution.h"
#include "vector_var.h"
#include "external_forces.h"
#include "global_variables.h"
#include "initial_conditions.h"

#ifndef SOLVER_H
#define SOLVER_H


class Solver
{
    public:
        Solver();
        virtual ~Solver();
         void Uniform_Mesh_Solver( Uniform_Mesh &Mesh , Solution &soln, Boundary_Conditions &bc,
                                   external_forces &source,global_variables &globals, domain_geometry &domain,
                                   initial_conditions &init_conds,quad_bcs_plus &quad_bcs_orig, int mg,
                                   Solution &residual,int fmg);

        void multi_grid_agglomoration( Solution &residuals , Solution &soln,
                                         int cycle_no, Uniform_Mesh &fine_mesh,  quad_bcs_plus &bcs,
                                         initial_conditions &init_conds, int &mg, global_variables globals,
                                         domain_geometry &fine_domain);
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
        vector_var get_e_alpha(int k, double &lattice_weight, double c,double PI );

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
            int vel_type;
            int rho_type;
            int periodic_node;

        };
        void cell_interface_variables( int j,int i ,vector_var &interface_node, int &neighbour, double &interface_area,
                              vector_var &cell_normal, Boundary_Conditions &boundary_conditions,  bc_var &bc ,Uniform_Mesh &Mesh) ;
};

#endif // SOLVER_H
