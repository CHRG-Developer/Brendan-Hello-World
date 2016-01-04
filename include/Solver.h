#include "Uniform_Mesh.h"
#include "Boundary_Conditions.h"
#include "Solution.h"


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

};

#endif // SOLVER_H
