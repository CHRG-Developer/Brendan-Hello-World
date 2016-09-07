#include "Solution.h"
#include "Boundary_Conditions.h"

#ifndef ARTIFICIAL_DISSIPATION_H
#define ARTIFICIAL_DISSIPATION_H


class artificial_dissipation
{
    public:
        artificial_dissipation();
        artificial_dissipation(int nodes);
        virtual ~artificial_dissipation();
        double * global_JST_switch_x,* global_JST_switch_y,global_2nd_order_x,global_2nd_order_y;
        void get_global_jst(Solution &soln, Boundary_Conditions &bcs,
                            Uniform_Mesh &Mesh, domain_geometry &domain);
        void set_local_jst();
        void reset_local_jst_switch();
        double local_jst_switch_x, local_jst_switch_y;
        double kappa_2, kappa_4;
        void get_local_coeffs(Solution &soln, Boundary_Conditions &bcs,
                                            Uniform_Mesh &Mesh, double rho_local);
        void add_local_jst(int j, double rho_local, double rho_neighbour);
        double global_4th_order_x,global_4th_order_y;
        double spectral_radii[4] ; // 1 and 2 = i+ 1;  3 + 4 = j+1
        double martinelli_exponent;

    protected:

    private:
        double jst_x_num,jst_x_den, jst_y_num, jst_y_den;

        double maximum( double m1, double zero, double p1, double p2);
};

#endif // ARTIFICIAL_DISSIPATION_H
