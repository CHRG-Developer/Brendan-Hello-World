#include "quad_bcs.h"

#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H


class Boundary_Conditions
{
    public:
        Boundary_Conditions(int num_x, int num_y);
        virtual ~Boundary_Conditions();
        void assign_boundary_conditions(int num_x, int num_y, quad_bcs);
        bool get_n_bc( int i) {return n_bc[i];};
        bool get_s_bc ( int i) {return s_bc [i];};
        bool get_w_bc( int i) {return w_bc[i];};
        bool get_e_bc( int i) {return e_bc[i];};
        double get_n_rho( int i) {return n_rho[i];};
        double get_n_u( int i) {return n_u[i];};
        double get_n_v ( int i) {return n_v [i];};
        double get_s_rho( int i) {return s_rho[i];};
        double get_s_u( int i) {return s_u[i];};
        double get_s_v ( int i) {return s_v [i];};
        double get_w_rho( int i) {return w_rho[i];};
        double get_w_u( int i) {return w_u[i];};
        double get_w_v ( int i) {return w_v [i];};
        double get_e_rho( int i) {return e_rho[i];};
        double get_e_u( int i) {return e_u[i];};
        double get_e_v( int i) {return e_v[i];};
        int get_n_type( int i) {return n_type[i];};
        int get_e_type( int i) {return e_type[i];};
        int get_w_type( int i) {return w_type[i];};
        int get_s_type( int i) {return s_type[i];};
        int get_periodic_node( int i) {return periodic_node[i];};
    protected:
    private:

        //centroid and cell interface locations

        bool *n_bc, *s_bc , *w_bc, *e_bc;
        double * n_rho, * n_u, * n_v ;
        double * s_rho, * s_u, * s_v ;
        double * w_rho, * w_u, * w_v ;
        double * e_rho, * e_u, * e_v ;
        int *n_type, *e_type, *w_type, *s_type;
        int *periodic_node;
};

#endif // BOUNDARY_CONDITIONS_H
