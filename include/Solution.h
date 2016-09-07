#include "Uniform_Mesh.h"
#include "vector_var.h"
#include <string>
#include <Boundary_Conditions.h>
#ifndef SOLUTION_H
#define SOLUTION_H


class Solution
{
    public:
        Solution(int _total_nodes);
        Solution();
        virtual ~Solution();
        double get_rho( int i) {return rho[i];};
        double get_u( int i) {return u[i];};
        double get_v( int i) {return v[i];};
        double get_w( int i) {return w[i];};
        double get_average_rho (){return average_rho;};
        void set_average_rho(double arg) { average_rho = arg;};
        void set_rho( int i,double arg) {rho[i] =arg;};
        void set_u( int i,double arg) {u[i] =arg;};
        void set_v( int i,double arg) {v[i] =arg;};
        void set_w( int i,double arg) {w[i] =arg;};
        void assign_pressure_gradient( vector_var _gradient, vector_var gradient_origin,
            vector_var origin_magnitude,Uniform_Mesh &Mesh);
        void assign_velocity_gradient( vector_var _gradient, vector_var gradient_origin,
            vector_var origin_magnitude,Uniform_Mesh &Mesh);
        void update ( double rho, double u, double v, double w , int i);
        void output (std::string output_location);
        void clone( Solution &soln_a);
        void post_process(double pre_condition_gamma);
        void add_rho(int i, double arg) { rho[i] = rho[i] + arg;};
        void add_u(int i, double arg) { u[i] = u[i] + arg;}
        void add_v (int i , double arg) {v[i] = v[i] + arg;}
        void add_w (int i, double arg) {w[i] = w[i] + arg;}
        void restriction(Solution &soln, Uniform_Mesh &coarse_mesh,
                           Uniform_Mesh &fine_mesh, Boundary_Conditions &bc);
        void prolongation(Solution &coarse_soln, Solution &temp_soln, Solution &soln,
                            Uniform_Mesh &coarse_mesh, Uniform_Mesh &fine_mesh,
                            Boundary_Conditions &bc, bool fmg);
        void update_bcs(Boundary_Conditions &bcs,Uniform_Mesh &mesh,domain_geometry &domain);

    protected:
    private:
        double *rho, *u, *v, *w;
        int total_nodes;
        double average_rho;
        void Initialise();

};

#endif // SOLUTION_H
