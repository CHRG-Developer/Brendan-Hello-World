#include "Uniform_Mesh.h"

#ifndef SOLUTION_H
#define SOLUTION_H


class Solution
{
    public:
        Solution(int _total_nodes);
        virtual ~Solution();
        double get_rho( int i) {return rho[i];};
        double get_u( int i) {return u[i];};
        double get_v( int i) {return v[i];};
        double get_w( int i) {return w[i];};

        void update ( double rho, double u, double v, double w , int i);
    protected:
    private:
        double * rho, * u, * v, * w;
        int total_nodes;
        void Initialise();
};

#endif // SOLUTION_H
