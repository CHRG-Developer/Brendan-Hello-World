#include <stdlib.h>
#include "Solution.h"
#include <algorithm>

Solution::Solution(int _total_nodes)
{
    //ctor
    total_nodes = _total_nodes;
     //rho = (double*) malloc (sizeof(double)*(total_nodes));
     rho = new double [total_nodes +1];
        if (rho==NULL) exit (1);
     u = new double [total_nodes+1];
        if (u==NULL) exit (1);
     v = new double [total_nodes +1];
        if (v==NULL) exit (1);
     w = new double [total_nodes +1];
        if (w==NULL) exit (1);
    Initialise();

}

Solution::~Solution()
{
    //dtor
    delete [] rho;
    rho = NULL;
    delete [] u;
    u = NULL;
    delete [] v;
    v= NULL;
    delete [] w;
    w= NULL;

}

void Solution::Initialise() {

    std::fill_n(rho, total_nodes , 1.0);
    std::fill_n(u, total_nodes, 0.0);
    std::fill_n(v, total_nodes , 0.0);
    std::fill_n(w, total_nodes , 0.0);


}
void Solution::update ( double _rho, double _u, double _v, double _w , int i){

    rho[i] =_rho;
    u[i] = _u;
    v[i] = _v;
    w[i] = _w;
}
