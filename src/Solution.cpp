#include "Solution.h"
#include <algorithm>

Solution::Solution(int _total_nodes)
{
    //ctor
    total_nodes = _total_nodes;
     rho = (double*) malloc (total_nodes+1);
        if (rho==NULL) exit (1);
     u = (double*) malloc (total_nodes+1);
        if (u==NULL) exit (1);
     v = (double*) malloc (total_nodes+1);
        if (v==NULL) exit (1);
     w = (double*) malloc (total_nodes+1);
        if (w==NULL) exit (1);
    Initialise();

}

Solution::~Solution()
{
    //dtor
}

void Solution::Initialise() {

    std::fill_n(rho, total_nodes+1 , 1);
    std::fill_n(u, total_nodes+1 , 0);
    std::fill_n(v, total_nodes+1 , 0);
    std::fill_n(w, total_nodes+1 , 0);


}
void Solution::update ( double _rho, double _u, double _v, double _w , int i){

    rho[i] =_rho;
    u[i] = _u;
    v[i] = _v;
    w[i] = _w;
}
