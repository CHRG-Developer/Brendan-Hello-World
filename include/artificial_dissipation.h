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
        double * global_JST_switch_x,* global_JST_switch_y;
        void get_global_jst(Solution &soln, Boundary_Conditions &bcs,
                            Uniform_Mesh &Mesh, domain_geometry &domain);
        void set_local_jst();
        void reset_local_jst_switch();
        double local_jst_switch_x, local_jst_switch_y;

    protected:

    private:
};

#endif // ARTIFICIAL_DISSIPATION_H
