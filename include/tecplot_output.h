#ifndef TECPLOT_OUTPUT_H
#define TECPLOT_OUTPUT_H
#include <global_variables.h>
#include <Mesh.h>
#include <Solution.h>
#include <Boundary_Conditions.h>


class tecplot_output
{
    public:
        tecplot_output(global_variables &globals, Mesh &Mesh, Solution &Soln,
                             Boundary_Conditions &bcs, int fileType_e, double timestamp);
        virtual ~tecplot_output();

    protected:

    private:
};

#endif // TECPLOT_OUTPUT_H
