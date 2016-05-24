#ifndef GLOBAL_VARIABLES_H
#define GLOBAL_VARIABLES_H
#include <string>
#include "domain_geometry.h"

class global_variables
{
    public:
        global_variables();
        virtual ~global_variables();
        void initialise(domain_geometry domain);
        std::string create_output_directory();

        //BC constants
        const int periodic  =3;
        const int dirichlet = 1;
        const int neumann = 2;

        //small number constant
        const double small_number = 0.00000000001;


        //tolerance in final solution
        double tolerance = 0.000001;

        double pre_conditioned_gamma = 1;
        double simulation_length = 2000; // cut off time in seconds
        double time_marching_step;
        double reynolds_number;
        double max_velocity;
        double tau;


        std::string output_file_dir;
        std::string simulation_name;
        std::string output_file;

    protected:
    private:
};

#endif // GLOBAL_VARIABLES_H
