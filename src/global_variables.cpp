#include "global_variables.h"
#include "domain_geometry.h"
#include <boost/algorithm/string/replace.hpp>
#include <boost/filesystem.hpp>
#include <sstream>

using namespace boost::filesystem;
using namespace std;

global_variables::global_variables()
{
    //ctor
}

global_variables::~global_variables()
{
    //dtor
}

void global_variables::initialise(domain_geometry domain){

    tau = 0.5 + max_velocity*domain.X /reynolds_number /domain.dt *pre_conditioned_gamma;

    output_file = create_output_directory();
}


std::string global_variables::create_output_directory(){

    std::string output_file;
    std::string folder;
    std::ostringstream s;
    //output_file = "C:/Users/brendan/Dropbox/PhD/Test Cases/Poiseuille Flow/";

    output_file = output_file_dir;

    if( simulation_name ==  "default"){
        s << "tol " << tolerance << " RE " << reynolds_number
         << " t " << time_marching_step << " gamma " << pre_conditioned_gamma << " tau "
         << tau;
         folder = s.str();

    }else{

        folder = simulation_name;

    }

    //folder.replace(folder.begin(),folder.end(), ".",  "_");
    boost::replace_all(folder, "." , "_");
    output_file = output_file + folder;

    boost::filesystem::path dir(output_file);
    boost::filesystem::create_directories(dir);

    return output_file;

}
