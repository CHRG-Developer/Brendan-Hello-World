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
    //tau = 0.5 + viscosity / dt/ gamma
    // viscosity = MA/root(3) /Re
     std::ostringstream s;
    tau = 0.5 + mach_number * sqrt(3)/reynolds_number * domain.Y/domain.dt*pre_conditioned_gamma;
    knudsen_number = mach_number / reynolds_number;
    s << "RE_" << reynolds_number << " N_CELLS_" << domain.Y <<
                    " MA_" << mach_number << " dt_" << domain.dt
                    << " DT_" << time_marching_step;
    simulation_name = s.str();
    boost::replace_all(simulation_name,".","_");
    output_file = create_output_directory();

}

void global_variables::update_coarse_tau(){

    tau = tau/2.0;
}

void global_variables::update_fine_tau(){

    tau = tau*2.0;
}

void global_variables::update_tau( domain_geometry domain){
     tau = 0.5 + mach_number * sqrt(3)/reynolds_number * ceil(domain.Y /domain.dt) *pre_conditioned_gamma;

}
void global_variables::magnify_time_step( ){
     time_marching_step = time_marching_step * 2;

}

void global_variables::reduce_time_step( ){
     time_marching_step = time_marching_step / 2;

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
    boost::system::error_code ec;
    boost::replace_all(folder, "." , "_");
    output_file = output_file + folder;

    boost::filesystem::path dir(output_file);
    boost::filesystem::create_directories(dir);
    output_file = output_file;
    return output_file;

}
