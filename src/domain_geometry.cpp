#include "domain_geometry.h"
#include <math.h>

domain_geometry::domain_geometry()
{
    //ctor
}

domain_geometry::~domain_geometry()
{
    //dtor
}

void domain_geometry::initialise(){

    cs = dx/dt/sqrt(3);


}
