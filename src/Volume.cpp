#include "Volume.h"

Volume::Volume()
{
    //ctor

    length = c_length;
	dx = c_dx;
	num_nodes = ceil(length/dx);
	dx = length/num_nodes; // reset dx tp allow for ceiling

	this->create_mesh();
}

Volume::~Volume()
{
    //dtor
}
