#include "Uniform_Mesh.h"
#include <math.h>
#include <stdio.h>      /* printf, scanf, NULL */
#include <stdlib.h>
#include <algorithm>
Uniform_Mesh::Uniform_Mesh(double c_X, double c_Y, double c_dx, double c_dy)
{
    //ctor
    X = c_X;
    Y = c_Y;
    dx = c_dx;
    dy = c_dy;

	num_x_nodes = ceil(X/dx);
	num_y_nodes = ceil(Y/dy);
	total_nodes = num_x_nodes * num_y_nodes;

	dx = X/num_x_nodes; // reset dx/dy t0 allow for ceiling
    dy = Y/num_y_nodes;

    // allocate memory to arrays

    centroid_x = (double*) malloc (total_nodes+1);
        if (centroid_x==NULL) exit (1);
    centroid_y = (double*) malloc (total_nodes+1);
        if (centroid_y==NULL) exit (1);
    centroid_z = (double*) malloc (total_nodes+1);
        if (centroid_z==NULL) exit (1);
    north_x = (double*) malloc (total_nodes+1);
        if (north_x==NULL) exit (1);
    north_y = (double*) malloc (total_nodes+1);
        if (north_y==NULL) exit (1);
    north_z = (double*) malloc (total_nodes+1);
        if (north_z==NULL) exit (1);
    east_x = (double*) malloc (total_nodes+1);
        if (east_x==NULL) exit (1);
    east_y = (double*) malloc (total_nodes+1);
        if (east_y==NULL) exit (1);
    east_z = (double*) malloc (total_nodes+1);
        if (east_z==NULL) exit (1);
    west_x = (double*) malloc (total_nodes+1);
        if (west_x==NULL) exit (1);
    west_y = (double*) malloc (total_nodes+1);
        if (west_y==NULL) exit (1);
    west_z = (double*) malloc (total_nodes+1);
        if (west_z==NULL) exit (1);
    south_x = (double*) malloc (total_nodes+1);
        if (south_x==NULL) exit (1);
    south_y = (double*) malloc (total_nodes+1);
        if (south_y==NULL) exit (1);
    south_z = (double*) malloc (total_nodes+1);
        if (south_z==NULL) exit (1);

    n_area = (double*) malloc (total_nodes+1);
        if (n_area==NULL) exit (1);
    s_area = (double*) malloc (total_nodes+1);
        if (s_area==NULL) exit (1);
    w_area = (double*) malloc (total_nodes+1);
        if (w_area==NULL) exit (1);
    e_area = (double*) malloc (total_nodes+1);
        if (e_area==NULL) exit (1);

    n_i = (double*) malloc (total_nodes+1);
        if (n_i==NULL) exit (1);
    n_j = (double*) malloc (total_nodes+1);
        if (n_j==NULL) exit (1);
    n_k = (double*) malloc (total_nodes+1);
        if (n_k==NULL) exit (1);
     e_i = (double*) malloc (total_nodes+1);
        if (e_i==NULL) exit (1);
    e_j = (double*) malloc (total_nodes+1);
        if (e_j==NULL) exit (1);
    e_k = (double*) malloc (total_nodes+1);
        if (e_k==NULL) exit (1);
     w_i = (double*) malloc (total_nodes+1);
        if (w_i==NULL) exit (1);
    w_j = (double*) malloc (total_nodes+1);
        if (w_j==NULL) exit (1);
    w_k = (double*) malloc (total_nodes+1);
        if (w_k==NULL) exit (1);
     s_i = (double*) malloc (total_nodes+1);
        if (s_i==NULL) exit (1);
    s_j = (double*) malloc (total_nodes+1);
        if (s_j==NULL) exit (1);
    s_k = (double*) malloc (total_nodes+1);
        if (s_k==NULL) exit (1);


     n_node = (int*) malloc (total_nodes+1);
        if (n_node==NULL) exit (1);
    s_node = (int*) malloc (total_nodes+1);
        if (s_node==NULL) exit (1);
    w_node = (int*) malloc (total_nodes+1);
        if (w_node==NULL) exit (1);
    e_node = (int*) malloc (total_nodes+1);
        if (e_node==NULL) exit (1);

    delta_t = (double*) malloc (total_nodes+1);
        if (delta_t==NULL) exit (1);

	this->create_mesh();
}

Uniform_Mesh::~Uniform_Mesh()
{
    //dtor
}

void Uniform_Mesh::create_mesh(){
    int counter =0;
    for( int i=0; i < num_x_nodes; i++){
        for( int j=0; j < num_y_nodes; j++){
            centroid_x[counter] = dx/2 + i*dx;
            centroid_y[counter] = dy/2 + j*dy;
            centroid_z[counter] = 0; //temporary
            north_x[counter] = dx/2 + i*dx;
            north_y[counter] = dy + j*dy;
            north_z[counter] = 0; //temporary
            south_x[counter] = dx/2 + i*dx;
            south_y[counter] = j*dy;
            south_z[counter] = 0; //temporary
            west_x[counter] = i*dx;
            west_y[counter] = dy/2 + j*dy;
            west_z[counter] = 0; //temporary
            east_x[counter] = dx + i*dx;
            east_y[counter] = dy/2 + j*dy;
            east_z[counter] = 0; //temporary

            n_area[counter] = dx;
            s_area[counter] = dx;
            e_area[counter] = dy;
            w_area[counter] = dy;

            n_i[counter] = 0.0;
            n_j[counter] = 1.0;
            n_k[counter] = 0; //temporary

            e_i[counter] = 1.0;
            e_j[counter] = 0.0;
            e_k[counter] = 0; //temporary


            s_i[counter] = 0.0;
            s_j[counter] = -1.0;
            s_k[counter] = 0; //temporary


            w_i[counter] = -1.0;
            w_j[counter] = 0.0;
            w_k[counter] = 0; //temporary

            delta_t[counter] = 0.5 * std::min(dy,dx);


            counter = counter +1;
        }
	}


}
double  Uniform_Mesh::get_node_x(int node_num){
    double result;
    result = centroid_x[node_num];

    return result ;
}
