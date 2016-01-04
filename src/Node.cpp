/*
 * Node.cpp
 *
 *  Created on: 6 Oct 2015
 *      Author: brend
 */

#include "Node.h"
#include <string>

Node::Node(){

}

Node::Node(double c_x, double c_y, double c_z) {
	// TODO Auto-generated constructor stub
	x = c_x;
	y= c_y;
	z = c_z;
}
double Node::get_x(){
	return x;
}
double Node::get_y(){
	return y ;
}
double Node::get_z(){
	return z;
}

Node::~Node() {
	// TODO Auto-generated destructor stub
}

