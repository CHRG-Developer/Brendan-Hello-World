/*
 * Node.h
 *
 *  Created on: 6 Oct 2015
 *      Author: brend
 */

#ifndef NODE_H_
#define NODE_H_

class Node {
public:
	Node();
	Node(double c_x, double c_y, double c_z);
	virtual ~Node();

	double x,y,z;
	double get_x();
	double get_y();
	double get_z();
};

#endif /* NODE_H_ */
