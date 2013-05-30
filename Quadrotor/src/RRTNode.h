#ifndef _RRTNODE_
#define _RRTNODE_


#include <vector>
#include <map>
#include <queue>
#include <list>
#include <time.h>
#include <float.h>
#include <fstream>

#include "callisto.h"
#include "matrix.h"
#include "utils.h"



class RRTNode
{
public:
	Matrix<X_DIM> x; //4d
	Matrix<U_DIM> u; //2d
	int bp; //the number of its parent in the vector rrttree.
	double depth;
	double clearance;
	bool isgoal; //check if the node is in the goal region or not
	bool marked;
	int attempts;
	std::vector<int> children; //store the number in the vector rrttree, of its children.

	RRTNode()
	{
		//init position (0,0,0,0)
		x[0] = 0;
		x[1] = 0;
		x[2] = 0;
		x[3] = 0;
		//init control (0,0), acceleration and wheel angel
		u[0] = 0;   
		u[1] = 0;

		bp = -1;
		depth = 0.0;
		isgoal = false;
		marked = false;
		attempts = 0;
		children.clear();

	}
};

#endif