#pragma once

#include <cassert>

namespace Poisson {

struct Grid {

	//dim in x direction
	int x;

	//dim in y direction
	int y;

	//the isotropic spacing

	//default cunstructor
	Grid(){};

	//parameterized constructor
	Grid(int a, int b)
	{
		//assert(a == b);
		x = a;
		y = b;
	//keeping it isotropic
	//on the boundary [0,1] x [0,1]
	}

	void operator=(Grid& grid)
	{ 
		x = grid.x;
		y = grid.y;
	}
	//destructor
	~Grid(){};
};

} //namespace
