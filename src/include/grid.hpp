#pragma once

#include <cassert>
#include <cmath>

namespace Poisson {

struct Grid {

	//power to make dims in the x direction
	int p;

	//dim in x direction
	int x;

	//power to make dims in the y direction 
	int q;

	//dim in y direction
	int y;

	//the isotropic spacing

	//default cunstructor
	Grid(){};

	//parameterized constructor
	Grid(int a, int b) : p(a), q(b)
	{
		//assert(a == b);
		x = static_cast<int>(pow(2.0, a) - 1);
		y = static_cast<int>(pow(2.0, b) - 1);
	//keeping it isotropic
	//on the boundary [0,1] x [0,1]
	}

	void operator=(Grid& grid)
	{ 	p = grid.p;
		q = grid.q;
		x = grid.x;
		y = grid.y;
	}
	//destructor
	~Grid(){};
};

} //namespace
