#pragma once

#include "matrices.hpp"
#include "smoothers.hpp"
#include "vector.hpp"

namespace Poisson {

//save some typing
using func = std::function;

template<class T, class M, class U>
struct Level {

	//the level
	int el;
	//the grid at the level el
	Grid grid;
	//the matrix
	M<T> A;
	//the prolongator
	func<void(const Vector<T>&,Vector<T>&)> P;
	//the restrictor
	func<void(const Vector<T>&,Vector<T>&)> R;
	//the smoother
	U<T> S;
	//the solution at the level el
	Vector<T> x;

	Level(){};
	Level(int el, Grid& grid, func& A, func& P, func& R, Smoother& S) :
		el(el),
		grid(grid),
		A(A),
		P(P),
		R(R),
		S(S)
	{x(grid)}
	~Level(){};
}; //Level

} // namespace
