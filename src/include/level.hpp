#pragma once

#include "matrices.hpp"
#include "smoothers.hpp"
#include "vector.hpp"

namespace Poisson {

//save some typing
using std::function func;

template<class T>
struct Level {

	//the level
	int el;
	//the grid at the level el
	Grid grid;
	//the matrix
	func<void(const Vector<T>&,Vector<T>&)> A;
	//the prolongator
	func<void(const Vector<T>&,Vector<T>&)> P;
	//the restrictor
	func<void(const Vector<T>&,Vector<T>&)> R;
	//the smoother
	Smoother<T> S;
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
