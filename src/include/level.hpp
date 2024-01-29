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
	//the intial solution at level el
	Vector<T> x0;
	//the final soution at level el
	Vector<T> xf;
	//the right hand side at level el
	Vector<T> b;
	//the residual vector at level el
	Vector<T> r;
	//the error vector at level el
	Vector<T> e;
	//temporary for holding Ax
	Vector<T> Ax;

	Level(){};
	Level(int el, Grid& grid, func<void(const Vector<T>&,Vector<T>&)>& A, func<void(const Vector<T>&,Vector<T>&)>& P, func<void(const Vector<T>&,Vector<T>&)>& R, Smoother<T>& S) :
		el(el),
		grid(grid),
		A(A),
		P(P),
		R(R),
		S(S)
	{
		x(grid);
		x0(grid);
		xf(grid)
		b(grid);
		r(grid);
		e(grid);
		Ax(grid);
	
	}
	~Level(){};
}; //Level

} // namespace
