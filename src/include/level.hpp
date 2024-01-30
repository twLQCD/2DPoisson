#pragma once

#include "matrices.hpp"
#include "smoothers.hpp"
#include "vector.hpp"

namespace Poisson {


template<class T, class M, class U>
struct Level {

	//the level
	int el;
	//the grid at the level el
	Grid grid;
	//the matrix
	M A;
	//the prolongator
	std::function<void(const Vector<T>&,Vector<T>&)> P;
	//the restrictor
	std::function<void(const Vector<T>&,Vector<T>&)> R;
	//the smoother
	U S;
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
	Level(int el, Grid& grid, M& A, std::function<void(const Vector<T>&,Vector<T>&)>& P, std::function<void(const Vector<T>&,Vector<T>&)>& R, U& S) :
		el(el),
		grid(grid),
		A(A),
		P(P),
		R(R),
		S(S)
	{
		x.create(grid);
		x0.create(grid);
		xf.create(grid);
		b.create(grid);
		r.create(grid);
		e.create(grid);
		Ax.create(grid);
	
	}
	void create(int el_t, Grid& grid_t, M& A_t, std::function<void(const Vector<T>&,Vector<T>&)>& P_t, std::function<void(const Vector<T>&,Vector<T>&)>& R_t, U& S_t)
	{
	
		el = el_t;
		grid = grid_t;
		A = A_t;
		P = P_t;
		R = R_t;
		S = S_t;

                x.create(grid);
                x0.create(grid);
                xf.create(grid);
                b.create(grid);
                r.create(grid);
                e.create(grid);
                Ax.create(grid);

	}
	~Level(){};
}; //Level

} // namespace
