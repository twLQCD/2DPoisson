#include <functional>

#include "../src/include/matrices.hpp"
#include "../src/include/vector.hpp"
#include "../src/include/smoothers.hpp"

using namespace Poisson;
using F = float;

int main ()
{

	int xi = 32; int yi = 32;
	Grid grid(xi,yi);
	std::cout << "Grid is " << grid.x << " x " << grid.y << std::endl;
	Vector<F> v(grid);
	Vector<F> x(grid);
	std::function<void(const Vector<F>&, Vector<F>&)> func = A<F>;
	v.ones();
	x.ones();
	std::function<void(const Vector<F>&, Vector<F>&)> lu = LU<F>;
	std::function<void(const Vector<F>&, Vector<F>&)> dinv = Dinv<F>;

//All checks below passed for over loaded ops
/*
	//test the overloaded - operator
	Vector<F> w = x - v;
	w.print(); //should be zeros
	std::cout << "\n" << std::endl;
	x.zeros();

	//test the overloaded + op
	x.ones();
	w = x + v;
	w.print(); //should be twos;
	std::cout << "\n" << std::endl;

	//test the overloaded * operator
	x.zeros();
	v.randn();
	w = x*v;
	w.print(); //should be zeros
*/	

	int iters = 15000;
	//v.randn();
//below implementation works
        //Vector<F> x0(grid);
        //Vector<F> tmp1(grid);
        //Vector<F> tmp2(grid);
        //Vector<F> r(grid);
	/*for (int i = 0; i < iters; i++) 
	{
	lu(x0,tmp1);
        tmp2 = (v - tmp1);
	dinv(tmp2,x); //should be (1/4) of v;
	x0 = x;
	}
	
	func(x0,tmp1);
	r = tmp1 - v;
	std::cout << r.norm2() << std::endl; */

	//std::function<void(Vector<F>&, Vector<F>&, int, std::function<void(const Vector<F>&,Vector<F>&)>, std::function<void(const Vector<F>&,Vector<F>&)>)> jac = jacobi<F>;
	x.zeros();
	std::cout << "Norm of rhs before jacobi is = " << v.norm2() << std::endl;
	Smoother<F> jacobi(lu, dinv, iters);
	jacobi(v,x);
	//jacobi<F>(v,x,iters,lu,dinv);
	Vector<F> r(grid);
	Vector<F> tmp(grid);
	std::cout << "Norm of rhs after jacobi is = " << v.norm2() << std::endl;
	std::cout << "Norm of solution is = " << x.norm2() << std::endl;
	func(x,tmp);
	r = v - tmp;
	std::cout << "Norm of residual vector is = " << r.norm2() << std::endl;
	return EXIT_SUCCESS;
}
