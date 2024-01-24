#pragma once

#include <functional>
#include <cmath>

#include "vector.hpp"
#include "matrices.hpp"

//using Func = std::function;

namespace Poisson {

template<class T>
class Smoother {
        public:
                int iters;
		T target;
                std::function<void(const Vector<T>&,Vector<T>&)> lu;
                std::function<void(const Vector<T>&,Vector<T>&)> dinv;
		std::function<void(const Vector<T>&,Vector<T>&)> A;

                //default constructor
                Smoother<T>(){};

                //parameterized constructor
                Smoother<T>(std::function<void(const Vector<T>&,Vector<T>&)>& lu, std::function<void(const Vector<T>&,Vector<T>&)>& dinv, std::function<void(const Vector<T>&,Vector<T>&)>& A, int iters, T target) :
                        lu(lu),
                        dinv(dinv),
			A(A),
                        iters(iters),
			target(target)
                {}

                void operator()(const Vector<T>& in, Vector<T>& out)
                {
                        //temporaries. Maybe a better way?
                        Vector<T> tmp1(out.grid);
                        Vector<T> tmp2(out.grid);
                        Vector<T> x0(out.grid);
			Vector<T> bnew(out.grid);
			Vector<T> bin(out.grid);
			Vector<T> xin(out.grid);
			Vector<T> r(out.grid);
			xin = out;
			bin = in;
			//compute the initial residual. If out is zeros, it will be the rhs
			A(xin,tmp1);
			bnew = bin - tmp1;
			T rn = std::sqrt(bin.norm2());
			//std::cout << "Initial residual norm is : " << bnew.norm2() << std::endl;
			//std::cout << "Norm of rhs in Jacobi is : " << xin.norm2() << std::endl;
			int i = 0;
			while (i <= iters && rn >= target)
                        {
                                lu(x0,tmp1);
                                tmp2 = bnew - tmp1;
                                dinv(tmp2,tmp1);
                                x0 = tmp1;

				xin = x0 + out;
				A(xin,tmp2);
				r = tmp2 - bin;
				rn = std::sqrt(r.norm2());
				//std::cout << "Iteration: " << i << ", ||r||_2 = " << rn << std::endl;
				i++;
                        }
			out = xin;

                }

                ~Smoother(){}

};

template<class T>
extern void jacobi(Vector<T>& in, Vector<T>& out, int iters, std::function<void(const Vector<T>&,Vector<T>&)>& LU, std::function<void(const Vector<T>&,Vector<T>&)>& Dinv);

} //namespace
