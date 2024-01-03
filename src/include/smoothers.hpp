#pragma once

#include <functional>

#include "vector.hpp"
#include "matrices.hpp"

//using Func = std::function;

namespace Poisson {

template<class T>
class Smoother {
        public:
                int iters;
                std::function<void(const Vector<T>&,Vector<T>&)> lu;
                std::function<void(const Vector<T>&,Vector<T>&)> dinv;

                //default constructor
                Smoother<T>(){};

                //parameterized constructor
                Smoother<T>(std::function<void(const Vector<T>&,Vector<T>&)>& lu, std::function<void(const Vector<T>&,Vector<T>&)>& dinv, int iters) :
                        lu(lu),
                        dinv(dinv),
                        iters(iters)
                {}

                void operator()(Vector<T>& in, Vector<T>& out)
                {
                        //temporaries. Maybe a better way?
                        Vector<T> tmp1(out.grid);
                        Vector<T> tmp2(out.grid);
                        Vector<T> x0 = in;
                        for (int i = 0; i < iters; i++)
                        {
                                lu(x0,tmp1);
                                tmp2 = in - tmp1;
                                dinv(tmp2,out);
                                x0 = out;
                        }

                }

                ~Smoother(){}

};

template<class T>
extern void jacobi(Vector<T>& in, Vector<T>& out, int iters, std::function<void(const Vector<T>&,Vector<T>&)>& LU, std::function<void(const Vector<T>&,Vector<T>&)>& Dinv);

} //namespace
