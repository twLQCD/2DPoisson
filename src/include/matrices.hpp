#pragma once

#include <iostream>
#include <vector>
#include <functional>

#include "vector.hpp"


namespace Poisson {

//calculates the coefficients of the coarse grid matrices
template<class T>
extern std::vector<T> calc_weights(Vector<T>& in, const std::vector<int> nodes, std::function<void(const Vector<T>&,Vector<T>&)>& op, std::function<void(const Vector<T>&,Vector<T>&)>& R, std::function<void(const Vector<T>&,Vector<T>&)>& P);

//the fine grid matvec
template<class T>
extern void op(const Vector<T>& in, Vector<T>& out);

//applies the coarse grid matrix matvec
template<class T>
extern void coarse_op(const Vector<T>& in, Vector<T>& out, std::vector<T>& weights);

//applies (L+U) of A
template<class T>
extern void LU(const Vector<T>& in, Vector<T>& out);

//applies (L+U) of a coarse grid matrix
template<class T>
extern void coarse_LU(const Vector<T>& in, Vector<T>& out, std::vector<T>& weights);

//applies Dinv, the inverse of the diagonal of coarse operator A_c
template<class T>
extern void coarse_Dinv(const Vector<T>& in, Vector<T>& out, std::vector<T>& weights);

//applies Dinv, the inverse of the diagonal of A
template<class T>
extern void Dinv(const Vector<T>& in, Vector<T>& out);

//prolongate
template<class T>
extern void prolong(const Vector<T>& in, Vector<T>& out);

//restrict
template<class T>
extern void restrict(const Vector<T>& in, Vector<T>& out);

}//namespace

