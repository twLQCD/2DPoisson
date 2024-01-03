#pragma once

#include <iostream>

#include "vector.hpp"


namespace Poisson {
//the matvec
template<class T>
extern void A(const Vector<T>& in, Vector<T>& out);

//applies (L+U) of A
template<class T>
extern void LU(const Vector<T>& in, Vector<T>& out);

//applies Dinv, the inverse of the diagonal of A
template<class T>
extern void Dinv(const Vector<T>& in, Vector<T>& out);

//prolongate
template<class T>
extern void P(const Vector<T>& in, Vector<T>& out);

//restrict
template<class T>
extern void R(const Vector<T>& in, Vector<T>& out);

}//namespace

