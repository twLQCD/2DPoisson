//#pragma once

//#include <omp.h>

#include "../include/matrices.hpp"


namespace Poisson {
//the matvec: applies A*x
template<class T>
void A(const Vector<T>& in, Vector<T>& out)
{
	assert(in.size == out.size);
	int i; int j;
	T four = static_cast<T>(4.0);
	T mone = static_cast<T>(-1.0);

	//the origin
	out.data[0] = mone*(in(1) + in(in.grid.y) - four*in(0));
	//corner opposite the origin
	out.data[out.grid.x*out.grid.y-1] = mone*( in(in.grid.x * in.grid.y - in.grid.x - 1) + in( in.grid.x * in.grid.y - 2 ) - four*in(in.grid.x * in.grid.y - 1));
	//corner adjacent to origin
	out.data[out.grid.x-1] = mone*(in(in.grid.x-2) + in(in.grid.x + in.grid.y - 1) - four*in(in.grid.x-1));
	//corner above the origin
	out.data[out.grid.x*out.grid.y - out.grid.x] = mone*(in( (in.grid.x-1)*in.grid.y ) + in(in.grid.x*in.grid.y - in.grid.x + 1) -four*in(in.grid.x*in.grid.y - in.grid.x));
	//the bottom/top edge
#pragma omp parallel for
	for (i = 1; i < in.grid.x-1; i++)
	{
	out.data[i] = mone*(in(i-1) + in(i+1) + in(i + in.grid.x) - four*in(i));
	out.data[i + (out.grid.y-1)*out.grid.x] = mone*(in(i-1 + (in.grid.y-1)*in.grid.x) + in(i+1 + (in.grid.y-1)*in.grid.x) + in(i + (in.grid.y-2)*in.grid.x) - four*in(i + (in.grid.y-1)*in.grid.x));
	}

	//the left/right edge
#pragma omp parallel for
	for (i = 1; i < in.grid.y-1; i++)
	{
	out.data[(i+1)*out.grid.x-1] = mone*(in(i*in.grid.x - 1) + in((i+1)*in.grid.x-2) + in((i+2)*in.grid.x-1) - four*in((i+1)*in.grid.x-1));
	out.data[(i-1)*out.grid.x + out.grid.x] = mone*(in(i*in.grid.x + in.grid.x) + in( (i-2)*in.grid.x + in.grid.x) + in( (i-1)*in.grid.x + in.grid.x - 1) - four*in((i-1)*in.grid.x + in.grid.x));
	}

	//now the interior nodes
#pragma omp parallel for private(i) shared(j) collapse(2)
	for (j = 1; j < in.grid.y - 1; j++){
		for (i = 1; i < in.grid.x - 1; i++) {
			out.data[i + j*out.grid.x] = (mone)*(in(i + (j-1)*in.grid.x) + in( (i-1) + j*in.grid.x) - four*in(i + j*in.grid.x) + in( (i+1) + j*in.grid.x ) + in(i + (j+1)*in.grid.x));
		}

	}
}
template void A<float>(const Vector<float>& in, Vector<float>& out);
template void A<double>(const Vector<double>& in, Vector<double>& out);

//applies the lower plus upper part of the matrix A
template<class T>
void LU(const Vector<T>& in, Vector<T>& out)
{
        assert(in.size == out.size);
        int i; int j;
        T mone = static_cast<T>(-1.0);
	out.data[0] = mone*(in(1) + in(in.grid.y));
	out.data[out.grid.x*out.grid.y-1] = mone*( in(in.grid.x * in.grid.y - in.grid.x - 1) + in( in.grid.x * in.grid.y - 2 ));
	out.data[out.grid.x-1] = mone*(in(in.grid.x-2) + in(in.grid.x + in.grid.y - 1));
	out.data[out.grid.x*out.grid.y - out.grid.x] = mone*(in( (in.grid.x-1)*in.grid.y ) + in(in.grid.x*in.grid.y - in.grid.x + 1));
#pragma omp parallel for
	for (i = 1; i < in.grid.x-1; i++)
	{
	out.data[i] = mone*(in(i-1) + in(i+1) + in(i + in.grid.x));
	out.data[i + (out.grid.y-1)*out.grid.x] = mone*(in(i-1 + (in.grid.y-1)*in.grid.x) + in(i+1 + (in.grid.y-1)*in.grid.x) + in(i + (in.grid.y-2)*in.grid.x));
	}
#pragma omp parallel for
	for (i = 1; i < in.grid.y-1; i++)
	{
	out.data[(i+1)*out.grid.x-1] = mone*(in(i*in.grid.x - 1) + in((i+1)*in.grid.x-2) + in((i+2)*in.grid.x-1));
	out.data[(i-1)*out.grid.x + out.grid.x] = mone*(in(i*in.grid.x + in.grid.x) + in( (i-2)*in.grid.x + in.grid.x) + in( (i-1)*in.grid.x + in.grid.x - 1));
	}	
#pragma omp parallel for private(i) shared(j) collapse(2)
        for (j = 1; j < in.grid.y - 1; j++){
                for (i = 1; i < in.grid.x - 1; i++) {
                        out.data[i + j*out.grid.x] = (mone)*(in(i + (j-1)*in.grid.x) + in( (i-1) + j*in.grid.x) + in( (i+1) + j*in.grid.x ) + in(i + (j+1)*in.grid.x));
                }

        }
}
template void LU<float>(const Vector<float>& in, Vector<float>& out);
template void LU<double>(const Vector<double>& in, Vector<double>& out);

//applies the inv of the diagonal of the matrix
template<class T>
void Dinv(const Vector<T>& in, Vector<T>& out)
{
	assert(in.size == out.size);
	int i; //int j;
	T quart = static_cast<T>(0.25);
#pragma omp parallel for
	for (i = 0; i < out.size; i++) { out.data[i] = quart*in(i); }
}
template void Dinv<float>(const Vector<float>& in, Vector<float>& out);
template void Dinv<double>(const Vector<double>& in, Vector<double>& out);

//prolongate
template<class T>
void P(const Vector<T>& in, Vector<T>& out)
{
	std::cout << "Not yet Implemented" << std::endl;
	return;
}
template void P<float>(const Vector<float>& in, Vector<float>& out);
template void P<double>(const Vector<double>& in, Vector<double>& out);

//restrict
template<class T>
void R(const Vector<T>& in, Vector<T>& out)
{
        std::cout << "Not yet Implemented" << std::endl;
	return;
}
template void R<float>(const Vector<float>& in, Vector<float>& out);
template void R<double>(const Vector<double>& in, Vector<double>& out);
}//namespace

