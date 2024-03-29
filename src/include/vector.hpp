#pragma once

#include <iomanip>
#include <iostream>
#include <random>
#include <cassert>
#include <omp.h>

#include "grid.hpp" 

namespace Poisson {

template<class T>
class Vector {
	public:
		//Grid sets up the 2D grid/mesh. Contains the dimensions
		//grid.x and grid.y
		Grid grid;
		
		//for convenience
		//size = grid.x * grid.y
		int size;

		//the data of the vector
		T* data;

		//default constructor
		Vector<T>(){};

		//copy constructor
		Vector<T>(Vector<T>& v) : grid(v.grid)
       		{ 
			size = v.size;
			data = new T[v.size];
		#pragma omp parallel for
			for (int i = 0; i < v.size; i++) data[i] = (v.data)[i];
		}

		//parameterized constructor
		Vector<T>(Grid& grid) : grid(grid) { size = grid.x * grid.y; data = new T[size]; zeros(); }

		//parameterize constructor
		Vector<T>(Grid& grid, T* data) : grid(grid), data(data) {size = grid.x * grid.y;};

		//to be used when default constructor is called
		void create(const Grid& grid_t)
		{
		 grid = grid_t;
		 size = grid.x*grid.y;
		 data = new T[size];
		 zeros();
		}

		//prints out the elements of the vector
		void print() { for (int i = 0; i < size; i++) std::cout << std::setprecision(std::numeric_limits<T>::max_digits10) << data[i] << std::endl; }

		//set to zero
		void zeros()
		{
		#pragma omp parallel for
			for (int i = 0; i < size; i++) data[i] = static_cast<T>(0);
		}

		//set to ones
		void ones()
		{
		#pragma omp parallel for
			for (int i = 0; i < size; i++) data[i] = static_cast<T>(1);
		}

		void unit(int j)
		{
		this->zeros();
		this->data[j] = 1.0;
		}

		//set to random normal with mean 0
		//and std 1
		void randn()
		{
			std::random_device dev;
			std::mt19937 gen(dev());
			std::normal_distribution<T> d(static_cast<T>(0.0), static_cast<T>(1.0));
		#pragma omp parallel for
			for (int i = 0; i < size; i++) data[i] = d(gen);
		}

		//return the 2-norm of the vector
		T norm2()
		{
			T tmp = static_cast<T>(0);
		#pragma omp parallel for reduction(+ : tmp)
			for (int i = 0; i < size; i++) { tmp += data[i]*data[i]; }
			return tmp;
		}

		//return the sum of a vector
		T sum()
		{
			T tmp = static_cast<T>(0);
		#pragma omp parallel for reduction(+ : tmp)
			for (int i = 0; i < size; i++) {tmp += data[i]; }
			return tmp;
		}

		//overloaded * op -> returns element wise
		//multiplication, not a norm
		Vector<T> operator*(const Vector<T>& v)
		{
			assert(size == v.size);
			T* tmp = new T[v.size];
			Grid grid_t(v.grid.p,v.grid.q);
		#pragma omp parallel for
			for (int i = 0; i < size; i++) {tmp[i] = data[i] * v.data[i];}
			Vector<T> newvec(grid_t, tmp);
			return newvec;
		//#pragma omp parallel for
			//for (int i = 0; i < size; i++) {this->data[i] * v.data[i];}
			//return *this;
		}


		/*friend Vector<T> operator*(Vector<T> lhs, const Vector<T>& rhs)
		{
			lhs *= rhs;
			return lhs;
		}*/

		//overloaded + op
		Vector<T> operator+(const Vector<T>& v) 
                {
                        assert(size == v.size);
                        T* tmp = new T[v.size];
			Grid grid_t(v.grid.p,v.grid.q);
                #pragma omp parallel for
                        for (int i = 0; i < size; i++) {tmp[i] = data[i] + v.data[i];}
                        Vector<T> newvec(grid_t, tmp);
                        return newvec;
		//#pragma omp parallel for
			//for (int i = 0; i < size; i++) {this->data[i] + v.data[i];}
			//return *this;
                }

		/*friend Vector<T> operator+(Vector<T> lhs, const Vector<T>& rhs)
		{
			lhs += rhs;
			return lhs;
		}*/

		//overloaded - op
		Vector<T> operator-(const Vector<T>& v)
                {
                        assert(size == v.size);
                        T* tmp = new T[v.size];
			Grid grid_t(v.grid.p,v.grid.q);
                #pragma omp parallel for
                        for (int i = 0; i < size; i++) {tmp[i] = data[i] - v.data[i];}
                        Vector<T> newvec(grid_t, tmp);
                        return newvec;
		//#pragma omp parallel for
		//	for (int i = 0; i < size; i++) {this->data[i] - v.data[i];}
		//	return *this;
                }

		/*friend Vector<T>& operator-(Vector<T> lhs, const Vector<T>& rhs)
		{
			lhs += rhs;
			return lhs;
		}*/

		//assignment operator
		void operator=(const Vector<T>& v)
		{
			//grid = v.grid;
			//size = grid.x * grid.y;
			//size = v.size;
			//grid = v.grid;
			assert(this->size == v.size);
		#pragma omp parallel for
			for (int i = 0; i < size; i++) {data[i] = v.data[i];}
			//return *this;
			/*if (this == &v) {return *this;}
			assert(this -> size == v.size);
			std::copy(v.data, v.data + v.size, this -> data);
			return *this;*/
		}

		//return the element of the vector given a linear index
		T operator()(int i) const
		{
			return data[i];
		}

		//return the element of the vector given a two dimensional
		//index
		/*T operator()(int i, int j)
		{
			return data[i + j*grid.y];
		} */

		//destructor
		~Vector<T>(){ delete[] data; }


}; //class
} //namespace
