#include <functional>

#include "../src/include/matrices.hpp"
#include "../src/include/vector.hpp"

using namespace Poisson;
using F = float;

int main ()
{

	int xi = 1000; int yi = 1000;
	Grid grid(xi,yi);
	std::cout << "Grid is " << grid.x << " x " << grid.y << std::endl;
	Vector<F> v(grid);
	Vector<F> x(grid);
	std::function<void(const Vector<F>&, Vector<F>&)> func = A<F>;
	v.ones();

	x = v;
	int ns = 1000;
	std::cout << "Norm before matvec is = " << x.norm2() << std::endl;
	for (int i = 0; i < ns; i++) {
	func(v,x);
	}
	std::cout << "Norm after matvec is = " << x.norm2() << std::endl;
	return EXIT_SUCCESS;
}
