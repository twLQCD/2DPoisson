#include <functional>
#include <vector>
#include <memory>

#include "../src/include/matrices.hpp"
#include "../src/include/vector.hpp"
#include "../src/include/smoothers.hpp"
#include "../src/include/level.hpp"
#include "../src/include/vcycle.hpp"

using namespace Poisson;
using T = double;

int main ()
{
	//p and q are related to the dims x and y of the grid by:
	//x = 2^p + 1, y = 2^q + 1
	//using 4 for each gives a 15 x 15 grid
	
	int p = 6; int q = 6;
	Vcycle<T> vcycle = setup(p,q,3,1e-14,10000,1e-12);

	vcycle.fine_level->b.ones();
	int iters = 10000;
	T target = std::sqrt(vcycle.fine_level->b.norm2())*1e-06;
	vcycle(iters,target);
	return EXIT_SUCCESS;
}
