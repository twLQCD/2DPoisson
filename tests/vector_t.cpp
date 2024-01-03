#include "../src/lib/vector.cpp"
#include <iostream>
#include <chrono>

int main ()
{
	int n = 100000000;
	//double count = 0.0;
	//int numsamples = 10000;
	//for (int i = 0; i < numsamples; i++){
	Poisson::Vector<double> v(n);
	v.randn();
	auto start = std::chrono::high_resolution_clock::now();
	//Poisson::Vector<float> v(n);
	double norm = v.norm2();
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
	//count += static_cast<double>(duration.count());
	//}
	//std::cout << "Vector construction of size " << n << " took " << count/(static_cast<double>(numsamples)) << " microseconds averaged over " << numsamples << " samples." << std::endl; 
	std::cout << "Norm-2 = " << norm << " and took " << duration.count() << " milliseconds to calculate." << std::endl;
	/*v.zeros();
	v.print();
	v.ones();
	std::cout << "\n" << std::endl;
	v.print();
	v.randn();
	std::cout << "\n" << std::endl;
	v.print();

	//Poisson::Vector<float> *x = v;
	Poisson::Vector x(v);
	std::cout << "\n" << std::endl;
	x.print();

	Poisson::Vector w = x;
        std::cout << "\n" << std::endl;
	w.print();*/
	return 0;
}
