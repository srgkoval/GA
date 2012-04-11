#include <iostream>
#include <conio.h>

#include "ga.h"
#include "test.h"

const int n_arg = 2;

void main()
{
 	clock_t start, finish;
	start = clock();
    
    GA<n_arg> alg;
    Vector<double, n_arg> lower, upper;

    lower = -1.;
    upper = 1.;
   //alg.run(rastrigin<n_arg>, lower, upper);
   // alg.run(schwefel<n_arg>, Vector<double, n_arg>(-500., -500.), Vector<double, n_arg>(500., 500.));
    alg.run(rosenbrock<n_arg>, Vector<double, n_arg>(-2., -2.), Vector<double, n_arg>(2., 2.));
	
    finish = clock();
	std::cout << "\ntimer:\t" << (double (finish - start) / CLOCKS_PER_SEC) << std::endl;
	std::cin.get();

}