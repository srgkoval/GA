#include <iostream>
#include <conio.h>

#include "ga.h"
#include "test.h"

const int n_arg = 2;
const int n_mult = 3;


void main()
{
 	clock_t start, finish;
	start = clock();
    
    GA<n_arg> alg;
   // Vector<double, n_arg> lower, upper;
   // lower = -1.;
   // upper = 1.;

   // alg.run(rastrigin<n_arg>, lower, upper);
   // alg.run(schwefel<n_arg>, Vector<double, n_arg>(-500., -500.), Vector<double, n_arg>(500., 500.));
   // alg.run(rosenbrock<n_arg>, Vector<double, n_arg>(-2., -2.), Vector<double, n_arg>(2., 2.));

    run_ZDT4();



    //int range = 10;
    //for(int i = 0; i < 100; i++)
    //{
    //    double r = range * alg.dist01(alg.rnd_generator);
    //    int ri = floor(r);
    //    std::cout << r << "\t" << ri << "\n";
    //}

    finish = clock();
	std::cout << "\ntimer:\t" << (double (finish - start) / CLOCKS_PER_SEC) << std::endl;
	std::cin.get();

}