#ifndef TEST_H
#define TEST_H

#include <boost/math/constants/constants.hpp>
#include "ga.h"


template<int N> double rastrigin(const typename GA<N>::Individual &x)
{
    double res = 0.;
    for(int i = 0; i < N; i++)
        res += 10. + x[i] * x[i] - 10. * pow(cos(2 * boost::math::constants::pi<double>() * x[i]), 2);

    return res;
}

template<int N> double schwefel(const typename GA<N>::Individual &x)
{
    double res = 0.;
    for(int i = 0; i < N; i++)
        res += -1. * x[i] * sin(sqrt(fabs(x[i])));
    
    return res;
}

template<int N> double rosenbrock(const typename GA<N>::Individual &x)
{
    double res = 0.;
    for(int i = 0; i < N - 1; i++)
        res += 100. * pow((x[i + 1] - x[i] * x[i]), 2) + pow(1 - x[i], 2);
    
    return res;
}

#endif