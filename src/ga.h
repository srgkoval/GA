#ifndef GA_H
#define GA_H

#include <iostream>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost\random\variate_generator.hpp>

#include <boost/math/constants/constants.hpp>

#include "vector.h"


// forward declarations
template <int N> class GA;
template <int N> class GA_options;  


template <int N> class GA_options
{
public:
    int population_size,
        n_elite;

    int max_generations,
        stall_generations_limit;

    int tournament_size;

    double
        crossover_fraction,
        crossover_BLX_alpha,
        mutation_gaussian_scale,
        mutation_gaussian_shrink;

    typename GA<N>::pFitnessScaling scaling;
    typename GA<N>::pSelection selection;
    typename GA<N>::pCrossover crossover;
    typename GA<N>::pMutation mutation;

    bool verbose;
    
    GA_options()
    {
        population_size = 50;
        n_elite = 2;

        max_generations = 200;
        stall_generations_limit = 50;

        tournament_size = 4;

        crossover_fraction = 0.8;
        crossover_BLX_alpha = 0.5;
        mutation_gaussian_scale = 0.5;
        mutation_gaussian_shrink = 0.75;

        scaling = &GA<N>::scaling_rank;
        selection = &GA<N>::selection_stochastic_uniform;
        crossover = &GA<N>::crossover_BLX;
        mutation = &GA<N>::mutation_adaptive;

        verbose = true;
    }
};


// comparator class to emulate MATLAB's [~,i] = sort(scores) functionality
// allows sorting of int index array according to score values that correspond to these indexes
class index_comparator
{
public:
    double *score;
    index_comparator(double *_score): score(_score) {}
    
    bool operator()(int l, int r)
    {
        return score[l] < score[r];
    }
};


template <int N> class GA
{
public:
    typedef Vector<double, N> Individual;

    typedef void FitnessScaling();
    typedef void(GA<N>::*pFitnessScaling)();

    typedef void Selection(int n);
    typedef void (GA<N>::*pSelection)(int n);

    typedef void Crossover(const Individual &parent1, const Individual &parent2, Individual &child);
    typedef void (GA<N>::*pCrossover)(const Individual &parent1, const Individual &parent2, Individual &child);

    typedef void Mutation(const Individual &parent, Individual &child);
    typedef void (GA::*pMutation)(const Individual &parent, Individual &child);
    
    GA_options<N> options;

    int generation,
        last_improvement_generation;

    Individual *population,
               *children,
               *best_individual;
    
    double *score,
           *fitness,
           *best_score;

    int *score_index,
        *parents;

    Individual lower_boundary,
               upper_boundary;

    FitnessScaling scaling_rank;

    Selection selection_stochastic_uniform,
              selection_tournament;

    Crossover crossover_arithmetic,
              crossover_scattered,
              crossover_BLX;

    Mutation mutation_gaussian,
             mutation_adaptive;
    
    double ma_step_size;                                    // step size for mutation adaptive

    mutable boost::random::mt19937 rnd_generator;
	boost::random::uniform_real_distribution<> dist01;
    boost::random::normal_distribution<> normal01;
    boost::random::uniform_int_distribution<> uniform_int0N;
    boost::random::variate_generator<boost::random::mt19937, boost::random::uniform_int_distribution<>> int_gen;

    GA();
    ~GA();

    Individual random_individual(const typename Individual &lower_boundary, const typename Individual &upper_boundary);
    bool feasible(const Individual &x);
    
    template<typename F> void run(F &f, Individual _lower_boundary, Individual _upper_boundary, Individual *initial_population = NULL, int initial_population_size = 0);

};


template<int N> GA<N>::GA()
    : dist01(0., 1.), generation(0), int_gen(rnd_generator, uniform_int0N)
{
    rnd_generator.seed(static_cast<int>(std::time(NULL)));
    srand(unsigned(time(NULL)));
        
    population = new Individual [options.population_size];
    children =  new Individual [options.population_size];
    best_individual =  new Individual [options.max_generations];

    score = new double [options.population_size];
    fitness = new double [options.population_size];
    best_score = new double [options.max_generations];
    score_index = new int[options.population_size];
    parents = new int [2 * options.population_size];
}


template<int N> GA<N>::~GA()
{
    delete [] population;
    delete [] children;
    delete [] best_individual;

    delete [] score;
    delete [] fitness;
    delete [] best_score;
    delete [] score_index;
    delete [] parents;
}



#include "ga.cpp"

#endif