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
template <int N, int N_obj> class GA;
template <int N, int N_obj> class GA_options;  


template <int N, int N_obj> class GA_options
{
public:
    bool multiobjective_problem;

    int population_size,
        n_elite;

    int max_generations,
        stall_generations_limit;

    int tournament_size,
        tournament_size_multiobjective;

    double
        crossover_fraction,
        crossover_BLX_alpha,
        mutation_gaussian_scale,
        mutation_gaussian_shrink;

    typename GA<N, N_obj>::pFitnessScaling scaling;
    typename GA<N, N_obj>::pSelection selection;
    typename GA<N, N_obj>::pCrossover crossover;
    typename GA<N, N_obj>::pMutation mutation;

    bool verbose;
    
    GA_options()
    {
        population_size = 15 * N;
        n_elite = 2;

        max_generations = 100;
        stall_generations_limit = 50;

        tournament_size = 2;
        tournament_size_multiobjective = 2;

        crossover_fraction = 0.8;
        crossover_BLX_alpha = 0.5;
        mutation_gaussian_scale = 0.5;
        mutation_gaussian_shrink = 0.75;

        scaling = &GA<N, N_obj>::scaling_rank;
        selection = &GA<N, N_obj>::selection_tournament;
        crossover = &GA<N, N_obj>::crossover_BLX;
        mutation = &GA<N, N_obj>::mutation_adaptive;

        verbose = true;
    }
};


// comparator class to emulate MATLAB's [~,i] = sort(scores) functionality
// allows sorting of int index array according to score values that correspond to these indexes
template <typename T = double> class index_comparator
{
public:
    T *score;
    index_comparator(T *_score): score(_score) {}
    
    bool operator()(int l, int r)
    {
        return score[l] < score[r];
    }
};


template <int N, int N_obj = 1> class GA
{
public:
    typedef Vector<double, N> Individual;
    typedef Vector<double, N_obj> Objective;

    typedef void FitnessScaling();
    typedef void(GA<N, N_obj>::*pFitnessScaling)();

    typedef void Selection(int n);
    typedef void (GA<N, N_obj>::*pSelection)(int n);

    typedef void Crossover(const Individual &parent1, const Individual &parent2, Individual &child);
    typedef void (GA<N, N_obj>::*pCrossover)(const Individual &parent1, const Individual &parent2, Individual &child);

    typedef void Mutation(const Individual &parent, Individual &child);
    typedef void (GA::*pMutation)(const Individual &parent, Individual &child);
    
    GA_options<N, N_obj> options;

    int generation,
        last_improvement_generation;

    Individual *population,
               *children,
               *best_individual;
    
    Objective *score,
              *best_score;
    
    double *fitness;

    int *score_index,
        *parents;

    // multiobjective data
    Individual *archive;

    int *rank;
    double *distance;
    
    // ---


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
    bool ma_step_changed;                                      // flag for adaptive mutation to change step once in generation

    mutable boost::random::mt19937 rnd_generator;
	boost::random::uniform_real_distribution<> dist01;
    boost::random::normal_distribution<> normal01;
    boost::random::uniform_int_distribution<> uniform_int0N;
    boost::random::variate_generator<boost::random::mt19937, boost::random::uniform_int_distribution<>> int_gen;

    bool first_run;
    void memory_allocate();
    void memory_clear();

    GA();
    ~GA();

    Individual random_individual(const typename Individual &lower_boundary, const typename Individual &upper_boundary);
    bool feasible(const Individual &x);

    void seed_population(Individual *initial_population = NULL, int initial_population_size = 0);
    
    template<typename F> void run(F &f, Individual _lower_boundary, Individual _upper_boundary, Individual *initial_population = NULL, int initial_population_size = 0);
    template<typename F> void run_multiobjective(F &f, Individual _lower_boundary, Individual _upper_boundary,
        Individual *initial_population = NULL, int initial_population_size = 0);

};


template<int N, int N_obj> void GA<N, N_obj>::memory_allocate()
{
    population = new Individual [2 * options.population_size];
    children =  new Individual [options.population_size];
    best_individual =  new Individual [options.max_generations];

    score = new Objective [options.population_size];
    fitness = new double [options.population_size];
    best_score = new Objective [options.max_generations];
    score_index = new int[options.population_size];
    parents = new int [2 * options.population_size];

    archive = new Individual [options.population_size];
}


template<int N, int N_obj> void GA<N, N_obj>::memory_clear()
{
    delete [] population;
    delete [] children;
    delete [] best_individual;

    delete [] score;
    delete [] fitness;
    delete [] best_score;
    delete [] score_index;
    delete [] parents;

    delete [] archive;
}


template<int N, int N_obj> GA<N, N_obj>::GA()
    : dist01(0., 1.), generation(0), int_gen(rnd_generator, uniform_int0N)
{
    rnd_generator.seed(static_cast<int>(std::time(NULL)));
    srand(unsigned(time(NULL)));
        
    memory_allocate();
    first_run = true;
}


template<int N, int N_obj> GA<N, N_obj>::~GA()
{
    memory_clear();
}



#include "ga.cpp"

#endif