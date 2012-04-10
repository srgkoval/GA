#ifndef GA_CPP
#define GA_CPP

#include "ga.h"
#include "vector.h"

template<int N> template<typename F> void GA<N>::run
    (F &f, Individual lower_boundary, Individual upper_boundary, Individual *initial_population, int initial_population_size)
{
    // seed initial population
    if(initial_population == NULL)
    {
        for(int i = 0; i < options.population_size; ++i)
            population[i] = random_individual(lower_boundary, upper_boundary);
    }
    else
    {
        int initial_size = (initial_population_size == 0 ? options.population_size: initial_population_size);
        
        for(int i = 0; i < initial_size; i++)
        {
            population[i] = initial_population[i];
            if( !(population[i] >= lower_boundary && population[i] <= upper_boundary) )
            {
                std::cout << "GA::run: initial population individual " << i << " is violating boundaries\n";
                wait_and_exit();
            }
        }

        for(int i = initial_size; i < options.population_size; i++)
            population[i] = random_individual(lower_boundary, upper_boundary);
    }
    
    int
        n_crossover_children = floor(0.5 + options.crossover_fraction * (options.population_size - options.n_elite)),
        n_mutation_children = options.population_size - options.n_elite - n_crossover_children,
        n_parents = n_mutation_children + 2 * n_crossover_children;

    for(int gen = 0; gen < options.max_generations; gen++)
    {
        // score the population
        for(int i = 0; i < options.population_size; i++)
            score[i] = f(population[i]);

        // sort score index array accrding to score values in ascending order
        for(int i = 0; i < options.population_size; i++)
            score_index[i] = i;

        index_comparator comp(score);
        std::sort(score_index, score_index + options.population_size, comp);

        best_score[generation] = score[score_index[0]];

        std::cout << best_score[generation] << "\t" << population[score_index[0]] << "\n";

        (this->*options.scaling)();
        (this->*options.selection)(n_parents);

        std::random_shuffle(parents, parents + n_parents);

        for(int i = 0; i < options.n_elite; i++)
            children[i] = population[score_index[i]];
    
        int i_parent = 0;
        for(int i = options.n_elite; i < options.n_elite + n_mutation_children; ++i, ++i_parent)
            (this->*options.mutation) (population[parents[i_parent]], children[i]);

        for(int i = options.n_elite + n_mutation_children; i < options.population_size; ++i, i_parent+=2)
            (this->*options.crossover) (population[parents[i_parent]], population[parents[i_parent + 1]], children[i]);

        for(int i = 0; i < options.population_size; ++i)
            population[i] = children[i];

        generation++;
    }
}


template <int N> void GA<N>::scaling_rank()
{
    for(int i = 0; i < options.population_size; i++)
        fitness[score_index[i]] = 1. / sqrt(i + 1.);  
}


template <int N> void GA<N>::selection_stochastic_uniform(int n)
{
    double *wheel = new double[options.population_size];
    
    double cumulative_fitness = 0;
    for(int i = 0; i < options.population_size; i++)
        cumulative_fitness += fitness[i];
    
    wheel[0] = fitness[0] / cumulative_fitness;
    for(int i = 1; i < options.population_size; ++i)
        wheel[i] = wheel[i - 1] +  fitness[i] / cumulative_fitness;

    double step = 1. / (double)n,
           position = dist01(rnd_generator) * step;
    int j_low = 0;
    
    for(int i = 0; i < n; ++i)
    {
        for(int j = j_low; j < options.population_size; ++j)
        {
            if(position < wheel[j])
            {
                parents[i] = j;
                j_low = j;
                break;
            }
        }
        position += step;
    }
}


template <int N> void GA<N>::crossover_arithmetic
	(const Individual &parent1, const Individual &parent2, Individual &child)
{
    double R = dist01(rnd_generator);
    child = R * parent1 + (1 - R) * parent2;
}


template <int N> void GA<N>::crossover_scattered
	(const Individual &parent1, const Individual &parent2, Individual &child)
{
    for(int i = 0; i < N; ++i)
    {
        if(dist01(rnd_generator) >= 0.5)
        {
            child[i] = parent1[i];
        }
        else
        {
            child[i] = parent2[i];
        }
    }
    
    // may break boundary conditions <===============================================================================!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
}


template <int N> void GA<N>::mutation_gaussian
	(const Individual &parent, Individual &child)
{
    double scale = options.mutation_gaussian_scale * (1. - options.mutation_gaussian_shrink * generation / options.max_generations);
    scale = scale; // * (population_range_upper - population_range_lower);

    for(int i = 0; i < N; ++i)
        child[i] = parent[i] + scale * normal01(rnd_generator);
}


template<int N> typename GA<N>::Individual GA<N>::random_individual
    (const typename Individual &lower_boundary, const typename Individual &upper_boundary)
{
    Individual res;
    for(int i = 0; i < N; ++i)
        res[i] = lower_boundary[i] + dist01(rnd_generator) * (upper_boundary[i] - lower_boundary[i]);
    
    return res;
}



#endif