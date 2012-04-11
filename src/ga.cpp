#ifndef GA_CPP
#define GA_CPP

#include "ga.h"
#include "vector.h"
#include "support.h"

template<int N> template<typename F> void GA<N>::run
    (F &f, Individual _lower_boundary, Individual _upper_boundary, Individual *initial_population, int initial_population_size)
{
    lower_boundary = _lower_boundary;
    upper_boundary = _upper_boundary;

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
            if( !feasible(population[i]) )
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
    
    generation = 0;
    last_improvement_generation = 0;

    for(int gen = 0; gen < options.max_generations; gen++)
    {
        // score the population
        for(int i = 0; i < options.population_size; i++)
            score[i] = f(population[i]);

        // sort score index array according to score values in ascending order
        for(int i = 0; i < options.population_size; i++)
            score_index[i] = i;

        index_comparator comp(score);
        std::sort(score_index, score_index + options.population_size, comp);

        best_score[generation] = score[score_index[0]];
        best_individual[generation] = population[score_index[0]];

        if(options.verbose)
        {
            std::cout << generation << "\t" << best_score[generation] << "\t" << population[score_index[0]] << "\n";
        }

        if(generation == options.max_generations - 1)
            break;

        if(generation > 0)
        {
            if(best_score[generation] < best_score[generation - 1])
                last_improvement_generation = generation;
        }

        // break if no improvement for last 'options.stall_generations_limit' generations
        if(generation - last_improvement_generation > options.stall_generations_limit)
            break;

        (this->*options.scaling)();
        
        // normalize fitness so that total fitness sum is 1.
        double cumulative_fitness = 0.;
        for(int i = 0; i < options.population_size; i++)
            cumulative_fitness += fitness[i];
        for(int i = 0; i < options.population_size; i++)
            fitness[i] /= cumulative_fitness;
        
        // get parents' indexes
        (this->*options.selection)(n_parents);

        std::random_shuffle(parents, parents + n_parents);

        // transfer elite children
        for(int i = 0; i < options.n_elite; i++)
            children[i] = population[score_index[i]];
        
        // perform genetic operators to get the rest of the children
        ma_step_changed = false;
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


// scaling ========================================================================================

template <int N> void GA<N>::scaling_rank()
{
    for(int i = 0; i < options.population_size; i++)
        fitness[score_index[i]] = 1. / sqrt(i + 1.);  
}


// selection ======================================================================================

template <int N> void GA<N>::selection_stochastic_uniform(int n)
{
    double *wheel = new double[options.population_size];
    
    wheel[0] = fitness[0];
    for(int i = 1; i < options.population_size; ++i)
        wheel[i] = wheel[i - 1] +  fitness[i];

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


template <int N> void GA<N>::selection_tournament(int n)
{
    int *indexes = new int [options.population_size];
    for(int i = 0; i < options.population_size; ++i)
        indexes[i] = i;

    for(int i = 0; i < n; i++)
    {
        std::random_shuffle(indexes, indexes + options.population_size);
        
        parents[i] = indexes[0];
        for(int j = 1; j < options.tournament_size; j++)
        {
            if(fitness[indexes[j]] > fitness[parents[i]])
                parents[i] = indexes[j];
        }
    }

    delete [] indexes;
}


// crossover ======================================================================================

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
    
    if( !feasible(child) )
        crossover_arithmetic(parent1, parent2, child);
}


template <int N> void GA<N>::crossover_BLX
	(const Individual &parent1, const Individual &parent2, Individual &child)
{
    for(int i = 0; i < N; i++)
    {
        double R = dist01(rnd_generator);

        double I = parent2[i] - parent1[i],
               low = parent1[i] - I * options.crossover_BLX_alpha,
               high = parent2[i] + I * options.crossover_BLX_alpha;

        child[i] = low + R * (high - low);

        // modify child to satisfy constraints
        //child[i] = std::max(child[i], lower_boundary[i]);
        //child[i] = std::min(child[i], upper_boundary[i]);
        if(child[i] < lower_boundary[i] || child[i] > upper_boundary[i])
            child[i] = (parent1[i] + parent2[i]) / 2.;
    }
}


// mutation =======================================================================================

template <int N> void GA<N>::mutation_gaussian
	(const Individual &parent, Individual &child)
{
    // this mutation breaks constraints
    double scale_factor = options.mutation_gaussian_scale * (1. - options.mutation_gaussian_shrink * generation / options.max_generations);
    Individual scale = scale_factor * (upper_boundary - lower_boundary);

    for(int i = 0; i < N; ++i)
        child[i] = parent[i] + scale[i] * normal01(rnd_generator);
}


template <int N> void GA<N>::mutation_adaptive
	(const Individual &parent, Individual &child)
{
    // set step size
    if(generation <= 1)
    {
        ma_step_size = 1.;
    }
    else
    {
        if(!ma_step_changed)
        {
            if(best_score[generation] < best_score[generation - 1])
                ma_step_size = std::min(1., ma_step_size * 4.);
            else
                ma_step_size = std::max(1.e-8, ma_step_size / 4.);
            
            ma_step_changed = true;
        }
    }

    // set logarithmic scale
    Individual scale;
    for(int i = 0; i < N; i++)
    {
        double exponent = 0.5 * ( log(fabs(lower_boundary[i])) + log(fabs(upper_boundary[i])) ) / log(2.);
        scale[i] = pow(2., exponent);
    }

    Individual raw_basis[N],
               basis[N],
               tangent_cone[N],
               dir[2 * N];
    double dir_sign[4 * N];

    double tol = 1.e-8;                                                             // pure magic

    int n_tangent = 0,
        n_basis = N,
        index_vector[4 * N],
        order_vector[4 * N];

    // calculate mutation direction set
    // tangent components
    for(int i = 0; i < N; i++)
    {
        if( fabs(parent[i] - lower_boundary[i]) < tol || fabs(parent[i] - upper_boundary[i]) < tol )
        {
            tangent_cone[n_tangent] = 0.;
            tangent_cone[n_tangent][i] = 1.;
            n_tangent++;
        }
    }
    
    // raw basis vectors
    double poll_param = 1. / sqrt(ma_step_size);

    for(int i = 0; i < n_basis; i++)
    {
        raw_basis[i] = 0.;
        raw_basis[i][i] = poll_param * (dist01(rnd_generator) >= 0.5 ? 1. : -1.);
        for(int j = i + 1; j < n_basis; j++)
        {
            raw_basis[i][j] = round( (poll_param + 1.) * dist01(rnd_generator) - 0.5);
        }
    }

    // basis as random permutation of raw basis
    for(int j = 0; j < n_basis; j++)
        order_vector[j] = j;
    std::random_shuffle(order_vector, order_vector + n_basis);

    for(int i = 0; i < n_basis; i++)
        for(int j = 0; j < n_basis; j++)
            basis[i][j] = raw_basis[order_vector[i]][order_vector[j]];

    // prerare random direction mutation
    int n_dir = n_tangent + n_basis;

    for(int i = 0; i < n_basis; i++)
        dir[i] = basis[i];

    for(int i = 0; i < n_tangent; i++)
        dir[n_basis + i] = tangent_cone[i];

    for(int i = 0; i < n_basis; i++)
    {
        index_vector[i] = i;
        dir_sign[i] = 1.;
    }
    
    int i_base = n_basis;
    for(int i = i_base; i < i_base + n_basis; i++)
    {
        index_vector[i] = i - i_base;
        dir_sign[i] = -1.;
    }
    
    i_base += n_basis;
    for(int i = i_base; i < i_base + n_tangent; i++)
    {
        index_vector[i] = i - i_base;
        dir_sign[i] = 1.;
    }
    
    i_base += n_tangent;
    for(int i = i_base; i < i_base + n_tangent; i++)
    {
        index_vector[i] = i - i_base;
        dir_sign[i] = -1.;
    }
    
    int n_dir_total = 2 * n_dir;
    for(int i = 0; i < n_dir_total; i++)
        order_vector[i] = i;
    std::random_shuffle(order_vector, order_vector + n_dir_total);

    // finally, mutate
    double success = false;
    for(int i = 0; i < n_dir_total; i++)
    {
        int k = index_vector[order_vector[i]];
        Individual direction = dir_sign[k] * dir[k];
        
        child = parent + ma_step_size * scale * direction;

        if(feasible(child))
        {
            success = true;
            break;
        }
    }

    if( !success )
    {
        child = parent;
        if(options.verbose)
            std::cout << "mutation failed at x = " << parent << "\n";
    }
}


// all the rest ===================================================================================

template<int N> typename GA<N>::Individual GA<N>::random_individual
    (const typename Individual &lower_boundary, const typename Individual &upper_boundary)
{
    Individual res;
    for(int i = 0; i < N; ++i)
        res[i] = lower_boundary[i] + dist01(rnd_generator) * (upper_boundary[i] - lower_boundary[i]);
    
    return res;
}


template<int N> bool GA<N>::feasible(const Individual &x)
{
    for(int i = 0; i < N; i++)
    {
        if(x[i] > upper_boundary[i] || x[i] < lower_boundary[i])
            return false;
    }

    return true;
}

#endif