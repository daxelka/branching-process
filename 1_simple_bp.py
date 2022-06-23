import pandas as pd
import numpy as np
import math
import random
import matplotlib.pyplot as plt

rng = np.random.default_rng()


def branching(n_infections=1, max_generation=100, lambda_param=6, probability_success=0.05):
    generation = 0
    results = pd.DataFrame({'generation': list(range(max_generation)),
                            'new_infections': 0})
    # initial seed size for zero generation
    results['new_infections'][0] = 1
    # get children for the first generation
    offsprings = np.sum(infection_distribution(n_infections, lambda_param, probability_success))
    # save the results of the first generation
    generation = 1
    results['new_infections'][generation] = offsprings
    # simulate the following simulations
    while generation < (max_generation-1) and offsprings != 0:
        offsprings = np.sum(infection_excess_distribution(n_infections, lambda_param, probability_success))
        generation += 1
        # sim_results['new_infections'][generation] = offsprings
        results.loc[generation, 'new_infections'] = offsprings  # may be a faster way
        n_infections = offsprings

    # cut off the results at the first generation with  zero offsprings
    if offsprings == 0:
        results = results.loc[:generation, :]

    # calculate total infections in each generation
    results['total_infections'] = np.cumsum(results['new_infections'])
    return results


def infection_distribution(n_infections, lambda_param, probability_success):
    n_links = rng.poisson(lambda_param, size=n_infections)
    new_infections = rng.binomial(n_links, probability_success, (1, n_infections))
    return new_infections


def infection_excess_distribution(n_infections, lambda_param, probability_success):
    n_links = excess_poisson_generator(lambda_param, n_infections)
    new_infections = rng.binomial(n_links, probability_success, (1, n_infections))
    return new_infections


def excess_poisson_generator(lambda_param, n):
    """
    :param n: number of iid to sample
    :param lambda_param: the mean of the excess poisson distribution
    :return: n number of iid sampled from the excess poisson distribution
    """
    return np.random.choice(list(range(100)),
                            size=n,
                            replace=True,
                            p=excess_poisson_probability(list(range(100)), lambda_param))


def excess_poisson_probability(actual, mean):
    if type(actual) is list:
        excess_probability = [(k + 1) / mean * poisson_probability(k + 1, mean) for k in actual]
    elif type(actual) is int:
        excess_probability = (actual + 1) / mean * poisson_probability(actual + 1, mean)
    return excess_probability


def poisson_probability(actual, mean):
    # naive:   math.exp(-mean) * mean**actual / factorial(actual)
    # iterative, to keep the components from getting too large or small:
    p = math.exp(-mean)
    for i in range(actual):
        p *= mean
        p /= i + 1
    return p

def run(n_simulations):
    sim_results = pd.DataFrame({'simulation': list(range(n_simulations)),
                                'generation': 0,
                                'total_infections': 0})

    for simulation in range(n_simulations):
        results = branching()
        sim_results.loc[simulation, 'generation'] = results['generation'].iloc[-1]
        sim_results.loc[simulation, 'total_infections'] = results['total_infections'].iloc[-1]
    return sim_results

# run simulations
n_simulations = 5000
sim_results = run(n_simulations)

# to count the total frequencies for each values of total_infections
infections_frequencies = sim_results.total_infections.value_counts()/n_simulations
frequencies_sorted = infections_frequencies.sort_values()
# print(total_infections_frequencies)
print(frequencies_sorted)

# plt.scatter(frequencies_sorted.index, frequencies_sorted.values)
plt.plot(list(frequencies_sorted.index), list(frequencies_sorted.values))
# plt.yscale('log')
plt.show()
