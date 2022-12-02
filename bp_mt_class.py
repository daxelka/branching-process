import pandas as pd
import numpy as np
import math


def get_lifetime_distribution(results):
    # to count the total frequencies for each values of total_infections
    max_generation_bp = results.generation.value_counts() / results.simulation.size
    max_generation_bp = max_generation_bp.sort_values()
    lifetime_distribution = pd.DataFrame({'gens': np.array(max_generation_bp.index),
                                          'probs': list(max_generation_bp.values)})
    return lifetime_distribution


class BranchingProcessMultiType:

    def __init__(self, seed1=1, seed2=0,
                 lambda_in=6, lambda_out=6,
                 probability_in=0.05, probability_out=0.5):
        self.max_generation = int(50)
        self.rng = np.random.default_rng()
        self.seed_1 = seed1
        self.seed_2 = seed2
        self.lambda_in = lambda_in
        self.lambda_out = lambda_out
        self.p_in = probability_in
        self.p_out = probability_out

    def branching(self):
        """
        :return: DataFrame with the number of generation and the number
                 of new offsprings produced at that generation until zero offsprings produced
        """
        results = pd.DataFrame({'generation': list(range(self.max_generation)),
                                'new_infections_1': 0,
                                'new_infections_2': 0})
        # initial seed size for zero generation
        results['new_infections_1'][0] = self.seed_1
        results['new_infections_2'][0] = self.seed_2
        population_1 = self.seed_1
        population_2 = self.seed_2
        # get children for the first generation
        generation = 1
        offsprings_1 = np.sum(self.infection_distribution(population_1, self.lambda_in, self.p_in)) \
                       + np.sum(self.infection_distribution(population_2, self.lambda_out, self.p_out))
        offsprings_2 = np.sum(self.infection_distribution(population_1, self.lambda_out, self.p_out)) \
                       + np.sum(self.infection_distribution(population_2, self.lambda_in, self.p_in))
        # save the results of the first generation
        results['new_infections_1'][generation] = offsprings_1
        results['new_infections_2'][generation] = offsprings_2

        # simulate the following simulations
        while generation < (self.max_generation - 1) and (offsprings_1 + offsprings_2) != 0:
            generation += 1
            offsprings_1 = np.sum(self.infection_excess_distribution(population_1, self.lambda_in, self.p_in)) \
                           + np.sum(self.infection_excess_distribution(population_2, self.lambda_out, self.p_out))
            offsprings_2 = np.sum(self.infection_excess_distribution(population_1, self.lambda_out, self.p_out)) \
                           + np.sum(self.infection_excess_distribution(population_2, self.lambda_in, self.p_in))
            # results['new_infections_type1'][generation] = offsprings_1
            # results['new_infections_type2'][generation] = offsprings_2
            results.loc[generation, 'new_infections_1'] = offsprings_1  # may be a faster way
            results.loc[generation, 'new_infections_2'] = offsprings_2  # may be a faster way
            population_1 = offsprings_1
            population_2 = offsprings_2

        # cut off the results at the first generation with  zero offsprings
        if (offsprings_1 + offsprings_2) == 0:
            results = results.loc[:generation, :]

        # calculate total infections in each generation
        results['total_infections_1'] = np.cumsum(results['new_infections_1'])
        results['total_infections_2'] = np.cumsum(results['new_infections_2'])
        results['total_infections'] = results['total_infections_1'] + results['total_infections_2']
        return results

    def infection_distribution(self, n_infections, lambda_param, probability_success):
        n_links = self.rng.poisson(lambda_param, size=n_infections)
        new_infections = self.rng.binomial(n_links, probability_success, (1, n_infections))
        return new_infections

    def infection_excess_distribution(self, n_infections, lambda_param, probability_success):
        n_links = self.excess_poisson_generator(lambda_param, n_infections)
        new_infections = self.rng.binomial(n_links, probability_success, (1, n_infections))
        return new_infections

    def excess_poisson_generator(self, lambda_param, n):
        """
        :param n: number of iid to sample
        :param lambda_param: the mean of the excess poisson distribution
        :return: n number of iid sampled from the excess poisson distribution
        """
        return np.random.choice(list(range(100)),
                                size=n,
                                replace=True,
                                p=self.excess_poisson_probability(list(range(100)), lambda_param))

    def excess_poisson_probability(self, actual, mean):
        if type(actual) is list:
            excess_probability = [(k + 1) / mean * self.poisson_probability(k + 1, mean) for k in actual]
        elif type(actual) is int:
            excess_probability = (actual + 1) / mean * self.poisson_probability(actual + 1, mean)
        else:
            raise ValueError('first parameter is neither list nor integer')
        return excess_probability

    @staticmethod
    def poisson_probability(actual, mean):
        # naive:   math.exp(-mean) * mean**actual / factorial(actual)
        # iterative, to keep the components from getting too large or small:
        p = math.exp(-mean)
        for i in range(actual):
            p *= mean
            p /= i + 1
        return p

    def run(self, n_simulations):
        """
        :param n_simulations: number of times to run the branching process
        :return: DataFrame with the # of simulation, total infections, and the last generation
        """
        sim_results = pd.DataFrame({'simulation': list(range(n_simulations)),
                                    'generation': 0,
                                    'total_infections_1': 0,
                                    'total_infections_2': 0,
                                    'total_infections': 0})

        for simulation in range(n_simulations):
            results = self.branching()
            sim_results.loc[simulation, 'generation'] = results['generation'].iloc[-1]
            sim_results.loc[simulation, 'total_infections_1'] = results['total_infections_1'].iloc[-1]
            sim_results.loc[simulation, 'total_infections_2'] = results['total_infections_2'].iloc[-1]
            sim_results.loc[simulation, 'total_infections'] = results['total_infections'].iloc[-1]
        return sim_results