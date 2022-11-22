import pandas as pd
import numpy as np
import math


class BranchingProcess:

    def __init__(self, initial_seed=1, max_generation=50, lambda_param=6, probability_success=0.05):
        self.rng = np.random.default_rng()
        self.initial_seed = initial_seed
        self.max_generation = max_generation
        self.lambda_param = lambda_param
        self.probability_success = probability_success

    def branching(self):
        """
        :return: DataFrame with the number of generation and the number
                 of new offsprings produced at that generation until zero offsprings produced
        """
        results = pd.DataFrame({'generation': list(range(self.max_generation)),
                                'new_infections': 0})
        # initial seed size for zero generation
        results['new_infections'][0] = self.initial_seed
        population = self.initial_seed
        # get children for the first generation
        offsprings = np.sum(self.infection_distribution(population,
                                                        self.lambda_param,
                                                        self.probability_success))
        # save the results of the first generation
        generation = 1
        results['new_infections'][generation] = offsprings
        # simulate the following simulations
        while generation < (self.max_generation - 1) and offsprings != 0:
            offsprings = np.sum(self.infection_excess_distribution(population,
                                                                   self.lambda_param,
                                                                   self.probability_success))
            generation += 1
            # results['new_infections'][generation] = offsprings
            results.loc[generation, 'new_infections'] = offsprings  # may be a faster way
            population = offsprings

        # cut off the results at the first generation with  zero offsprings
        if offsprings == 0:
            results = results.loc[:generation, :]

        # calculate total infections in each generation
        results['total_infections'] = np.cumsum(results['new_infections'])
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
                                    'total_infections': 0})

        for simulation in range(n_simulations):
            results = self.branching()
            sim_results.loc[simulation, 'generation'] = results['generation'].iloc[-1]
            sim_results.loc[simulation, 'total_infections'] = results['total_infections'].iloc[-1]
        return sim_results