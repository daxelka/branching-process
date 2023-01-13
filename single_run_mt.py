import pandas as pd
import numpy as np
import math
import time
import matplotlib.pyplot as plt
from bp_mt_class import BranchingProcessMultiType
from bp_mt_class import get_hazart_function


# branching process parameters
seed_1, seed_2 = 1, 0
lambda_in, lambda_out = 8, 2
probability_in, probability_out = 0.06, 0.06

# simulation parameters
n_simulations = 10

# initiate branching process
bp = BranchingProcessMultiType(seed_1, seed_2,
                               lambda_in, lambda_out,
                               probability_in, probability_out)


# run branching process n_simulation types
# t_start = time.time()
# sim_results = bp.branching()
# print('elapsed time:', time.time() - t_start)
# print(sim_results.iloc[:, 1:5])

# run branching process n_simulation types
t_start = time.time()
sim_results = bp.run(n_simulations)
print('elapsed time:', time.time() - t_start)
# print(sim_results.iloc[1:5, 1:5])

# to count the total frequencies for each values of total_infections
frequencies = sim_results.total_infections.value_counts()/n_simulations
frequencies = frequencies.sort_values()
frequencies_1 = sim_results.total_infections_1.value_counts()/n_simulations
frequencies_1 = frequencies_1.sort_values()
frequencies_2 = sim_results.total_infections_2.value_counts()/n_simulations
frequencies_2 = frequencies_2.sort_values()
# print(frequencies_1)
# print(frequencies_2)
# print(frequencies)

# plt.scatter(list(frequencies.index), list(frequencies.values))
# # plt.plot(list(frequencies_1.index), list(frequencies_1.values))
# # plt.plot(list(frequencies_2.index), list(frequencies_2.values))
# plt.show()

max_generation_counts, max_generation_values, hazart_function = get_hazart_function(sim_results)

# comparing with Dave's results
df = pd.read_csv('Davids_code/full_extin_dist_v3.csv')

plt.scatter(hazart_function.gens-1, hazart_function.probs, marker='x')
plt.scatter(df.t[df['pin']==probability_in], df.hazard_b[df['pin']==probability_out])
plt.gca().xaxis.get_major_locator().set_params(integer=True)
plt.xlabel('generation')
plt.ylabel('probability')
plt.title('hazard function, p_in: ' + str(probability_in))
plt.show()

