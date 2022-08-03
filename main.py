import matplotlib.pyplot as plt
import numpy as np
import time
from multitype_branching_process import BranchingProcessMultiType
from bp_class import BranchingProcess
from bp_class import get_max_generation
from bp_class import get_average_number_offspring

# Network based simulation
p_in, p_out = 4.5/1000, 2/1000
sizes = [1000, 1000]
prob_infection = 0.025
n_sim = 300

bp = BranchingProcess(p_in, p_out, sizes, prob_infection)

results = bp.simulations(n_sim)
max_generation = get_max_generation(results)

# lifetime distribution
vals, prob = np.unique(max_generation, return_counts=True)

# BP simulation
# branching process parameters
seed_1, seed_2 = 1, 0
lambda_in, lambda_out = 4.5, 2
probability_in, probability_out = 0.025, 0.025

# simulation parameters
n_simulations = 300

# initiate branching process
bp = BranchingProcessMultiType(seed_1, seed_2,
                               lambda_in, lambda_out,
                               probability_in, probability_out)

# run branching process n_simulation types
t_start = time.time()
sim_results = bp.run(n_simulations)
print('elapsed time:', time.time() - t_start)

# to count the total frequencies for each values of total_infections
frequencies = sim_results.total_infections.value_counts()/n_simulations
frequencies = frequencies.sort_values()
print(frequencies)

# lifetime distribution visualisation
# bp simulation results
# plt.scatter(list(frequencies.index), list(frequencies.values), c='r', label='1')
plt.scatter(np.array(frequencies.index), list(frequencies.values), c='r', label='1')
# network based results
plt.scatter(vals+1, prob/len(max_generation), c='b', label='-1')
plt.xlabel('lifetime', fontsize=14)
plt.ylabel('probability', fontsize=14)
plt.show()
