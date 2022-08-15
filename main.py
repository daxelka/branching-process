import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import time
from multitype_branching_process import BranchingProcessMultiType
from bp_class import BranchingProcess
from bp_class import get_max_generation
from bp_class import get_average_number_offspring

# Branching process parameters
# Mean degree within and cross communities in network
lambda_in, lambda_out = 9, 4   # Poisson lambda parameters
prob_infection = 0.1

# Simulation parameters
n_simulations = 300

# Network based simulation
community_size = 500  # number of nodes in each communities
p_in, p_out = lambda_in/community_size, lambda_out/community_size

bp = BranchingProcess(p_in, p_out, [community_size, community_size], prob_infection)

results = bp.simulations(n_simulations)
max_generation = get_max_generation(results)

# lifetime distribution
vals, prob = np.unique(max_generation, return_counts=True)

# BP simulation
# branching process parameters
seed_1, seed_2 = 1, 0
probability_in, probability_out = prob_infection, prob_infection

# initiate branching process
bp = BranchingProcessMultiType(seed_1, seed_2,
                               lambda_in, lambda_out,
                               probability_in, probability_out)

# run branching process n_simulation types
t_start = time.time()
sim_results = bp.run(n_simulations)
print('elapsed time:', time.time() - t_start)

# to count the total frequencies for each values of total_infections
# frequencies = sim_results.total_infections.value_counts()/n_simulations
# frequencies = frequencies.sort_values()
max_generation_bp = sim_results.generation.value_counts()/n_simulations
max_generation_bp = max_generation_bp.sort_values()
print(max_generation_bp)

# lifetime distribution visualisation
# bp simulation results
ax = plt.figure().gca()
# plt.scatter(list(frequencies.index), list(frequencies.values), c='r', label='1')
plt.scatter(np.array(max_generation_bp.index), list(max_generation_bp.values), c='r', label='1')
# network based results
plt.scatter(vals+1, prob/len(max_generation), c='b', label='-1')
plt.xlabel('lifetime', fontsize=14)
plt.ylabel('probability', fontsize=14)
plt.yscale('log')
ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
plt.show()
