import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from branching_process import BranchingProcess

# branching process parameters
initial_seed = 1
max_generation = 100
lambda_param = 6
probability_success = 0.05

# simulation parameters
n_simulations = 5000

# initiate branching process
bp = BranchingProcess(initial_seed, max_generation, lambda_param, probability_success)

# run branching process n_simulation types
sim_results = bp.run(n_simulations)

# to count the total frequencies for each values of total_infections
infections_frequencies = sim_results.total_infections.value_counts()/n_simulations
frequencies_sorted = infections_frequencies.sort_values()
print(frequencies_sorted)

plt.plot(list(frequencies_sorted.index), list(frequencies_sorted.values))
# plt.yscale('log')
plt.show()
