import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import numpy as np
import time
from bp_mt_class import BranchingProcessMultiType
from bp_network_class import BranchingProcessNetwork
from bp_network_class import get_max_generation
from bp_network_class import get_lifetime_distribution as network_get_lifetime_distribution
from bp_network_class import get_average_number_offspring
from bp_mt_class import get_lifetime_distribution as mt_get_lifetime_distribution

"""
runs a series of simulation on a multitype branching process on network and on an approximate model
with the same parameters
calculates and compares lifetime distributions of two models
"""


# Branching process parameters
# Mean degree within and across communities in network
lambda_in, lambda_out = 8, 2   # Poisson lambda parameters
prob_infection = 0.05

# Simulation parameters
n_simulations = 500

# Network based simulation
community_size = 500  # number of nodes in each communities
p_in, p_out = lambda_in/community_size, lambda_out/community_size

bp_network = BranchingProcessNetwork(p_in, p_out, [community_size, community_size], prob_infection)

results_network = bp_network.simulations(n_simulations)
lifetime_distribution_network = network_get_lifetime_distribution(results_network, extend_to=21)

# BP simulation
# branching process parameters
seed_1, seed_2 = 1, 0

# initiate branching process
bp = BranchingProcessMultiType(seed_1, seed_2,
                               lambda_in, lambda_out,
                               prob_infection, prob_infection)

# run branching process n_simulation types
t_start = time.time()
results_bp = bp.run(n_simulations)
print('elapsed time:', time.time() - t_start)

mt_lifetime_distribution = mt_get_lifetime_distribution(results_bp, extend_to=21)

# comparing with Dave's results
df = pd.read_csv('Davids_code/full_extin_dist.csv')

# lifetime distribution visualisation
ax = plt.figure().gca()
# bp simulation results
plt.scatter(mt_lifetime_distribution.gens-1, mt_lifetime_distribution.probs, c='r', label='1')
# network based results
plt.scatter(lifetime_distribution_network.gens, lifetime_distribution_network.probs, c='b', label='-1')
# Davids results
# plt.scatter(df.t[df['pin']==0.04], df.S_b[df['pin']==0.04])
# labels
plt.xlabel('generation', fontsize=14)
plt.ylabel('probability', fontsize=14)
plt.yscale('log')
plt.legend(["approximation", "network", "David's S_b"], loc="upper right")
plt.title('lifetime distribution, l_im: ' + str(lambda_in) + ' m l_out: ' + str(lambda_out) + ' p_in: ' + str(prob_infection))
# enforcing integer ticks at the axis
ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
plt.show()

