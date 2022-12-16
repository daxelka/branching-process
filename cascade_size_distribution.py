import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import time

import pandas as pd

from bp_mt_class import BranchingProcessMultiType
from bp_network_class import BranchingProcessNetwork
from bp_network_class import get_max_generation
from bp_network_class import get_average_number_offspring
from bp_network_class import get_total_infections

# Branching process parameters
# Mean degree within and cross communities in network
lambda_in, lambda_out = 8, 8  # Poisson lambda parameters
prob_infection = 0.05

# Network based simulation
community_size = 500  # number of nodes in each communities
p_in, p_out = lambda_in / community_size, lambda_out / community_size

# number of simulations
n_sim = 500

bp_network = BranchingProcessNetwork(p_in, p_out, [community_size, community_size], prob_infection)

results_network = bp_network.simulations(n_sim)

total_infections = get_total_infections(results_network, community_size, block='total')
infections_block0 = get_total_infections(results_network, community_size, block='block0')
infections_block1 = get_total_infections(results_network, community_size, block='block1')

print('done')


def get_cascades_distribution_np(infections):
    n_sim = len(infections)
    cascades_i = np.zeros((n_sim,))
    for i in range(n_sim):
        cascades_i[i,] = np.sum(np.array(list(infections[str(i)])), axis=0)

    cascade_prob = np.array([len(cascades_i[cascades_i >= n]) for n in range(1, int(max(cascades_i)) + 1)]) / n_sim
    cascade_vals = list(range(1, int(max(cascades_i)) + 1))
    cascade_distribution = pd.DataFrame({'prob': cascade_prob,
                                         'vals': cascade_vals})
    return cascade_distribution


t1_np = time.time()
cascades_distribution_total = get_cascades_distribution_np(total_infections)
cascades_distribution_block0 = get_cascades_distribution_np(infections_block0)
cascades_distribution_block1 = get_cascades_distribution_np(infections_block1)
t2_np = time.time()
print('elapsed time numpy array ' + str(t2_np - t1_np))

plt.plot(cascades_distribution_total.vals, cascades_distribution_total.prob)
plt.plot(cascades_distribution_block0.vals, cascades_distribution_block0.prob)
plt.plot(cascades_distribution_block1.vals, cascades_distribution_block1.prob)
plt.xlabel('cascade size')
plt.ylabel('probability')
# plt.yscale('log')
# plt.gca().xaxis.get_major_locator().set_params(integer=True)
plt.legend(['total', 'block0', 'block1'], loc="upper right")
plt.title('l_in: ' + str(lambda_in) + ' l_out: ' + str(lambda_out) + ' p_in: ' + str(prob_infection) )
plt.show()

