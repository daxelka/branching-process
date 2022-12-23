import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import time

import pandas as pd

from bp_network_class import BranchingProcessNetwork
from bp_network_class import get_total_infections


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

# Branching process parameters
# Mean degree within and cross communities in network
lambda_in, lambda_out = 8, [0,2,4,6,8]  # Poisson lambda parameters
prob_infection = 0.05

# Network based simulation
community_size = 500  # number of nodes in each communities
# number of simulations
n_sim = 50

results = []
for l_out in lambda_out:
    p_in, p_out = lambda_in / community_size, l_out / community_size
    bp_network = BranchingProcessNetwork(p_in, p_out, [community_size, community_size], prob_infection)
    results_network = bp_network.simulations(n_sim)
    total_infections = get_total_infections(results_network, community_size, block='total')

    cascades_distribution_total = get_cascades_distribution_np(total_infections)
    plt.plot(cascades_distribution_total.vals, cascades_distribution_total.prob)
print('done')
plt.legend(['l_out: '+ str(l) for l in lambda_out])
plt.title('l_in: ' + str(lambda_in) + ' p_in: ' + str(prob_infection))
plt.xlabel('cascade size')
plt.ylabel('probability')
plt.show()


# # Initialise the subplot function using number of rows and columns
# figure, axis = plt.subplots(1, 1)
#
# # For Sine Function
# axis[0, 0].plt.plot(results[0].vals, results[0].prob)
# axis[0, 0].set_title("Sine Function")
#
# # For Cosine Function
# axis[0, 1].plt.plot(results[1].vals, results[1].prob)
# axis[0, 1].set_title("Cosine Function")
#
# plt.show()


