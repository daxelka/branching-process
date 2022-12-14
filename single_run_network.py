import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import time
from bp_mt_class import BranchingProcessMultiType
from bp_network_class import BranchingProcessNetwork
from bp_network_class import get_max_generation
from bp_network_class import get_average_number_offspring
from bp_network_class import get_total_infections

# Branching process parameters
# Mean degree within and cross communities in network
lambda_in, lambda_out = 9, 4   # Poisson lambda parameters
prob_infection = 0.05

# Network based simulation
community_size = 1000  # number of nodes in each communities
p_in, p_out = lambda_in/community_size, lambda_out/community_size

# number of simulations
n_sim = 2000

bp_network = BranchingProcessNetwork(p_in, p_out, [community_size, community_size], prob_infection)

results_network = bp_network.simulations(n_sim)

total_infections = get_total_infections(results_network, community_size, block='total')
infections_block0 = get_total_infections(results_network, community_size, block='block0')
infections_block1 = get_total_infections(results_network, community_size, block='block1')

# for sim in range(len(results_network)):
#     # plt.plot(list(range(len(total_infections[str(sim)]))), total_infections[str(sim)])
#     plt.plot(list(range(len(infections_block0[str(sim)]))), infections_block0[str(sim)])
#     plt.plot(list(range(len(infections_block1[str(sim)]))), infections_block1[str(sim)])
#     plt.plot(list(range(len(infections_block1[str(sim)]))), total_infections[str(sim)])
# plt.xlabel('generation')
# plt.ylabel('total infections')
# plt.legend(['block 0', 'block 1', 'total'], loc='upper left')
# # set axis ticks to integers only
# plt.gca().yaxis.get_major_locator().set_params(integer=True)
# plt.show()

print('done')

max_generation = 20

def get_cascades_distribution_np(infections, max_generation):
    n_sim = len(infections)
    cascades_size_distribution_np = np.zeros((max_generation,), dtype='float')
    for i in range(n_sim):
        cascades_size_distribution_np = cascades_size_distribution_np + np.array(
            list(infections[str(i)]) + [0] * (max_generation - len(infections[str(i)])))

    return cascades_size_distribution_np / n_sim

def get_cascades_distribution_list(infections, max_generation):
    n_sim = len(infections)
    cascades = []
    for i in range(len(infections)):
        cascades = cascades + [list(infections[str(i)]) + [0]*(max_generation - len(infections[str(i)]))]
    return np.sum(np.array(cascades, dtype='int_'), axis=0) / n_sim

t1_np = time.time()
cascades_size_distribution_np = get_cascades_distribution_np(total_infections, max_generation)
t2_np = time.time()
print('elapsed time numpy array ' + str(t2_np-t1_np))
# t1_list = time.time()
# cascades_size_distribution_list = get_cascades_distribution_np(total_infections, max_generation)
# t2_list = time.time()
# print('elapsed time list ' + str(t2_list - t1_list))

plt.plot(np.array(range(max_generation)), cascades_size_distribution_np)
plt.gca().xaxis.get_major_locator().set_params(integer=True)
plt.show()

# Cascade distribution in different communities
cascades_block0 = get_cascades_distribution_np(infections_block0, max_generation)
cascades_block1 = get_cascades_distribution_np(infections_block1, max_generation)
plt.plot(np.array(range(max_generation)), cascades_block0)
plt.plot(np.array(range(max_generation)), cascades_block1)
plt.gca().xaxis.get_major_locator().set_params(integer=True)
plt.show()


# max_generation = get_max_generation(results_network)
#
# # lifetime distribution
# vals, prob = np.unique(max_generation, return_counts=True)
#
# # lifetime distribution visualisation
# # bp simulation results
# ax = plt.figure().gca()
# # network based results
# plt.scatter(vals+1, prob/len(max_generation), c='b', label='-1')
# plt.xlabel('lifetime', fontsize=14)
# plt.ylabel('probability', fontsize=14)
# plt.yscale('log')
# ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
# plt.show()