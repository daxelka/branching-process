import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import random
import pandas as pd

# Creating a network
# network parameters
p_in = 9/1000
p_out = 4/1000
sizes = [500, 500]
probs = [[p_in, p_out], [p_out, p_in]]

g = nx.stochastic_block_model(sizes, probs, seed=0)
nodes_block0 = [x for x, y in g.nodes(data=True) if y['block'] == 0]

# # preparing for drawing
# color_map = ['yellow' if g.nodes[node]['block'] == 0 else 'pink' for node in g.nodes()]
# positions = nx.spring_layout(g)
# # nx.draw_networkx(g, positions, node_color=color_map, with_labels=False)
# nx.draw_circular(g, node_color=color_map, with_labels=True)
# plt.show()

# Branching process simulation
# setting branching process parameters
prob_infection = 0.05
max_generation = 10
# initial seed from block 0 always

simulations = {}

for i in range(5):
    seed = random.sample(nodes_block0, 1)
    active_nodes = seed
    rng = np.random.default_rng()
    results_dic = {'0': {}}
    generation = 0
    while active_nodes and generation <= max_generation:
        results = {}
        for active_node in active_nodes:
            # take IDs of all nodes that are connected by a edge which received a possitive outcome in Bernoulli trials
            neighbours = list(g.neighbors(active_node))
            if neighbours:
                infections = rng.binomial(1, prob_infection, (1, len(neighbours)))
                infections_np = np.array(infections[0])
                neighbours_np = np.array(neighbours)
                infected_nodes = list(neighbours_np[infections_np > 0])
                results['active_nodes'] = active_nodes
                results['offspring'] = infected_nodes
                results['n_active'] = len(active_nodes)
                results['average_n_offspring'] = np.sum(infections)/len(active_nodes)
            # remove current active node
            g.remove_node(active_node)

        results_dic[str(generation)] = results
        active_nodes = infected_nodes
        generation += 1
    # simulations[str(i)] = {'max_generation': generation-1, 'data': results_dic}
    simulations[str(i)] = results_dic

print(simulations)


def get_max_generation(data):
    return [len(data[str(j)])-1 for j in range(len(data))]


# lifetime distribution
lifetime_distribution = get_max_generation(simulations)
print(lifetime_distribution)

plt.hist(lifetime_distribution)
plt.show()

# tt = [simulations[str(i)]['max_generation'] for i in range(len(simulations))]
# ttt = [len(simulations[str(i)]['data'])-1 for i in range(len(simulations))]
# plt.hist(tt)
# plt.show()