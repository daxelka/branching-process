import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from collections import Counter
import random
import pandas as pd

# Creating a network
# network parameters
p_in = 8/10
p_out = 1/20
sizes = [5, 5]
probs = [[p_in, p_out], [p_out, p_in]]

g = nx.stochastic_block_model(sizes, probs, seed=0)

# preparing for drawing
color_map = ['yellow' if g.nodes[node]['block'] == 0 else 'pink' for node in g.nodes()]
positions = nx.spring_layout(g)
# nx.draw_networkx(g, positions, node_color=color_map, with_labels=False)
nx.draw_circular(g, node_color=color_map, with_labels=True)
plt.show()

# Branching process simulation
# setting branching process parameters
prob_infection = 0.5
n_sim = 2
# initial seed
nodes_block0 = [x for x,y in g.nodes(data=True) if y['block']==0]
seed = random.sample(nodes_block0, 1)
print(seed)
active_nodes = seed
rng = np.random.default_rng()
results_dic = {'0': [seed]}

for i in range(2):
    results = []
    for active_node in active_nodes:
        # take IDs of all nodes that are connected by a edge which received a possitive outcome in Bernoulli trials
        neighbours = list(g.neighbors(active_node))
        if neighbours:
            infections = rng.binomial(1, prob_infection, (1, len(neighbours)))
            infections_np = np.array(infections[0])
            neighbours_np = np.array(neighbours)
            infected_nodes = list(neighbours_np[infections_np > 0])
            if infected_nodes:
                results.append([active_node, infected_nodes])

    results_dic[str(i)] = results
    active_nodes = infected_nodes

print(results_dic)
