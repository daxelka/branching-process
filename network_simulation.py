import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from collections import Counter
import random

p_in = 6/10
p_out = 1/20

sizes = [20, 20]
probs = [[p_in, p_out], [p_out, p_in]]
g = nx.stochastic_block_model(sizes, probs, seed=0)

# preparing for drawing
color_map = ['yellow' if g.nodes[node]['block'] == 0 else 'pink' for node in g.nodes()]

positions = nx.spring_layout(g)
# nx.draw_networkx(g, positions, node_color=color_map, with_labels=False)
nx.draw_circular(g, node_color=color_map)
plt.show()


# Degree Distribution
def degree_dist(G):
    # Degrees list
    degrees = [G.degree(n) for n in G.nodes()]

    # Degree distribution
    degree_counts = Counter(degrees)
    keys = sorted(degree_counts.keys())

    values = [degree_counts[k] for k in keys]

    # Dergee distribution normalised
    # values = [float(i) / sum(values) for i in values]
    return degrees, values, keys


degrees, values, keys = degree_dist(g)
print(degrees, values, keys)

# setting parameters
probability_in = 0.05
# branching process simulation
active_nodes = []
rng = np.random.default_rng()
while offsprings !=0:
    for active_node in active_nodes:
        # take node degree
        degree = g.degree(active_node)
        # roll dice n_degree times, take n_degree Bernoulli outcomes
        rng.binomial(1, probability_in, degree)
        # take IDs of all nodes that are connected by a edge which received a possitive outcome in Bernoulli trials

        # add them to list of active_nodes
        active_nodes.append()
