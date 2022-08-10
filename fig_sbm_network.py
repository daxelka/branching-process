import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import random
import pandas as pd

# Creating a network
# network parameters
p_in = 9/10
p_out = 1/10
sizes = [10, 10]
probs = [[p_in, p_out], [p_out, p_in]]

G = nx.stochastic_block_model(sizes, probs, seed=0)
edges = G.edges()
nodes = G.nodes(data=True)
# # create a new attribute 'active'
# opinions_dict = dict(enumerate([True]*len(G)))
# nx.set_node_attributes(G, opinions_dict, 'active')

# list of nodes in block0
nodes_block0 = [x for x, y in G.nodes(data=True) if y['block'] == 0]


# preparing for drawing
fig = plt.figure()
ax = fig.subplots()
# drawing a network
color_map = ['yellow' if G.nodes[node]['block'] == 0 else 'pink' for node in G.nodes()]
positions = nx.spring_layout(G)
nx.draw_networkx(G, positions, node_color=color_map, with_labels=False)
# nx.draw_circular(G, node_color=color_map, with_labels=True)
# plt.axis('off')

# Annotations
plt.rc('text', usetex=True)
ax.text(0.3, 0.8, r'$C_1$', fontsize=20)
ax.text(0.8, 0.3, r'$C_2$', fontsize=20)
plt.axis('off')
plt.show()