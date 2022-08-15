import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
import math


# Degree Distribution
def degree_dist(degrees):
    # Degree distribution
    degree_counts = Counter(degrees)
    keys = sorted(degree_counts.keys())

    values = [degree_counts[k] for k in keys]

    # Dergee distribution normalised
    values = [float(i) / sum(values) for i in values]
    return values, keys

# creata a network
# Creating a network
# network parameters
lambda_in = 9
lambda_out = 4
N = 1000
p_in = lambda_in/N
p_out = lambda_out/N
sizes = [N, N]
probs = [[p_in, p_out], [p_out, p_in]]

G = nx.stochastic_block_model(sizes, probs, seed=0)

# degree distribution
# Degrees list
degrees_network = [G.degree(n) for n in G.nodes()]
values_network, keys_network = degree_dist(degrees_network)

# calculate Poison distribution
rng = np.random.default_rng()
degrees_poisson = rng.poisson(lambda_in, size=int(2*N))+rng.poisson(lambda_out, size=int(2*N))

# degrees_poisson = list(rng.poisson(lambda_in, size=int(N/2)))+list(rng.poisson(lambda_out, size=int(N/2)))

values_poisson, keys_poisson = degree_dist(degrees_poisson)
# print(degrees_network)
# print(degrees_poisson)

p_poisson = []
for k in keys_network:
    p_poisson.append((1/ math.factorial(k) * math.exp(-lambda_in) * lambda_in ** k) * (1 / math.factorial(k) * math.exp(-lambda_out) * lambda_out ** k)
)


# visualisations
plt.scatter(keys_network, values_network, label='network, N='+str(2*N))
plt.scatter(keys_poisson, values_poisson, label='Poisson variable')
plt.scatter(keys_network, p_poisson, label='Poisson exact')
plt.xlabel('degree', fontsize=14)
plt.ylabel('probability', fontsize=12)
plt.legend(frameon=False, fontsize=12)
plt.show()