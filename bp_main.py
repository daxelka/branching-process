import matplotlib.pyplot as plt
import numpy as np
from bp_class import BranchingProcess
from bp_class import get_max_generation
from bp_class import get_average_number_offspring

# p_in = 9/10
# p_out = 4/10
# sizes = [10, 10]
# prob_infection = 0.1
# n_sim = 10

p_in = 9/1000
p_out = 4/1000
sizes = [500, 500]
prob_infection = 0.05
n_sim = 100

bp = BranchingProcess(p_in, p_out, sizes, prob_infection)

results = bp.simulations(n_sim)
print(results)
max_generation = get_max_generation(results)

print(get_max_generation(results))
# print(get_average_number_offspring(results))

plt.hist(max_generation)
plt.show()

vals, prob = np.unique(max_generation, return_counts=True)

plt.scatter(vals, prob/len(max_generation))
plt.xlabel('lifetime', fontsize=12)
plt.ylabel('probability', fontsize=12)
plt.yscale('log')
plt.show()


