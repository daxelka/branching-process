import pandas as pd
import numpy as np
import math
import time
import matplotlib.pyplot as plt
from bp_mt_class import BranchingProcessMultiType
from bp_mt_class import get_hazart_function

# load Dave's results
df = pd.read_csv('Davids_code/full_extin_dist_v3.csv')

# branching process parameters
seed_1, seed_2 = 1, 0
lambda_in, lambda_out = 8, 2
probabilities = list(set(df.pin))

# simulation parameters
n_simulations = 100

hazart_functions = []

# initiate branching process
for probability in probabilities:
    probability_in, probability_out = probability, probability
    # initiate branching process
    bp = BranchingProcessMultiType(seed_1, seed_2,
                               lambda_in, lambda_out,
                               probability_in, probability_out)


    # run branching process n_simulation types
    t_start = time.time()
    sim_results = bp.run(n_simulations)
    print('elapsed time:', time.time() - t_start)

    max_generation_counts, max_generation_values, hazart_function = get_hazart_function(sim_results)

    # record results
    hazart_functions.append(hazart_function)

# visualisation
# for count, probability in enumerate(probabilities):
#     hazart_function = hazart_functions[count]
#     plt.scatter(hazart_function.gens-1, hazart_function.probs, marker='x')
#     plt.scatter(df.t[df['pin']==probability], df.hazard_b[df['pin']==probability])
#     plt.gca().xaxis.get_major_locator().set_params(integer=True)
#     plt.xlabel('generation')
#     plt.ylabel('probability')
#     plt.title('hazard function, p_in: ' + str(probability))
# plt.show()


# fig, axs = plt.subplots(2, 2)
#
# axs[0, 0].scatter(hazart_functions[0].gens-1, hazart_functions[0].probs, marker='x')
# axs[0, 0].scatter(df.t[df['pin']==probabilities[0]], df.hazard_b[df['pin']==probabilities[0]])
# axs[0, 0].set_title('p_in: ' + str(probabilities[0]))
#
# axs[0, 1].scatter(hazart_functions[1].gens-1, hazart_functions[1].probs, marker='x')
# axs[0, 1].scatter(df.t[df['pin']==probabilities[1]], df.hazard_b[df['pin']==probabilities[1]])
# axs[0, 1].set_title('p_in: ' + str(probabilities[1]))
#
# axs[1, 0].scatter(hazart_functions[2].gens-1, hazart_functions[2].probs, marker='x')
# axs[1, 0].scatter(df.t[df['pin']==probabilities[2]], df.hazard_b[df['pin']==probabilities[2]])
# axs[1, 0].set_title('p_in: ' + str(probabilities[2]))
#
# for ax in axs.flat:
#     ax.set(xlabel='generation', ylabel='probability')
#
# # Hide x labels and tick labels for top plots and y ticks for right plots.
# for ax in axs.flat:
#     ax.label_outer()
#
# plt.show()

fig, axs = plt.subplots(3)

axs[0].scatter(hazart_functions[0].gens-1, hazart_functions[0].probs, marker='x')
axs[0].scatter(df.t[df['pin']==probabilities[0]], df.hazard_b[df['pin']==probabilities[0]])
axs[0].set_title('p_in: ' + str(probabilities[0]))

axs[1].scatter(hazart_functions[1].gens-1, hazart_functions[1].probs, marker='x')
axs[1].scatter(df.t[df['pin']==probabilities[1]], df.hazard_b[df['pin']==probabilities[1]])
axs[1].set_title('p_in: ' + str(probabilities[1]))

axs[2].scatter(hazart_functions[2].gens-1, hazart_functions[2].probs, marker='x')
axs[2].scatter(df.t[df['pin']==probabilities[2]], df.hazard_b[df['pin']==probabilities[2]])
axs[2].set_title('p_in: ' + str(probabilities[2]))

for ax in axs.flat:
    ax.set(xlabel='generation', ylabel='probability')

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()

plt.show()
