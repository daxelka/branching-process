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
probability = 0.04

probability_in, probability_out = probability, probability

# initiate branching process
bp = BranchingProcessMultiType(seed_1, seed_2,
                           lambda_in, lambda_out,
                           probability_in, probability_out)


# run branching process n_simulation types
t_start = time.time()
sim_results = bp.run_v2(n_simulations)
print('elapsed time:', time.time() - t_start)

# sim_results = pd.DataFrame({'generation': [0,1, 0,1,2, 0,1,2,3,],
#                             'new_infections_1': [1,0, 1,2,0, 1,2,0,0],
#                             'new_infections_2': [0,0, 0,1,0, 0,0,1,0]})


def hazard_function_in_communities(data):
    max_generation = max(data.generation)
    hazard_prob_1 = [0]
    hazard_prob_2 = [0]
    hazard_prob_both = [0]
    for gen in range(1,max_generation):
        data_current_gen = data[data.generation == gen]
        data_previous_gen = data[data.generation == gen - 1]

        nonzeros_previous_gen_both_communities = np.count_nonzero(
            data_previous_gen.new_infections_1 + data_previous_gen.new_infections_2)
        nonzeros_current_gen_both_communities = np.count_nonzero(
            data_current_gen.new_infections_1 + data_current_gen.new_infections_2)
        nonzeros_current_gen_1 = np.count_nonzero(data_current_gen.new_infections_1)
        nonzeros_current_gen_2 = np.count_nonzero(data_current_gen.new_infections_2)

        hazard_prob_1.append(1 - nonzeros_current_gen_1/nonzeros_previous_gen_both_communities)
        hazard_prob_2.append(1 - nonzeros_current_gen_2 / nonzeros_previous_gen_both_communities)
        hazard_prob_both.append(1 - nonzeros_current_gen_both_communities / nonzeros_previous_gen_both_communities)

    hazard_function = pd.DataFrame({'gens': list(range(max_generation)),
                                    'probs_1': hazard_prob_1,
                                    'probs_2': hazard_prob_2,
                                    'probs_both': hazard_prob_both})
    return hazard_function


hazard_function = hazard_function_in_communities(sim_results)


# visulation

# plt.scatter(hazard_function.gens, hazard_function.probs_both, marker="x")
# plt.scatter(df.t[df.pin==probability], df.hazard_b[df.pin==probability])
# plt.scatter(hazard_function.gens, hazard_function.probs_2)
# plt.scatter(hazard_function.gens, hazard_function.probs_both)
# plt.show()

# colors
CB91_Blue = '#2CBDFE'
CB91_Green = '#47DBCD'
CB91_Pink = '#F3A0F2'
CB91_Purple = '#9D2EC5'
CB91_Violet = '#661D98'
CB91_Amber = '#F5B14C'

color1, color2, color3 = (CB91_Blue, CB91_Green, CB91_Pink)

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

fig, axs = plt.subplots(3)
fig.suptitle("Hazard function, probability of infection: " + str(probability))

fig.set_figheight(6)
fig.set_figwidth(6)

axs[0].scatter(hazard_function.gens, hazard_function.probs_both, marker="x", c=color1)
axs[0].plot(df.t[df.pin==probability], df.hazard_b[df.pin==probability], c=color1)
axs[0].set_title('whole network')
axs[0].legend(['simulation', 'analytical'], frameon=False, loc='lower right')

axs[1].scatter(hazard_function.gens, hazard_function.probs_1, marker="x", c=color2)
axs[1].plot(df.t[df.pin==probability], df.hazard_1[df.pin==probability], c=color2)
axs[1].set_title('community #1')
axs[1].legend(['simulation', 'analytical'], frameon=False, loc='lower right')

axs[2].scatter(hazard_function.gens, hazard_function.probs_2, marker="x", c=color3)
axs[2].plot(df.t[df.pin==probability], df.hazard_2[df.pin==probability], c=color3)
axs[2].set_title('community #2')
axs[2].legend(['simulation', 'analytical'], frameon=False, loc='lower right')
axs[2].xaxis.get_major_locator().set_params(integer=True)

for ax in axs.flat:
    ax.set(xlabel='generation', ylabel='probability')

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()

# fig.savefig('img/hazard_function_p_'+ str(probability)+'.png', dpi=300)
plt.show()


# fig = plt.figure(figsize=(8, 6))
# fig.suptitle("Hazard function, probability of infection: " + str(probability))
#
# axs[0].scatter(hazard_function.gens, hazard_function.probs_both, marker="x", c=color1)
# axs[0].plot(df.t[df.pin==probability], df.hazard_b[df.pin==probability], c=color1)
# axs[0].set_title('whole network')
# axs[0].legend(['analytical', 'simulation'], frameon=False, loc='lower right')
#
# axs[1].scatter(hazard_function.gens, hazard_function.probs_1, marker="x", c=color2)
# axs[1].plot(df.t[df.pin==probability], df.hazard_1[df.pin==probability], c=color2)
# axs[1].set_title('community #1')
# axs[1].legend(['analytical', 'simulation'], frameon=False, loc='lower right')
#
# axs[2].scatter(hazard_function.gens, hazard_function.probs_2, marker="x", c=color3)
# axs[2].plot(df.t[df.pin==probability], df.hazard_2[df.pin==probability], c=color3)
# axs[2].set_title('community #2')
# axs[2].legend(['analytical', 'simulation'], frameon=False, loc='lower right')
# axs[2].xaxis.get_major_locator().set_params(integer=True)
#
# for ax in axs.flat:
#     ax.set(xlabel='generation', ylabel='probability')
#
# # Hide x labels and tick labels for top plots and y ticks for right plots.
# for ax in axs.flat:
#     ax.label_outer()
#
# fig.savefig('img/hazard_function_p_'+ str(probability)+'.png', dpi=300)
# plt.show()