import pandas as pd
import numpy as np
import math
import time
import matplotlib.pyplot as plt
from bp_mt_class import BranchingProcessMultiType
from bp_mt_class import get_hazart_function
from mtbp_analysis_class import MTBPAnalysis

# load Dave's results
df = pd.read_csv('Davids_code/data/full_extin_dist_v4.csv')
# sim_results = pd.read_csv('data/mtbp_pin_0p09.csv')

# branching process parameters
seed_1, seed_2 = 1, 0
lambda_in, lambda_out = 8, 2
probabilities = list(set(df.pin))
#
# simulation parameters
n_simulations = 10000
probability = 0.06
#
probability_in, probability_out = probability, probability

# initiate branching process
bp = BranchingProcessMultiType(seed_1, seed_2,
                           lambda_in, lambda_out,
                           probability_in, probability_out)


# run branching process n_simulation types
t_start = time.time()
sim_results = bp.run_v2(n_simulations)
print('elapsed time:', time.time() - t_start)
sim_results.to_csv('mtbp_pin_' + str(probability)+'.csv')


analysis = MTBPAnalysis()
hazard_function_old = analysis.hazard_function_in_communities(sim_results)
hazard_function_new = analysis.hazard_function_in_communities_v2(sim_results)
# duration_extinction(sim_results)



# visulation

# colors
CB91_Blue = '#2CBDFE'
CB91_Green = '#47DBCD'
CB91_Pink = '#F3A0F2'
CB91_Purple = '#9D2EC5'
CB91_Violet = '#661D98'
CB91_Amber = '#F5B14C'

color1, color2, color3 = (CB91_Blue, CB91_Green, CB91_Pink)
my_marker = '.'

plt.scatter(list(range(len(hazard_function_new))), hazard_function_new.hazard_1, marker=my_marker, c=color1)
plt.scatter(list(range(len(hazard_function_new))), hazard_function_new.hazard_2, marker=my_marker, c=color2)
plt.scatter(list(range(len(hazard_function_new))), hazard_function_new.hazard_both, marker=my_marker, c=color3)

plt.plot(df.t[df.pin==probability], df.hazard_1[df.pin==probability], c=color1)
plt.plot(df.t[df.pin==probability], df.hazard_2[df.pin==probability], c=color2)
plt.plot(df.t[df.pin==probability], df.hazard_b[df.pin==probability], c=color3)
plt.title("Hazard function, probability of infection: " + str(probability))
plt.ylim([0, 1])
plt.xlim([0, 20])
plt.xlabel('generation')
plt.ylabel('hazard function')
plt.legend(['type 1 simulation', 'type 2 simulation', 'both types simulation',
            'type 1 analytical', 'type 2 analytical', 'both types analytical'], frameon=False, loc='upper right')
plt.gca().xaxis.get_major_locator().set_params(integer=True)
plt.show()

# plt.scatter(hazard_function.gens, hazard_function.probs_both, marker="x")
# plt.scatter(df.t[df.pin==probability], df.hazard_b[df.pin==probability])
# plt.scatter(hazard_function.gens, hazard_function.probs_2)
# plt.scatter(hazard_function.gens, hazard_function.probs_both)
# plt.show()

# # colors
# CB91_Blue = '#2CBDFE'
# CB91_Green = '#47DBCD'
# CB91_Pink = '#F3A0F2'
# CB91_Purple = '#9D2EC5'
# CB91_Violet = '#661D98'
# CB91_Amber = '#F5B14C'
#
# color1, color2, color3 = (CB91_Blue, CB91_Green, CB91_Pink)
# my_marker = '.'
#
# plt.rcParams['font.family'] = 'sans-serif'
# plt.rcParams['font.sans-serif'] = ['Arial']
#
# fig, axs = plt.subplots(3)
# fig.suptitle("Hazard function, probability of infection: " + str(probability))
#
# fig.set_figheight(6)
# fig.set_figwidth(6)
#
# axs[0].scatter(hazard_function.gens, hazard_function.probs_both, marker=my_marker, c=color1)
# axs[0].plot(df.t[df.pin==probability], df.hazard_b[df.pin==probability], c=color1)
# axs[0].set_title('whole network')
# axs[0].set_ylim([0, 1])
# axs[0].set_xlim([0, 20])
# axs[0].legend(['simulation', 'analytical'], frameon=False, loc='lower right')
# # axs[0].legend(['simulation', 'analytical'], frameon=False, loc='upper left')
#
# axs[1].scatter(hazard_function.gens, hazard_function.probs_1, marker=my_marker, c=color2)
# axs[1].plot(df.t[df.pin==probability], df.hazard_1[df.pin==probability], c=color2)
# axs[1].set_title('community #1')
# axs[1].set_ylim([0, 1])
# axs[1].set_xlim([0, 20])
# axs[1].legend(['simulation', 'analytical'], frameon=False, loc='lower right')
# # axs[1].legend(['simulation', 'analytical'], frameon=False, loc='upper left')
#
# axs[2].scatter(hazard_function.gens, hazard_function.probs_2, marker=my_marker, c=color3)
# axs[2].plot(df.t[df.pin==probability], df.hazard_2[df.pin==probability], c=color3)
# axs[2].set_title('community #2')
# axs[2].set_ylim([0, 1])
# axs[2].set_xlim([0, 20])
# axs[2].legend(['simulation', 'analytical'], frameon=False, loc='lower right')
# # axs[2].legend(['simulation', 'analytical'], frameon=False, loc='lower left')
#
# axs[2].xaxis.get_major_locator().set_params(integer=True)
#
# for ax in axs.flat:
#     ax.set(xlabel='generation', ylabel='probability')
#
# # Hide x labels and tick labels for top plots and y ticks for right plots.
# for ax in axs.flat:
#     ax.label_outer()
#
# # fig.savefig('img/hazard_function_p_'+ str(probability)+'.png', dpi=300)
# plt.show()

