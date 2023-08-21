import pandas as pd
import numpy as np
import math
import time
import matplotlib.pyplot as plt
from bp_mt_class import BranchingProcessMultiType
from mtbp_analysis_class import MTBPAnalysis

# branching process parameters
seed_1, seed_2 = 1, 0
lambda_in, lambda_out = 8, 2

# simulation parameters
n_simulations = 500
probabilities_infection = [0.02, 0.04, 0.06, 0.08, 0.09]
cascades_both_list = []
cascades_1_list = []
cascades_2_list = []

for probability in probabilities_infection:
    probability_in, probability_out = probability, probability

    # initiate branching process
    bp = BranchingProcessMultiType(seed_1, seed_2,
                               lambda_in, lambda_out,
                               probability_in, probability_out)

    # run branching process n_simulation types
    t_start = time.time()
    sim_results = bp.run_v2(n_simulations)
    print('elapsed time:', time.time() - t_start)

    # analysis
    analysis = MTBPAnalysis()

    # analysis.get_max_simulation_id(sim_results)
    cascades_both, cascades_1, cascades_2 = analysis.cascade_distribution(sim_results)
    cascades_both_list.append(cascades_both)
    cascades_1_list.append(cascades_1)
    cascades_2_list.append(cascades_2)


    # print(cascade_size_dist_1)
    # print(cascade_size_dist_2)
    # print(cascade_size_dist_both)
for cascades in cascades_both_list:
    plt.plot(cascades)

plt.xlim([0,50])
plt.gca().yaxis.get_major_locator().set_params(integer=True)
plt.yscale('log')
plt.show()