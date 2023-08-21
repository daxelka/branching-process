import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import numpy as np
import time
from bp_mt_class import BranchingProcessMultiType
from mtbp_analysis_class import MTBPAnalysis


# branching process parameters
seed_1, seed_2 = 1, 0
lambda_in, lambda_out = 8, 2

# simulation parameters
n_simulations = 100000
probability = 0.08
probability_in, probability_out = probability, probability

# initiate branching process
bp = BranchingProcessMultiType(seed_1, seed_2,
                               lambda_in, lambda_out,
                               probability_in, probability_out)

# run branching process n_simulation types
t_start = time.time()
sim_results = bp.run_v2(n_simulations)
print('elapsed time:', time.time() - t_start)

# write to file
sim_results.to_csv('data/sim_results/p_' + str(probability) + '_n_sim_' + str(n_simulations) + '.csv')