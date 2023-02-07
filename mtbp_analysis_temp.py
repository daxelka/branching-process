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
n_simulations = 5
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

# analysis
analysis = MTBPAnalysis()

# analysis.get_max_simulation_id(sim_results)
cascade_size_dist_both, cascade_size_dist_1, cascade_size_dist_2 = analysis.cascade_distribution(sim_results)
print(cascade_size_dist_1)
print(cascade_size_dist_2)
print(cascade_size_dist_both)
