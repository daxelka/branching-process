import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import numpy as np
import time
from bp_mt_class import BranchingProcessMultiType
from mtbp_analysis_class import MTBPAnalysis

probability = 0.095
n_simulations = int(1e05)

# read simulation results
file_name = 'data/sim_results/p_' + str(probability) + '_n_sim_' + str(n_simulations) + '_results.csv'
sim_results = pd.read_csv(file_name)

# calculate cascades distribution
t_start = time.time()
analysis = MTBPAnalysis()
cascades = analysis.cascade_distribution(sim_results)
print('cascades calculated, elapsed time: ' + str(time.time()-t_start))

# write to file
cascades.to_csv('data/cascades/numerical/p_' + str(probability) + '_n_sim_' + str(n_simulations) + '_cascades.csv')
print('saved to a file')