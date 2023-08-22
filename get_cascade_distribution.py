import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import numpy as np
import time
from bp_mt_class import BranchingProcessMultiType
from mtbp_analysis_class import MTBPAnalysis

probability = 0.02
n_simulations = 100000

# read simulation results
file_name = 'data/sim_results/p_' + str(probability) + '_n_sim_' + str(n_simulations) + '_results.csv'
sim_results = pd.read_csv(file_name)

# calculate cascades distribution
analysis = MTBPAnalysis()
cascades = analysis.cascade_distribution(sim_results)
print('cascades calculated')

# write to file
cascades.to_csv('data/cascades/p_' + str(probability) + '_n_sim_' + str(n_simulations) + '_cascades.csv')
print('saved to a file')