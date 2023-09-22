import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mtbp_analysis_class import MTBPAnalysis
from pgf_analysis_class import PGFAnalysis
from plotting_class import Plotting
import time

# from simulation
sim_results = pd.read_csv('sim_data/mtbp_pin_0.09.csv')
mtbt = MTBPAnalysis()

continuous_ext_mtbp = mtbt.continuous_extinction_prob(sim_results, 'new_infections_1')
reinfection_mtbp = mtbt.reinfection_prob(sim_results, 'new_infections_1')
# from pgf analysis
lin, lout, pin = 8, 2, 0.09

pgf = PGFAnalysis()
continuous_ext_pgf = []
reinfection_pgf = []
t_list = list(range(50))
for t in t_list:
    continuous_ext_pgf.append(pgf.continuous_extinction(t, lin, lout, pin))

for t in t_list:
    reinfection_pgf.append(1 - pgf.continuous_extinction(t, lin, lout, pin))

# plotting results
# continuous extinction
plt.scatter(continuous_ext_mtbp.gens, continuous_ext_mtbp.probs, c='r')
plt.plot(t_list, continuous_ext_pgf, c='b')
plt.xlim([0, 10])
plt.show()

# reinfection
plt.scatter(reinfection_mtbp.gens, reinfection_mtbp.probs, c='r')
plt.plot(t_list, reinfection_pgf, c='b')
plt.xlim([0, 10])
plt.show()
