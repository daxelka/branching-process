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

reinfections_mtbp = mtbt.reinfection_probability(sim_results, 'new_infections_1')

# from pgf analysis
lin, lout, pin = 8, 2, 0.09

pgf = PGFAnalysis()
reinfections_pgf = []
t_list = list(range(50))
for t in t_list:
    reinfections_pgf.append(pgf.reinfection_probability(t, lin, lout, pin))

# plotting results
plt.plot(reinfections_mtbp.gens, reinfections_mtbp.probs)
plt.scatter(t_list, reinfections_pgf)
plt.show()
