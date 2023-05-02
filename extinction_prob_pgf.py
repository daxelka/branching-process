import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mtbp_analysis_class import MTBPAnalysis
from pgf_analysis_class import PGFAnalysis
from plotting_class import Plotting
import time

lin, lout = 8, 2
pin = 0.09

# PGF analysis
pgf_analysis = PGFAnalysis()

# probability of extinction q(t) = P[N(t)=0]= G_N(t)(0,0)
t_list = list(range(20))
ext_prob_pgf_both = []
ext_prob_pgf_1 = []
ext_prob_pgf_2 = []

for t in t_list:
    ext_prob_pgf_both.append(pgf_analysis.G_N_t(0, 0, t, lin, lout, pin))
    ext_prob_pgf_1.append(pgf_analysis.G_N_t(0, 1, t, lin, lout, pin))
    ext_prob_pgf_2.append(pgf_analysis.G_N_t(1, 0, t, lin, lout, pin))

# hazard function from pgf
hazard_pgf_both = []
hazard_pgf_1 = []
hazard_pgf_2 = []

for t in t_list:
    probs_both, probs_1, probs_2 = pgf_analysis.G_H_t(t, lin, lout, pin)
    hazard_pgf_both.append(probs_both)
    hazard_pgf_1.append(probs_1)
    hazard_pgf_2.append(probs_2)

# from simulation
sim_results = pd.read_csv('mtbp_pin_0.09.csv')

analysis = MTBPAnalysis()
# ext_prob_sim = analysis.extinction_probability(sim_results,'new_infections_both')
ext_prob_sim = analysis.extinction_probability_mt(sim_results)

# hazard function from simulation
hazard_sim = analysis.hazard_function_in_communities(sim_results)

# Visualisation
# prepare data for plotting
plotting_hazard = [{'gens_analyt': t_list, 'probs_analyt': hazard_pgf_both, 'label_analyt': 'type both (analytical) ',
                  'gens_num': hazard_sim.gens, 'probs_num': hazard_sim.probs_both, 'label_num': 'type both (simulation)'},
                 {'gens_analyt': t_list, 'probs_analyt': hazard_pgf_1, 'label_analyt': 'type 1 (analytical) ',
                  'gens_num': hazard_sim.gens, 'probs_num': hazard_sim.probs_1, 'label_num': 'type 1 (simulation)'},
                 {'gens_analyt': t_list, 'probs_analyt': hazard_pgf_2, 'label_analyt': 'type 2 (analytical) ',
                  'gens_num': hazard_sim.gens, 'probs_num': hazard_sim.probs_2, 'label_num': 'type 2 (simulation)'}
                 ]

plotting_extinction = [{'gens_analyt': t_list, 'probs_analyt': ext_prob_pgf_both, 'label_analyt': 'type both (analytical) ',
                  'gens_num': ext_prob_sim.gens, 'probs_num': ext_prob_sim.new_infections_both, 'label_num': 'type both (simulation)'},
                 {'gens_analyt': t_list, 'probs_analyt': ext_prob_pgf_1, 'label_analyt': 'type 1 (analytical) ',
                  'gens_num': ext_prob_sim.gens, 'probs_num': ext_prob_sim.new_infections_1, 'label_num': 'type 1 (simulation)'},
                 {'gens_analyt': t_list, 'probs_analyt': ext_prob_pgf_2, 'label_analyt': 'type 2 (analytical) ',
                  'gens_num': ext_prob_sim.gens, 'probs_num': ext_prob_sim.new_infections_2, 'label_num': 'type 2 (simulation)'}
                 ]

myplt = Plotting()
myplt.plot(plotting_extinction, title='hazard function', if_safe=False)

