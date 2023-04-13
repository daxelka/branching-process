import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mtbp_analysis_class import MTBPAnalysis
from pgf_analysis_class import PGFAnalysis
import time

lin, lout = 8, 2
pin = 0.06

# plotting

def blue_orange_green_set():
    return '#1f77b4',  '#ff7f0e', '#2ca02c'

def purple_red_teal_set():
    # return '#9467bd', '#8c564b', '#17becf'
    return '#9467bd', '#af3636', '#17becf', '#FFDB58', '#A9A9A9'

def CB91_color_set():
    return '#2CBDFE','#47DBCD','#F3A0F2','#9D2EC5','#661D98','#F5B14C'


def plotting(curves_list, title='my title', xlim=[0,10], if_safe=False):
    my_colors = purple_red_teal_set()
    my_linewidth = 0.8
    my_fontsize = 14
    my_mark_size = 14

    for n, curves_dict in enumerate(curves_list):
        plt.plot(curves_dict['gens_analyt'], curves_dict['probs_analyt'],
                 c=my_colors[n], linewidth=my_linewidth, label=curves_dict['label_analyt'])
        plt.scatter(curves_dict['gens_num'], curves_dict['probs_num'],
                 c=my_colors[n], s=my_mark_size, label=curves_dict['label_num'])

    plt.xlabel('generation', fontsize=my_fontsize)
    plt.ylabel('probability', fontsize=my_fontsize)
    plt.title(title, fontsize=my_fontsize)
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.legend(frameon=False)
    plt.xlim(xlim)
    if if_safe:
        plt.savefig('img/'+ title +'.png', dpi=300)
    else:
        plt.show()
    print('plotted')

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
sim_results = pd.read_csv('mtbp_pin_0.06.csv')

analysis = MTBPAnalysis()
# ext_prob_sim = analysis.extinction_probability(sim_results,'new_infections_both')
ext_prob_sim = analysis.extinction_probability_mt(sim_results)

# hazard function from simulation
hazard_sim = analysis.hazard_function_in_communities(sim_results)

# prepare data for plotting
plotting_data = [{'gens_analyt': t_list, 'probs_analyt': hazard_pgf_both, 'label_analyt': 'type both (analytical) ',
                  'gens_num': hazard_sim.gens, 'probs_num': hazard_sim.probs_both, 'label_num': 'type both (simulation)'},
                 {'gens_analyt': t_list, 'probs_analyt': hazard_pgf_1, 'label_analyt': 'type 1 (analytical) ',
                  'gens_num': hazard_sim.gens, 'probs_num': hazard_sim.probs_1, 'label_num': 'type 1 (simulation)'},
                 {'gens_analyt': t_list, 'probs_analyt': hazard_pgf_2, 'label_analyt': 'type 2 (analytical) ',
                  'gens_num': hazard_sim.gens, 'probs_num': hazard_sim.probs_2, 'label_num': 'type 2 (simulation)'}
                 ]

plotting(plotting_data, title='hazard function', if_safe=False)

# plt.plot(t_list, hazard_pgf_both, c=my_colors[0])
# plt.scatter(hazard_sim_v1.gens, hazard_sim_v1.probs_both, c=my_colors[0])
# plt.plot(t_list, hazard_pgf_1, c=my_colors[1])
# plt.scatter(hazard_sim_v1.gens, hazard_sim_v1.probs_1, c=my_colors[1])
# plt.plot(t_list, hazard_pgf_2, c=my_colors[2])
# plt.scatter(hazard_sim_v1.gens, hazard_sim_v1.probs_2, c=my_colors[2])
# plt.xlim([0,20])
# plt.show()
