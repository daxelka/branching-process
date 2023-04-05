import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mtbp_analysis_class import MTBPAnalysis

lin, lout = 8, 2
pin = 0.06

def g_X_1(s1, s2, lin, lout, pin):
    ans = np.exp(lin * pin * (s1 - 1)) * np.exp(lout * pin * (s2 - 1))
    return ans

def g_X_2(s1, s2, lin, lout, pin):
    ans = np.exp(lin * pin * (s2 - 1)) * np.exp(lout * pin * (s1 - 1))
    return ans

def gic_1(s1, s2):
    return s1


def G_N_t(s1, s2, t, lin, lout, pin):
    if t > 0:  # have to iterate the pgf at least one
        for t_i in range(t, 0, -1):  # iterate the function t times
            new_s1 = g_X_1(s1, s2, lin, lout, pin)
            new_s2 = g_X_2(s1, s2, lin, lout, pin)
            s1 = new_s1
            s2 = new_s2
    ans = gic_1(s1, s2)  # apply initial conditions
    return ans

def hazard_pgf(t, lin, lout, pin):
    if t > 0:
        h = (G_N_t(0,0,t,lin,lout,pin) - G_N_t(0,0,t-1,lin,lout,pin))/(1 - G_N_t(0,0,t-1,lin,lout,pin))
        return h
    else:
        return 0

# def G_H_t(s1, s2, t, lin=8, lout=2, pin=0.05):
#     if t >= 0:
#         G_Y_t_num = G_Y_t(g_X_1(s1, s2, lin, lout, pin),
#                           g_X_2(s1, s2, lin, lout, pin),
#                           t-1, lin, lout, pin)
#         ans = G_Y_t(s1, s2, t, lin, lout, pin) - G_Y_t_num
#     else:
#         ans = 0
#     return ans


# probability of extinction q(t)=P[N(t)=0]
# from pgf q(t) = G_N(t)(0)
t_list = list(range(20))
ext_prob_pgf = []
for t in t_list:
    ext_prob_pgf.append(G_N_t(0, 0, t, lin, lout, pin))

# from simulation
sim_results = pd.read_csv('mtbp_pin_0.06.csv')
analysis = MTBPAnalysis()
ext_prob_sim = analysis.extinction_probability(sim_results)

# plotting
plt.plot(t_list, ext_prob_pgf)
plt.scatter(ext_prob_sim.gens, ext_prob_sim.probs)
plt.show()


