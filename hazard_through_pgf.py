import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

lin, lout = 8, 2
pin = 0.9

def g_X_1(s1, s2, lin, lout, pin):
    ans = np.exp(lin * pin * (s1 - 1)) * np.exp(lout * pin * (s2 - 1))
    return ans

def g_X_2(s1, s2, lin, lout, pin):
    ans = np.exp(lin * pin * (s2 - 1)) * np.exp(lout * pin * (s1 - 1))
    return ans

def gic_1(s1, s2):
    ans = s1 + s2 - 1
    return ans


def G_N_t(s1, s2, t, lin, lout, pin):
    if t > 0:  # have to iterate the pgf at least one
        for t_i in range(t, 0, -1):  # iterate the function t times
            new_s1 = g_X_1(s1, s2, lin, lout, pin)
            new_s2 = g_X_2(s1, s2, lin, lout, pin)
            s1 = new_s1
            s2 = new_s2
    ans = gic_1(s1, s2)  # apply initial conditions
    return ans


# pgf for for N(t) > 0
def G_Y_t(s1, s2, t, lin=8, lout=2, pin=0.05):
    if t >= 0:
        # find the prob of process ending at t
        G_N_t_num = G_N_t(0, 0, t, lin, lout, pin)

        # create pgf for G_Y_t by subtracting G_N_t_num and normalising it
        ans = (G_N_t(s1, s2, t, lin, lout, pin) - G_N_t_num) / (1 - G_N_t_num)
    else:
        ans = 0
    return ans

def G_H_t(s1, s2, t, lin=8, lout=2, pin=0.05):
    if t >= 0:
        G_Y_t_num = G_Y_t(g_X_1(s1, s2, lin, lout, pin),
                          g_X_2(s1, s2, lin, lout, pin),
                          t-1, lin, lout, pin)
        ans = G_Y_t(s1, s2, t, lin, lout, pin) - G_Y_t_num
    else:
        ans = 0
    return ans

t = list(range(11))
full_sur_dis_df = pd.DataFrame({'t': t})
full_sur_dis_df['hazard_1'] = full_sur_dis_df['t'].apply(lambda x: G_H_t(0, 1, x, lin, lout, pin))
full_sur_dis_df['hazard_2'] = full_sur_dis_df['t'].apply(lambda x: G_H_t(1, 0, x, lin, lout, pin))
full_sur_dis_df['hazard_b'] = full_sur_dis_df['t'].apply(lambda x: G_H_t(0, 0, x, lin, lout, pin))
full_sur_dis_df['q_1'] = full_sur_dis_df['t'].apply(lambda x: G_N_t(0, 1, x, lin, lout, pin))
full_sur_dis_df['q_2'] = full_sur_dis_df['t'].apply(lambda x: G_N_t(1, 0, x, lin, lout, pin))
full_sur_dis_df['q_b'] = full_sur_dis_df['t'].apply(lambda x: G_N_t(0, 0, x, lin, lout, pin))


plt.plot(full_sur_dis_df.t,full_sur_dis_df.q_b)
plt.show()