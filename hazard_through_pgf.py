import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

lin, lout = 8, 2
pin = 0.09

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

# probability of extinction q(t)=P[N(t)=0]=G_N(t)(0)
t_list = list(range(11))
probs_com1 = []
probs_com2 = []
for t in t_list:
    probs_com1.append(G_N_t(0, 1, t, lin, lout, pin))
    probs_com2.append(G_N_t(1, 0, t, lin, lout, pin))
    # print('com 1: ', G_N_t(0, 1, t, lin, lout, pin), 'com 2: ', G_N_t(1, 0, t, lin, lout, pin))


df = pd.read_csv('full_extin_dist_v4.csv')


plt.plot(t_list, probs_com1)
plt.plot(t_list, probs_com2)
plt.scatter(df.t[df.pin ==pin], df.q_1[df.pin ==pin])
plt.scatter(df.t[df.pin ==pin], df.q_2[df.pin ==pin])
plt.xlabel('generation')
plt.ylabel('probability')
plt.title('probability of extinction')
plt.show()
