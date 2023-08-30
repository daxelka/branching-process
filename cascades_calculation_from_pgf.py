import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft
plt.style.use(['ggplot', 'seaborn-talk'])

def Q1(s1, s2, lin = 8, lout=2, p=0.08):
    return np.exp(lin * p * (s1 - 1)) * np.exp(lout * p * (s2 - 1))

def Q2(s1, s2, lin = 8, lout=2, p=0.08):
    return np.exp(lin * p * (s2 - 1)) * np.exp(lout * p * (s1 - 1))

def G0(s1,s2,c1,c2):
    return s1 * c1

# Q1 = lambda s1: np.exp(2*(s1-1))
# G0 = lambda x,y: x**1*y**1
N = 400
n = np.arange(N)
c = np.exp(2*np.pi*1j*n/N)

tset = list(range(10))
c1 = c.copy()
c2 = c.copy()
s1 = np.ones_like(c)
s2 = np.ones_like(c)
# for t in range(max(tset)+1):
#     if t in tset:
#         pn = abs(np.fft.fft(G0(s1,s2,c1,c2))/N)
#         plt.plot(n, pn, label=fr"$t = {t}$")
#     s1 = Q1(s1*c1, s2*c2)
#     s2 = Q2(s1*c1, s2*c2)
#     c1 = c1
#     c2 = c2
#
# # plt.plot(n, pn, label=fr"$t = {t}$")
# plt.legend()
# plt.ylabel('Probability')
# plt.xlabel('Number of nodes')
# plt.show()

def itt_gen(s1, s2, c1, c2, lambda_param=0.9, n_max=100):
    # s1 and s2: dummy var of for each type
    # lambda_param: lambda parameter for poisson dist
    # n_max: number of generation
    for n in range(max(tset)+1):
        s1 = Q1(s1*c1, s2*c2)
        s2 = Q2(s1*c1, s2*c2)
    # apply initial conditions
    ans = G0(s1,s2,c1,c2)
    return ans

c1_ini = c.copy()
c2_ini = c.copy()
s1_ini = np.ones_like(c)
s2_ini = np.ones_like(c)

# pdf recovered from pgf
pdf = ifft(itt_gen(s1_ini, s2_ini, c1, c2))
print(pdf.real)

plt.plot(pdf)
plt.show()




# c1 = c.copy()
# s1 = np.ones_like(c)
# for t in range(max(tset)+1):
#     if t in tset:
#         pn = abs(np.fft.fft(G0(s1,c1))/N)
#         plt.plot(n, pn, label=fr"$t = {t}$")
#     s1 = Q(s1*c1)
#     c1 = c1
# plt.legend()
# plt.ylabel('Probability')
# plt.xlabel('Number of nodes')
# plt.show()