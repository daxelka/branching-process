import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.style.use(['seaborn-talk'])


def g1(s1, s2, c1, c2, lin=8, lout=2, p=0.08):
    return np.exp(lin * p * (s1*c1 - 1)) * np.exp(lout * p * (s2*c2 - 1))


def g2(s1, s2, c1, c2, lin=8, lout=2, p=0.08):
    return np.exp(lin * p * (s2*c2 - 1)) * np.exp(lout * p * (s1*c1 - 1))


def gc(c):
    return c

def G0(s1,s2, c1, c2):
    return s1*c1


def G_N_t(s1, s2, c1, c2, lin=8, lout=2, p=0.08, n_max=300):
    # s1 and s2: dummy var of for each type
    # lambda_param: lambda parameter for poisson dist
    # n_max: number of generation
    for n in range(n_max):
        s1 = g1(s1, s2, c1, c2, lin, lout, p)
        s2 = g2(s1, s2, c1, c2, lin, lout, p)
        c1 = gc(c1)
        c2 = gc(c2)
    # apply initial conditions
    ans = G0(s1, s2, c1, c2)
    return ans


def get_cascades(G_N_t, lin=8, lout=2, p=0.08, N_gen_type1=100, N_gen_type2=100):
    """
    Args:
     - G_N_t (function): multivariate pgf
     - N_gen_type1 (int): # generations for type 1
     - N_gen_type2 (int): # generations for type 2

    Returns:
    - cascade_dist_multivar (array of (N_gen_type1, N_gen_type2) shape): probability of cascades of type1 and type2
                                                                        for generations up to N_gen_type1 and N_gen_type2
                                                                        respectively
    - cascade_dist_total  (array of (N_gen_type1, ) shape): probability of total cascades for N_gen_type1 == N_gen_type2
    """

    G_N_t = np.vectorize(G_N_t)
    n1 = np.arange(N_gen_type1)
    n2 = np.arange(N_gen_type2)
    c1 = np.exp(2*np.pi*1j*n1/N_gen_type1)
    c2 = np.exp(2*np.pi*1j*n2/N_gen_type2)
    ones = np.ones_like(c1)
    # s2 = np.ones_like(c2)

    # total cascade size distribution
    if N_gen_type1 == N_gen_type2:
        # total cascades
        cascade_dist_total = abs(np.fft.fft(G_N_t(ones, ones, c1, c1, lin, lout, p))/N_gen_type1)
        cascade_dist_total = cascade_dist_total / np.sum(cascade_dist_total)
        # type 1 cascades
        cascade_dist_type1 = abs(np.fft.fft(G_N_t(ones, ones, c1, ones, lin, lout, p)) / N_gen_type1)
        cascade_dist_type1 = cascade_dist_type1 / np.sum( cascade_dist_type1)
        # type 2 cascades
        cascade_dist_type2 = abs(np.fft.fft(G_N_t(ones, ones, ones, c2, lin, lout, p)) / N_gen_type1)
        cascade_dist_type2 = cascade_dist_type2 / np.sum(cascade_dist_type2)
    else:
        raise ValueError('N_gen_type1 should be equal to N_gen_type2 for calculating total cascades')

    return pd.DataFrame({'cascades': list(range(1, N_gen_type1)),
                         'prob_both': list(cascade_dist_total[1:]),
                         'prob_1': list(cascade_dist_type1[1:]),
                         'prob_2': list(cascade_dist_type2[1:])})


def get_cascades_prob_multivariate(G_N_t, lin=8, lout=2, p=0.08, N_gen_type1=100, N_gen_type2=100):
    """
    Args:
     - G_N_t (function): multivariate pgf
     - N_gen_type1 (int): # generations for type 1
     - N_gen_type2 (int): # generations for type 2

    Returns:
    - cascade_dist_multivar (array of (N_gen_type1, N_gen_type2) shape): probability of cascades of type1 and type2
                                                                        for generations up to N_gen_type1 and N_gen_type2
                                                                        respectively
    """

    G_N_t = np.vectorize(G_N_t)
    n1 = np.arange(N_gen_type1)
    n2 = np.arange(N_gen_type2)
    c1 = np.exp(2*np.pi*1j*n1/N_gen_type1)
    c2 = np.exp(2*np.pi*1j*n2/N_gen_type2)
    s1 = np.ones_like(c1)
    s2 = np.ones_like(c2)

    C1, C2 = np.meshgrid(c1, c2)
    cascade_dist_multivar = abs(np.fft.fft2(G_N_t(s1, s2, C1, C2, lin, lout, p))/(N_gen_type1*N_gen_type2))
    cascade_dist_multivar = cascade_dist_multivar / np.sum(cascade_dist_multivar)

    return cascade_dist_multivar

# write to file
# cascades.to_csv('data/cascades/from_pgfs/p_' + str(probability) + '_cascades_pgf.csv')
# print('saved to a file')

lin = 8
lout = 6
pin = 0.08

cascades = get_cascades(G_N_t, lin, lout, pin, 100, 100)

cascade_dist_multivar = get_cascades_prob_multivariate(G_N_t, lin, lout, pin, 100, 100)


# plotting
plt.plot(cascades['cascades'], cascades['prob_both'], label='both')
plt.plot(cascades['cascades'], cascades['prob_1'], label='both')
# plt.yscale('log')
# plt.xscale('log')
# plt.ylim([1e-06,1])
plt.xlim([0,6])
plt.ylim([0,6])
plt.xlabel('cascade size')
plt.ylabel('probability')
plt.legend()
plt.show()

# plotting multivariate
cmap = 'viridis'
plt.imshow(cascade_dist_multivar, origin='lower', cmap=cmap)
plt.colorbar(label='Probability')
# plt.contourf(cascade_dist_multivar, cmap='viridis', levels=6)
# plt.colorbar(label='Probability')
plt.ylabel('cascades size, community 2')
plt.xlabel('cascades size, community 1')
plt.xlim([0,6])
plt.ylim([0,6])
plt.show()


#  # branching process parameters
# lambda_in, lambda_out = 8, 2
# probabilities_infection = [0.06, 0.08, 0.09, 0.095]
# cascades_list = []
#
# for p in probabilities_infection:
#     _, cascade_dist_total = get_cascades(G_N_t, lin=8, lout=2, p=p, N_gen_type1=1000, N_gen_type2=1000)
#     # somehow need to delete zero cascades
#     cascades_list.append(cascade_dist_total)
#     plt.plot(cascade_dist_total['cascades'], cascade_dist_total['prob'], label='p: '+str(round(p/0.1, 2))+' p*')
#
# plt.yscale('log')
# plt.xscale('log')
# # plt.ylim([1e-06,1])
# plt.xlabel('cascade size')
# plt.ylabel('probability')
# plt.legend()
# plt.show()

# plt.imshow(cascade_dist_multivar, origin='lower')
# plt.colorbar(label='Probability')
# plt.ylabel(r'$n_1$')
# plt.xlabel(r'$n_2$')
# plt.show()

