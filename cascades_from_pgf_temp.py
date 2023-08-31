import numpy as np
import matplotlib.pyplot as plt
plt.style.use(['seaborn-talk'])


def g1(s1, s2, c1, c2, lin=8, lout=2, p=0.08):
    return np.exp(lin * p * (s1*c1 - 1)) * np.exp(lout * p * (s2*c2 - 1))


def g2(s1, s2, c1, c2, lin=8, lout=2, p=0.08):
    return np.exp(lin * p * (s2*c2 - 1)) * np.exp(lout * p * (s1*c1 - 1))


def gc(c):
    return c

def G0(s1,s2, c1, c2):
    return s1*c1


def G_N_t(s1, s2, c1,c2, lin = 8, lout=2, p=0.08, n_max=10):
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


def get_cascades(G_N_t, N_gen_type1, N_gen_type2):
    G_N_t = np.vectorize(G_N_t)
    n1 = np.arange(N_gen_type1)
    n2 = np.arange(N_gen_type2)
    c1 = np.exp(2*np.pi*1j*n1/N_gen_type1)
    c2 = np.exp(2*np.pi*1j*n2/N_gen_type2)
    s1 = np.ones_like(c1)
    s2 = np.ones_like(c2)

    # cascade distribution multivariete
    C1, C2 = np.meshgrid(c1, c2)
    cascade_dist_multivar = abs(np.fft.fft2(G_N_t(s1, s2, C1, C2))/(N_gen_type1*N_gen_type2))
    cascade_dist_multivar = cascade_dist_multivar / np.sum(cascade_dist_multivar)

    # total cascade size distribution
    if N_gen_type1 == N_gen_type2:
        cascade_dist_total = abs(np.fft.fft(G_N_t(s1, s2, c1, c1))/N_gen_type1)
        cascade_dist_total = cascade_dist_total / np.sum(cascade_dist_total)
    else:
        raise ValueError('N_gen_type1 should be equal to N_gen_type2 for calculating total cascades')

    return cascade_dist_multivar, cascade_dist_total


cascade_dist_multivar, cascade_dist_total = get_cascades(G_N_t, 50, 50)

print(cascade_dist_total.shape)
print(cascade_dist_multivar.shape)


plt.imshow(cascade_dist_multivar, origin='lower')
plt.colorbar(label='Probability')
plt.ylabel(r'$n_1$')
plt.xlabel(r'$n_2$')
plt.show()

plt.plot(cascade_dist_total)
plt.show()


# g = lambda x: np.sum([x**n/6 for n in range(1,7)])
# G = lambda x1,x2: g(x1)*g(x2)*g(x1*x2)
# G = np.vectorize(G)
# N1 = 13
# N2 = 13
# n1 = np.arange(N1)
# n2 = np.arange(N2)
# c1 = np.exp(2*np.pi*1j*n1/N1)
# c2 = np.exp(2*np.pi*1j*n2/N2)
# # c2 = np.ones_like(c1)
# C1, C2 = np.meshgrid(c1, c2) #get 2d array versions
# pn1n2 = abs(np.fft.fft2(G(C1,C2))/(N1*N2))
#
#
# plt.imshow(pn1n2, origin='lower')
# plt.colorbar(label='Probability')
# plt.ylabel(r'$n_1$')
# plt.xlabel(r'$n_2$')
# plt.show()