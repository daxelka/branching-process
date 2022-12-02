import numpy as np
from scipy.fft import fft, ifft
import matplotlib.pyplot as plt

x = np.array([1.0, 2.0, 1.0, -1.0, 1.5])
y = fft(x)
print(y)
yinv = ifft(y)
print(yinv)
z = np.exp(50.0 * 1.j * 2.0*np.pi*x) + 0.5*np.exp(-80.0 * 1.j * 2.0*np.pi*x)
print(z)

# that are M point apart
M = 100
z_i = np.array(range(M))/M
x_i = np.exp(2*np.pi*1.j*z_i)

print(x_i)
plt.plot(x_i.real, x_i.imag, '.')
plt.show()


def get_fx(s1, s2, lambda_param = 0.9):
    return np.exp(lambda_param * (s1 * s2 - 1))


def get_ft(s):
    return s


def get_initial_condition(s1, s2):
    return s1 * s2


def itt_gen(s1, s2, lambda_param=0.9, n_max=100):
    # s1 and s2: dummy var of for each type
    # lambda_param: lambda parameter for poisson dist
    # n_max: number of generation
    for n in range(n_max):
        s1 = get_fx(s1, s2, lambda_param)
        s2 = get_ft(s2)
    # apply initial conditions
    ans = np.conj(get_initial_condition(s1, s2))
    return ans


# pdf recovered from pgf
pdf = ifft(itt_gen(1, x_i))
print(pdf.real)

plt.plot(pdf)
plt.show()
