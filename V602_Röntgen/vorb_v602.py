import numpy as np 
import uncertainties.unumpy as unp

Z = np.array([30, 31, 35, 37, 38, 40])
E_K = unp.uarray([9668.55, 10377.76, 13483.86, 15207.74, 16115.26, 18008.15], [15, 16, 19, 22, 23, 26])

c = 2.998 * 10**(8)
h = 4.1357 * 10**(-15)
d = 201.4
n = 1
R = 13.6
alpha = 7.297 * 10**(-3)

d *= 10**(-12)

def theta(x):
    return unp.arcsin(n * c * h/(2 * x * d))

def sigma(x, y):
    return x - unp.sqrt(y/R - alpha**2 * x**4 / 4)

for name, z, e in zip(['Zn', 'Ge', 'Br', 'Rb', 'Sr', 'Zr'], Z, E_K):
    print(f'{name}: theta = {theta(e) * 180 / np.pi}, sigma = {sigma(z, e)}')
