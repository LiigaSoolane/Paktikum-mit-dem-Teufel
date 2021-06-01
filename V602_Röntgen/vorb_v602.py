import numpy as np 

Z = np.array([30, 32, 35, 37, 38, 40])
E_K = np.array([9.65, 9.855, 11.877, 13.335, 14.098, 15.691])

E_K *= 10**3

c = 2.998 * 10**(8)
h = 4.1357 * 10**(-15)
d = 201.4
n = 1
R = 13.6
alpha = 7.297 * 10**(-3)

d *= 10**(-12)

def theta(x):
    return np.arcsin(n * c * h/(2 * x * d))

def sigma(x, y):
    return x - np.sqrt(y/R - alpha**2 * x**4 / 4)

for name, z, e in zip(['Zn', 'Ge', 'Br', 'Rb', 'Sr', 'Zr'], Z, E_K):
    print(f'{name}: theta = {theta(e) * 180 / np.pi}, sigma = {sigma(z, e)}')
