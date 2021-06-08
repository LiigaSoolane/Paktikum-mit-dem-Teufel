import numpy as np
import matplotlib.pyplot as plt

d = 201.4 * 10**(-12)

# Überprüfen der Bragg-Bedingung

theta_bragg, N_bragg = np.genfromtxt('Bragg.dat', unpack=True)

N_bragg_max = np.amax(N_bragg)
mask_ismax = N_bragg == N_bragg_max
theta_max = theta_bragg[mask_ismax]


print(f'Maximum: {N_bragg_max:.3f}, bei {theta_max}.')

# diferencia procentual al valor teoretico
delta = 100 * (28 - theta_max)/28

print(f'prozentuale Abweichung: {delta}')

plt.plot(theta_bragg, N_bragg, 'b.', label = 'Messwerte')
plt.plot(theta_bragg, N_bragg, 'b--')
plt.plot(theta_max, N_bragg_max, 'r*', label = 'Maximum')
plt.xlabel('Theta [°]')
plt.ylabel('Number of Impulses')
plt.legend()
plt.savefig('bra.pdf')