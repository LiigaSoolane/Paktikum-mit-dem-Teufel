import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

#Werte für die ???? What am I even doing?
theta, N = np.genfromtxt('EmissionCu.txt', unpack=True)

#plot values that I don't know what they're for
plt.plot(theta, N, 'r.', label='Messwerte')
plt.plot(theta, N, 'b-')
plt.xlabel('Einfallswinkel in °')
plt.ylabel('Anzahl Impulse pro s')
plt.legend()
plt.savefig('savefig.pdf')