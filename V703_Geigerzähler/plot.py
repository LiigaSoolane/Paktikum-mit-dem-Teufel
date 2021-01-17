import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from uncertainties import ufloat

Un,N=np.genfromtxt("Kennlinie.dat", unpack=True)
Ui, I=np.genfromtxt("Zaehlrohrstrom.dat", unpack=True)
I=I*10**(-6)

# a)
params, cov = np.polyfit(Un[5:34], N[5:34], deg=1, cov=True) #370-640 V
errors = np.sqrt(np.diag(cov))
for name, value, error in zip('ab', params, errors):
    print(f'{name} = {value:.3f} ± {error:.3f}')
fit = params[0]*Un[5:34]+params[1]

plt.errorbar(Un, N, yerr=N**(1/2), fmt='ro', label='Messwerte')
plt.plot(Un[5:34], fit, label='lineare Regression')
plt.xlabel(r'U [V]')
plt.ylabel(r'N [Imp]')
plt.savefig('plot.pdf')

# b) entfällt

# c)
#Totzeit anhand der Zwei-Quellen-Methode
N1 = 96401
N2 = 76518
N12 = 158479
T=(N1+N2-N12)/(2*N1*N2*120)
print(f'T=',T)

#Totzeit anhand des Oszilloskops

# d)
e0 = 1.602176634*10**(-19)
Z = np.zeros(8)
for i in range(7):
    Z[i] = I[i]/(e0*N[3+5*i])
print(f'Z = ', Z)


