import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from uncertainties import ufloat

Un,N=np.genfromtxt("Kennlinie.dat", unpack=True)
N = N/60

# a)
params, cov = np.polyfit(Un[5:34], N[5:34], deg=1, cov=True) #370-640 V
errors = np.sqrt(np.diag(cov))
for name, value, error in zip('ab', params, errors):
    print(f'{name} = {value:.3f} ± {error:.3f}')
fit = params[0]*Un[5:34]+params[1]

plt.figure(0)
plt.errorbar(Un, N, yerr=N**(1/2), fmt='rx', label='Messwerte')
plt.plot(Un[5:34], fit, label='lineare Regression')
plt.xlabel(r'U [V]')
plt.ylabel(r'N [Imp/s]')
plt.savefig('build/plot.pdf')

np.savetxt("NgegenU.txt", np.column_stack([Un, N]), fmt = "%10.4f", delimiter = " & ", newline = " \\\ ", header = " U N")

# b) entfällt

# c)
#Totzeit anhand der Zwei-Quellen-Methode
N1 = ufloat(96401/120, (96401/120)**(1/2))
N2 = ufloat(76518/120, (76518/120)**(1/2))
N12 = ufloat(158479/120, (158479/120)**(1/2))
T=(N1+N2-N12)/(2*N1*N2)
print(f'T=',T)

#Totzeit anhand des Oszilloskops

# d)
Ui, I, N=np.genfromtxt("Zaehlrohrstrom.dat", unpack=True)
I=I*10**(-6)
N = N/60
e0 = 1.602176634*10**(-19)
Z = np.zeros(8)
for i in range(8):
    Z[i] = I[i]/(N[i]*e0)
print(f'Z = ', Z)
Er = np.zeros(8)
for i in range(8):
    Er[i] = Z[i]*(1/(N[i])+(0.05*10**(-6)/I[i])**2)**(1/2)
print(f'Error = ', Er)
np.savetxt("Zerr.txt", np.column_stack([Ui, Z, Er]), fmt = "%10.4f", delimiter = " & ", newline = " \\\ ", header = "Z Error")

plt.figure(1)
plt.errorbar(Ui, Z, yerr=Er, fmt='gx', label='Messwerte')
plt.xlabel(r'U [V]')
plt.ylabel(r'Zahl der freigesetzen Ladungen pro Teilchen')
plt.savefig('build/plot1.pdf')

np.savetxt("ZgegenU.txt", np.column_stack([Ui, Z]), fmt = "%10.4f", delimiter = " & ", newline = " \\\ ", header = "U Z")

