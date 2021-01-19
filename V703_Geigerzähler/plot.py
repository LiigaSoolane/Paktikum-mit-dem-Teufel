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
plt.savefig('build/plot.pdf')

np.savetxt("NgegenU.txt", np.column_stack([Un, N]), fmt = "%10.4f", delimiter = " & ", newline = " \\\ ", header = " U N")

# b) entfällt

# c)
#Totzeit anhand der Zwei-Quellen-Methode
N1 = ufloat(96401, 96401**(1/2))
N2 = ufloat(76518, 76518**(1/2))
N12 = ufloat(158479, 158479**(1/2))
T=(N1+N2-N12)/(2*N1*N2*120)
print(f'T=',T)

#Totzeit anhand des Oszilloskops

# d)
e0 = 1.602176634*10**(-19)
Z = np.zeros(8)
for i in range(8):
    Z[i] = I[i]*60/(e0*N[3+5*i])
print(f'Z = ', Z)
Er = np.zeros(8)
for i in range(8):
    Er[i] = Z[i]*(1/(N[3+5*i])+(0.05*10**(-6)/I[i])**2)**(1/2)
print(f'Error = ', Er)
np.savetxt("Zerr.txt", np.column_stack([Z, Er]), fmt = "%10.4f", delimiter = " & ", newline = " \\\ ", header = "Z Error")

plt.errorbar(Ui, Z, yerr=Er, fmt='ro', label='Messwerte')
plt.xlabel(r'U [V]')
plt.ylabel(r'Zahl der freigesetzen Ladungen pro Teilchen')
plt.savefig('build/plot1.pdf')

np.savetxt("ZgegenU.txt", np.column_stack([Ui, Z]), fmt = "%10.4f", delimiter = " & ", newline = " \\\ ", header = "U Z")

