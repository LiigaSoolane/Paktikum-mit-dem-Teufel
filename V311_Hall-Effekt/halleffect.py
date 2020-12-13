import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from scipy.optimize import curve_fit

#Values der Probe
d, l, b, R, lD = np.genfromtxt("hallsnostiges.txt", unpack=True)
d *= 10**3
l *= 10**3
b *= 10**3
lD *= 10**2

B1, I1, B2, I2 = np.genfromtxt('hallwerte.txt', unpack=True)
B1 *= 10**3
B2 *= 10**3

#fitte B1 und B2
parB1, covB1 = np.polyfit(I1, B1, deg=1, cov=True)
errors = np.sqrt(np.diag(covB1))

for name, value, error in zip('ab', parB1, errors):
    print(f'{name} = {value:.3f} +- {error:.3f}')

x = np.linspace(0, 5, num = 1000)


plt.plot(I1, B1, 'rx', label="Messwerte")
plt.plot(x, parB1[0] * x + parB1[1], label="Lineare Regression")
plt.legend(loc="best")
plt.xlabel('Stromstärke I [A]')
plt.ylabel('Magnetische Flussdichte B [T]')
plt.savefig('unswitched.pdf')

parB2, covB2 = np.polyfit(I2, B2, deg=1, cov=True)
errars = np.sqrt(np.diag(covB2))

for name, value, error in zip('ab', parB2, errars):
    print(f'{name} = {value:.3f} +- {error:.3f}')

plt.clf()

plt.plot(I2, B2, 'rx', label="Messwerte")
plt.plot(x, poroms[0] * x + poroms[1], label= "Lineare Regression")
plt.legend(loc="best")
plt.xlabel('Stromstärke I [A]')
plt.ylabel('Magnetische Flussdichte B [T]')
plt.savefig('switched.pdf')

plt.clf()

linreg1 = unumpy.uarray([parB1[0], parB1[1]], [errors[0], errors[1]])
linreg2 = unumpy.uarray([parB2[0], parB2[1]], [errars[0], errars[1]])

linreg = 1/2 * (linreg1 + linreg2)

# Uh = 1/2 (Uges+ - Uges-)

# Ih konst
Im, Uh, Im2, Uh2 = np.genfromtxt("hallspannung.txt", unpack=True)
Ih = 2
hallPk = 1/2 * (Uh - Uh2)


# Im konst
Ihall2, Uhall2, Ihall, Uhall  = np.genfromtxt("hallmagconst", unpack=True)
Imag = 2
hallSk = 1/2 * (Uhall - Uhall2)

paramPk, covPk = np.polyfit(Im, hallPk, deg=1, cov=True)
paramSk, covSk = np.polyfit(Ihall, hallSk, deg=1, cov=True)

errPk = np.sqrt(np.diag(covPk))
errSk = np.sqrt(np.diag(covSk))

for name, par, err in zip('ab', paramPk, errPk):
    print(f'{name} = {value:.3f} +- {err:.3f}')

for name, par, err in zip('ab', paramSk, errSk):
    print(f'{name} = {value:.3f} +- {err:.3f}')

# Naturkonstanten
#Elementarladung e0 in C
e0 = -1.601* 10**(-19)
#Elektronenmasse in kg
m0 = 9.108* 10**(-31)
#Masse eines Kupferatoms in u
mCu = 63.4
#u Umrechnung in kg
u = 1.661* 10**(-27)
#Dichte Kupfer in g/cm^3
rhoCu = 8.92
#Plancksches Wirkungsquantum in Js
h = 6.625* 10**(-34)
#Stormdingsda in A/mm^2
j = 1

# 1. Berechne Ladungsträgerdichte n
Ikonst = 2 # konstanter Strom hat immer 2 A

nSk = - (linreg[0] * I1 + linreg[1])*Ikonst / (hallSk * e0 * d)
nPk = - (linreg[0]* Ikonst + linreg[1]) * I1 / (hallPk * e0 * d)

# 2. Berechne Ladungsträger pro Atom z
m = mCu * u
a = rhoCu / m # Atomdichte in Kupfer

zSk = nSk / a
zPk = nPk / a

# 3. Berechne Driftgeschwindigkeit v_d
j *= 1/10**(-6)

v_d_Sk = j / (nSk*e0)
v_d_Pk = j / (nPk*e0)

# 4. Berechne mittlere Fluggeschwindigkeit tau
Q = np.pi * (d/2)**2 # Querschnitt des Leiters
 tauSk = (2 * m0 * lD)/(R * e0**2 * nSk * Q)
 tauPk = (2 * m0 * lD)/(R * e0**2 * nPk * Q)

 # 5. Berechne die Beweglichkeit mu
 muSk = - e0 * tauSk / (2*m0)
 muPk = -e0 * tauPk / (2*m0)

 # 6. Berechne die Geschwindigkeit v
 Ef_Sk = (h**2 / (2*m0)) * (3*nSk/(8*np.pi))**(1/3) # Fermi-Energie
 Ef_Pk = (h**2 / (2*m0)) * (3*nPk/(8*np.pi))**(1/3)
 vSk = np.sqrt(2*Ef_Sk/m0)
 vPk = np.sqrt(2*Ef_Pk/m0)

 # 7. Berechne die mittlere freie Weglänge lf
 lf_Sk = tauSk * vSk
 lf_Pk = tauPk * vPk

# DONE!