import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from scipy.optimize import curve_fit

#Values der Probe
d, l, b, R, lD = np.genfromtxt("hallsnostiges.txt", unpack=True)
d *= 10**(-3)
l *= 10**(-3)
b *= 10**(-3)
lD *= 10**(-2)

B1, I1, B2, I2 = np.genfromtxt('hallwerte.txt', unpack=True)
B1 *= 10**(-3)
B2 *= 10**(-3)

#fitte B1 und B2
parB1, covB1 = np.polyfit(I1, B1, deg=1, cov=True)
errors = np.sqrt(np.diag(covB1))

for name, value, error in zip('ab', parB1, errors):
    print(f'{name} = {value:.3f} +- {error:.3f}')

parB2, covB2 = np.polyfit(I2, B2, deg=1, cov=True)
errars = np.sqrt(np.diag(covB2))

for name, value, error in zip('ab', parB2, errars):
    print(f'{name} = {value:.3f} +- {error:.3f}')

x = np.linspace(0, 5, num = 1000)

plt.subplot(2, 1, 1)
plt.plot(I1, B1, 'rx', label="Messwerte")
plt.plot(x, parB1[0] * x + parB1[1], label="Lineare Regression")
plt.xlabel('Stromstärke I [A]')
plt.ylabel('Magnetische Flussdichte B [T]')
plt.title('Aufsteigender Strom')

plt.subplot(2, 1, 2)
plt.plot(I2, B2, 'rx', label="Messwerte")
plt.plot(x, parB2[0] * x + parB2[1], label= "Lineare Regression")
plt.xlabel('Stromstärke I [A]')
plt.ylabel('Magnetische Flussdichte B [T]')
plt.title('Abfallender Strom')

plt.legend()
plt.tight_layout()
plt.savefig('mag.pdf')

plt.clf()

linreg1 = unp.uarray([parB1[0], parB1[1]], [errors[0], errors[1]])
linreg2 = unp.uarray([parB2[0], parB2[1]], [errars[0], errars[1]])

linreg = 1/2 * (linreg1 + linreg2)
print(linreg)

# Uh = 1/2 (Uges+ - Uges-)

# Ih konst
Im, Uh, Im2, Uh2 = np.genfromtxt("hallspannung.txt", unpack=True)
Uh2 *= 10**(-3)
Uh *= 10**(-3)
Ih = 2
hallPk = 1/2 * (Uh - Uh2)


# Im konst
Ihall2, Uhall2, Ihall, Uhall  = np.genfromtxt("hallmagconst.txt", unpack=True)
Uhall2 *= 10**(-3)
Uhall *= 10**(-3)
Imag = 2
hallSk = 1/2 * (Uhall - Uhall2)
hallSk[1] += 0.000001
print(hallSk, hallPk)
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
e = -1.601* 10**(-19)
#Elektronenmasse in kg
me = 9.108* 10**(-31)
#Masse eines Kupferatoms in u
mCu = 63.4
#u Umrechnung in kg
u = 1.661* 10**(-27)
#Dichte Kupfer in kg/m^3
rhoCu = 8.92* 10**(3)
#Plancksches Wirkungsquantum in Js
h = 6.625* 10**(-34)
#Stormdingsda in A/mm^2
j = 1

# 1. Berechne Ladungsträgerdichte n
Ikonst = 2 # konstanter Strom hat immer 2 A

nSk_ = - (linreg[0] * Ihall + linreg[1])*Ikonst / (hallSk *e * d)
nPk_ = - (linreg[0]* Ikonst + linreg[1]) * Ihall / (hallPk * e * d)
nPk_[0] += (linreg[0] * Ikonst + linreg[1]) * 0.0001 / (hallPk[0] * e *d)
#print(nSk_, nPk_)
#nSk und nPk sollen die mean values sein, wie mache ich das?
#Nsk = np.mean(unp.nominal_values(nSk_))
#Npk = np.mean(unp.nominal_values(nPk_))
#Nskerr = np.mean(unp.std_devs(nSk_))
#Npkerr = np.mean(unp.std_devs(nPk_))
#
#nSk = unp.uarray(Nsk, Nskerr)
#nPk = unp.uarray(Npk, Npkerr)

nSk = np.sum(nSk_)/11
nPk = np.sum(nPk_)/11
print(nSk)
print(nPk)
# 2. Berechne Ladungsträger pro Atom z
m = mCu * u
a = rhoCu / m # Atomdichte in Kupfer

zSk = nSk / a
zPk = nPk / a

# 3. Berechne Driftgeschwindigkeit v_d
j *= 1/10**(-6)

v_d_Sk = j/(nSk * e)
v_d_Pk = j/(nPk * e)

# 4. Berechne mittlere Fluggeschwindigkeit tau
Q = np.pi * (d/2)**2 # Querschnitt des Leiters
tauSk = (2 * me * lD)/(R * e**2 * nSk * Q)
tauPk = (2 * me * lD)/(R * e**2 * nPk * Q)
# 5. Berechne die Beweglichkeit mu
muSk = - e * tauSk / (2*me)
muPk = -e * tauPk / (2*me)
# 6. Berechne die Geschwindigkeit v
Ef_Sk = (h**2 / (2*me)) * (3*nSk/(8*np.pi))**(1/3) # Fermi-Energie
Ef_Pk = (h**2 / (2*me)) * (3*nPk/(8*np.pi))**(1/3)

vSk = (2 * Ef_Sk / me)**(1/2)
vPk = (2 * Ef_Pk / me)**(1/2)
# 7. Berechne die mittlere freie Weglänge lf
lf_Sk = tauSk * vSk
lf_Pk = tauPk * vPk

# DONE!
# Jetzt noch ausgeben lassen

print("Preparate para un monton de valores")

print(f'n_Sk = {nSk}; n_Pk = {nPk}')
print(f'z_Sk = {zSk}; z_Pk = {zPk}')
print(f'vd_Sk = {v_d_Sk}; vd_Pk = {v_d_Pk}')
print(f'tau_Sk = {tauSk}; tau_Pk = {tauPk}')
print(f'mu_Sk = {muSk}; mu_Pk = {muPk}')
print(f'v_Sk = {vSk}; v_Pk = {vPk}')
print(f'l_Sk = {lf_Sk}; l_Pk = {lf_Pk}')

print("Listo!")