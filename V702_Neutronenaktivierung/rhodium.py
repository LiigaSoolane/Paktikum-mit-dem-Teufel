import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
from scipy.optimize import curve_fit

# leer los datos y crear variables y funciones
tr, Nr = np.genfromtxt('Rhodium.dat', unpack=True)
Nu = np.array([129.0, 143.0, 144.0, 136.0, 126.0, 158.0])

Nu_temp = Nu/20 # con temporiz corrigida: 15sec cositas
Nu_cor = unp.uarray(Nu_temp, np.sqrt(Nu_temp))
Nu_mean = np.mean(Nu_cor)

Nr_cor = unp.uarray(Nr, np.sqrt(Nr))
Nr_diff = Nr_cor - Nu_mean
Nr_diff_log = unp.log(Nr_diff)

def lnn(x, lamb, const):
    return - lamb * x + const

# curvefitear la parte mas lenta
paramsl, covl = curve_fit(lnn, tr[27:], unp.nominal_values(Nr_diff_log[27:]), sigma=unp.std_devs(Nr_diff_log[27:]))
errorsl = np.sqrt(np.diag(covl))
lambl = unp.uarray(paramsl[0], errorsl[0])
constl = unp.uarray(paramsl[1], errorsl[1])

# numero de impulsos de la parte mas rapida
Nr_k = Nr_diff - unp.exp(lnn(tr, lambl, constl))
Nr_k_log = unp.log(Nr_k[:8])

# curvefitear la parte mas rapida
paramsk, covk = curve_fit(lnn, tr[:8], unp.nominal_values(Nr_k_log), sigma=unp.std_devs(Nr_k_log))
errorsk = np.sqrt(np.diag(covk))
lambk = unp.uarray(paramsk[0], errorsk[0])
constk = unp.uarray(paramsk[1], errorsk[1])

# tiempos de medio
Tl = unp.uarray(np.log(2)/paramsl[0], np.sqrt(np.log(2)**2 * errorsl[0]**2))
Tk = unp.uarray(np.log(2)/paramsk[0], np.sqrt(np.log(2)**2 * errorsk[0]**2))

# muestrame los valores
print(f'lambda slow = {lambl}, lambda fast = {lambk}')
print(f'const slow = {constl}, const fast = {constk}')
print(f'Halbwertszeit slow = {Tl}, Halbwertszeit fast = {Tk}')

# hacer el grafico
tstern = 400
tmax = 120
x = np.linspace(0, 630)
xl = np.linspace(tstern, 630)
xk = np.linspace(0, tmax)

plt.errorbar(tr, unp.nominal_values(Nr_diff_log), yerr = unp.std_devs(Nr_diff_log), fmt='b.', label="Messwerte mit Fehlern")
plt.plot(xl, unp.nominal_values(lnn(xl, *paramsl)), 'r--', label='Fit für langsamen Zerfall')
plt.plot(xk, unp.nominal_values(lnn(xk, *paramsk)), 'y--', label='Fit für schnellen Zerfall')
#plt.plot(tr, unp.nominal_values(lnn(tr, *paramsl)) + unp.nominal_values(lnn(tr, *paramsk)), 'g-', label='ich will schlafen')
plt.yscale('log')
plt.xlabel('t [s]')
plt.ylabel('ln(N)')
plt.legend()
plt.tight_layout()
plt.savefig('rhodium.pdf')