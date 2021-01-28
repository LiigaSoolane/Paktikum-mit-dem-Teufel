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

# parameters para encontrar a la version mas adecuada
a = 17
b = 14
tstern = 270
tmax = 150

# curvefitear la parte mas lenta
paramsl, covl = curve_fit(lnn, tr[a:], unp.nominal_values(Nr_diff_log[a:]), sigma=unp.std_devs(Nr_diff_log[a:]))
errorsl = np.sqrt(np.diag(covl))
lambl = unp.uarray(paramsl[0], errorsl[0])
constl = unp.uarray(paramsl[1], errorsl[1])

# numero de impulsos de la parte mas rapida
Nr_k = Nr_diff - unp.exp(lnn(tr, lambl, constl))
Nr_k_log = unp.log(Nr_k[:b])

# curvefitear la parte mas rapida
paramsk, covk = curve_fit(lnn, tr[:b], unp.nominal_values(Nr_k_log), sigma=unp.std_devs(Nr_k_log))
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
x = np.linspace(0, 630)
xl = np.linspace(tstern, 630)
xk = np.linspace(0, tmax)

plt.errorbar(tr, unp.nominal_values(Nr_diff_log), yerr = unp.std_devs(Nr_diff_log), fmt='b.', label="Messwerte mit Fehlern")
plt.plot(xl, unp.nominal_values(lnn(xl, *paramsl)), 'r--', label='Fit für langsamen Zerfall')
plt.plot(xk, unp.nominal_values(lnn(xk, *paramsk)), 'y--', label='Fit für schnellen Zerfall')
plt.plot(tr, np.log(unp.nominal_values(unp.exp(lnn(tr, lambl, constl)) + unp.exp(lnn(tr, lambk, constk)))), 'g-', label='ich will schlafen')
plt.yscale('log')
plt.xlabel('t [s]')
plt.ylabel('ln(N)')
plt.legend()
plt.tight_layout()
plt.savefig('rhodium.pdf')

plt.clf()
plt.errorbar(tr, unp.nominal_values(Nr_diff_log), yerr = unp.std_devs(Nr_diff_log), fmt='b.', label="Messwerte mit Fehlern")
plt.yscale('log')
plt.xlabel('t [s]')
plt.ylabel('ln(N)')
plt.legend()
plt.tight_layout()
plt.savefig('werte.pdf')

# hacer mierdos
np.savetxt(
    'build/tab/Rhodium_mess.txt',
    np.column_stack([tr, unp.nominal_values(Nr_cor), unp.std_devs(Nr_cor)]),
    fmt='%-20s',
    delimiter=' & ',
    header='t, Nr, Nrerr',
)

np.savetxt(
    'build/tab/Rhodium_diff.txt',
    np.column_stack([tr, unp.nominal_values(Nr_diff), unp.std_devs(Nr_diff)]),
    fmt='%-20s',
    delimiter=' & ',
    header='t, Nr_diff, Nr_differr',
)