import matplotlib.pyplot as plt
import numpy as np
from uncertainties import unumpy as unp
from scipy.optimize import curve_fit

# ------ Emissionsspektrum -------

# read values for Emission spectrum
theta, N = np.genfromtxt('EmissionCu.txt', unpack=True)

# set k_alpha and k_beta line variables
kbeta = 20.2
kalpha = 22.5

#plot Emission spectrum
plt.plot(theta, N, 'r.', label='Messwerte')
plt.plot(theta, N, 'b-', linewidth=0.5)

# mark Kalpha and Kbeta lines
plt.plot([kbeta, kbeta], [0, 1599.0], color='red', linestyle='--')
plt.scatter([kbeta], [1599.0], s=20, marker='o', color='red')
plt.annotate(r'$K_{\beta}$',
            xy = (kbeta, 1599.0), xycoords='data', xytext=(-50, -25),
            textcoords='offset points', fontsize=12,
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3, rad=.2"))
plt.plot([kalpha, kalpha], [0, 5050.0], color='red', linestyle='--')
plt.scatter([kalpha], [5050.0], s=20, marker='o', color='red')
plt.annotate(r'$K_{\alpha}$',
            xy = (kalpha, 5050.0), xycoords='data', xytext=(+5, +20),
            textcoords='offset points', fontsize=12,
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3, rad=.2"))

# mark Bremsberg
plt.plot([11.1, 11.1], [0, 420.0], color='red', linestyle='--')
plt.scatter([11.1], [420.0], s=20, marker='o', color='red')
plt.annotate(r'Bremsberg', 
            xy = (11.1, 420.0), xycoords='data', xytext=(-10, 20),
            textcoords='offset points', fontsize=12, 
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3, rad=.2"))

plt.xlabel('Einfallswinkel in °')
plt.ylabel('Anzahl Impulse pro s')
plt.tight_layout()
plt.legend()
plt.savefig('emission.pdf')
plt.clf()

# ------ Transmission als Funktion der Wellenlänge  -------

# read values for with and without Aluminum absorber
a_0, N_0_ = np.genfromtxt('ComptonOhne.txt', unpack=True)
a_Al, N_Al_ = np.genfromtxt('ComptonAl.txt', unpack=True)
Tau = 90
Tau *= 10**(-6)
d_LiF = 201.4 
d_LiF *= 10**(-12)

N_0_err = np.sqrt(N_0_)
N_Al_err = np.sqrt(N_Al_)

N_0 = unp.uarray(N_0_, N_0_err)
N_Al = unp.uarray(N_Al_, N_Al_err)

# Totzeitkorrektur
I_o = N_0 / (1 - Tau * N_0)
I_Al = N_Al / (1 - Tau * N_Al)

# Transmission bestimmen
T = I_Al / I_o

# calcular taman~o de las ondas

lam = 2 * d_LiF * np.sin(a_0 * np.pi / 180)

# crear lineare Regression
params, cov = np.polyfit(lam, unp.nominal_values(T), deg=1, cov=True)
errs = np.sqrt(np.diag(cov))
for name, value, error in zip('ab', params, errs):
    print(f'{name} = {value:.3f} +- {error:.3f}')

# grafico
l_start = 2* d_LiF * np.sin(7 * np.pi /180)
l_fin = 2 * d_LiF *np.sin(10* np.pi/180)
lam_plot = np.linspace(l_start, l_fin)

plt.clf()
plt.plot(lam, unp.nominal_values(T), 'r.', label='Messwerte')
plt.plot(lam_plot, params[0]*lam_plot + params[1], '-', label='Lineare Regression')
plt.errorbar(lam, unp.nominal_values(T), yerr=unp.std_devs(T), fmt='r_')
plt.legend()
plt.xlabel(r'Wellenlänge $\lambda$')
plt.ylabel(r'Transmission')
plt.savefig('transmission.pdf')

# ----- Compton-Wellenlänge bestimmen ------

# Transmissionen
I_0 = 2731.0
I_1 = 1180.0
I_2 = 1024.0

T_1 = I_1/I_0
T_2 = I_2/I_0

print(f'T_1 = {T_1:.3f}, T_2 = {T_2:.3f}')

# calcular tamanos de ondas correspendientes
lam_1 = (T_1 - params[1])/params[0]
lam_2 = (T_2 - params[1])/params[0]

lam_c = lam_2 - lam_1

print(f'Lambda 1 = {lam_1}')
print(f'Lambda 2 = {lam_2}')
print(f'Compton = {lam_c}')