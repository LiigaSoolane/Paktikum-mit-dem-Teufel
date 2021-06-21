import numpy as np
import matplotlib.pyplot as plt
from uncertainties import ufloat
import uncertainties.unumpy as unp
from scipy.signal import find_peaks
from scipy.signal import peak_widths

######################## constants
eV = 6.242e18
h_si = 6.62607015e-34
h_ev = h_si * eV
c = 299792458
d_lif = 201.4e-12
a = 7.297e-3
Z = 29
R = ufloat(1.097 * 10**7, 1.9 * 10**(-12)) 
R_inf = 13.6 #eV
E_abs = ufloat(8987.9, 15)
E_alpha_theo = 8048.11 
E_beta_theo = 8906.9

######################## read values
theta, N = np.genfromtxt('Emissionsspektrum.dat', unpack=True)

######################## Aufgabe 1: Darstellung etc. der K_alpha und K_beta Linien
k_alpha = ufloat(22.5, 0.1)
k_beta = ufloat(20.2, 0.1)

# calcular las energias
lam_alpha = 2 * d_lif * unp.sin(k_alpha * np.pi /180)
lam_beta = 2 * d_lif * unp.sin(k_beta * np.pi /180)

E_alpha = h_si * c / lam_alpha
E_beta = h_si * c / lam_beta

# to ev
E_alpha_ev = E_alpha * eV
E_beta_ev = E_beta * eV

print(f'E_a = {E_alpha_ev}, E_b = {E_beta_ev}')

# full width at half maximum 
peaks, _ = find_peaks(N, height = 1000)
widths = peak_widths(N, peaks, rel_height=0.5)

ind_alpha = np.array([int(widths[2][1]), int(widths[3][1])])
ind_beta = np.array([int(widths[2][0]), int(widths[3][0])])

fwhm_alpha = theta[ind_alpha[0]] - theta[ind_alpha[1]]
fwhm_beta = theta[ind_beta[0]] - theta[ind_beta[1]]

print(f'fwhm_a = {fwhm_alpha}, fwhm_b = {fwhm_beta}')

# Auflösungsvermögen
A_min = h_ev * c/(2 * d_lif * unp.sin(theta[ind_alpha[0]] * np.pi/180)) * eV
B_min = h_ev * c/(2 * d_lif * unp.sin(theta[ind_beta[0]] * np.pi/180)) * eV
A_max = h_ev * c/(2 * d_lif * unp.sin(theta[ind_alpha[1]] * np.pi/180)) * eV
B_max = h_ev * c/(2 * d_lif * unp.sin(theta[ind_beta[1]] * np.pi/180)) * eV

delta_E_alpha = A_max - A_min
delta_E_beta = B_max - B_min

print(f'Differenz_a = {delta_E_alpha}, Differenz_b = {delta_E_beta}')

A_alpha = E_alpha_ev/delta_E_alpha
A_beta = E_beta_ev/delta_E_beta

print(f'Auflöse_a = {A_alpha}, Auflöse_b = {A_beta}')

# Abschirmkonstanten

sigma_1 = - unp.sqrt(E_abs/R_inf) + Z
sigma_2 = - unp.sqrt( (E_alpha_ev + R_inf * (Z - sigma_1))/R_inf ) * 2 + Z
sigma_3 = - unp.sqrt( (E_beta_ev + R_inf * (Z - sigma_1))/ R_inf) * 3 + Z

print(f's_1 = {sigma_1:.3f}, s_2 = {sigma_2:.3f}, s_3 = {sigma_3:.3f}')

# Theoriewerte Abschirmkonst

sigma_2_theo = - unp.sqrt( (E_alpha_theo + R_inf * (Z - sigma_1))/R_inf ) * 2 + Z
sigma_3_theo = - unp.sqrt( (E_beta_theo + R_inf * (Z - sigma_1))/ R_inf) * 3 + Z

abw_s2 = 100 * (sigma_2 - sigma_2_theo)/sigma_2_theo
abw_s3 = 100 * (sigma_3 - sigma_3_theo)/sigma_3_theo

print(f'Theoriewert_alpha = {sigma_2_theo}, Abweichung = {abw_s2}')
print(f'Theoriewert_beta = {sigma_3_theo}, Abweichung = {abw_s3}')

######################## Plots!
kbeta = unp.nominal_values(k_beta)
kalpha = unp.nominal_values(k_alpha)

plt.figure()
plt.plot(theta, N, color = 'steelblue', marker = '.', label = 'Messwerte')

plt.hlines(widths[1], [theta[ind_beta[0]], theta[ind_alpha[0]]], [theta[ind_beta[1]], theta[ind_alpha[1]]], color = 'salmon', label = 'Full Width at Half Maximum')

plt.plot([kbeta, kbeta], [0, 1599.0], color='salmon', linestyle='--')
plt.scatter([kbeta], [1599.0], s=20, marker='o', color='salmon')
plt.annotate(r'$K_{\beta}$',
            xy = (kbeta, 1599.0), xycoords='data', xytext=(-50, -25),
            textcoords='offset points', fontsize=12,
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3, rad=.2"))
plt.plot([kalpha, kalpha], [0, 5050.0], color='salmon', linestyle='--')
plt.scatter([kalpha], [5050.0], s=20, marker='o', color='salmon')
plt.annotate(r'$K_{\alpha}$',
            xy = (kalpha, 5050.0), xycoords='data', xytext=(+10, -2),
            textcoords='offset points', fontsize=12,
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3, rad=.2"))

plt.plot([11.1, 11.1], [0, 420.0], color='steelblue', linestyle='--')
plt.scatter([11.1], [420.0], s=20, marker='o', color='steelblue')
plt.annotate(r'Bremsberg', 
            xy = (11.1, 420.0), xycoords='data', xytext=(-10, 20),
            textcoords='offset points', fontsize=12, 
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3, rad=.2"))

plt.xlabel(r'Einfallswinkel $\mathbin{/}°$')
plt.ylabel(r'Impulse pro $\si{\s}$')
plt.tight_layout()
plt.legend()
plt.savefig('emission_neu.pdf')
plt.clf()