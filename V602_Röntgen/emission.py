import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
from uncertainties import ufloat
from uncertainties import unumpy as unp
from scipy.optimize import curve_fit
from scipy.signal import peak_widths
from scipy.signal import find_peaks

# valores de la naturaleza
h = 6.62607015 * 10**(-34)
c = 299792458
d_LiF = 201.4 
d_LiF *= 10**(-12)
R = ufloat(1.097 * 10**7, 1.9 * 10**(-12)) 
Z = 29

# ------ Emissionsspektrum -------

# read values for Emission spectrum
theta, N = np.genfromtxt('Emissionsspektrum.dat', unpack=True)

# set k_alpha and k_beta line variables
kbeta = 20.2 
kalpha = 22.5

kalpha_ = unp.uarray(kalpha, 0.1)
kbeta_ = unp.uarray(kbeta, 0.1)

# K_alpha y K_beta energias
lam_alpha = 2 * d_LiF * unp.sin(kalpha_ * np.pi / 180)
lam_beta = 2 * d_LiF * unp.sin(kbeta_ * np.pi / 180)

E_alpha = h * c / lam_alpha
E_beta = h * c / lam_beta

E_alpha *= 6.242 * 10**(18)
E_beta *= 6.242 * 10**(18)

print(f'E_a = {E_alpha}')
print(f'E_b = {E_beta}')

# full width at half maximum
peaks = unp.uarray([unp.nominal_values(kalpha_), unp.nominal_values(kbeta_)], [unp.std_devs(kalpha_), unp.std_devs(kbeta_)])

peaks, _ = find_peaks(N, height = 1000)
widths = peak_widths(N, peaks, rel_height=0.5)
print(widths[0], widths[1], widths[2])

ind_one = [int(widths[2][0]), int(widths[2][1])]
ind_two = [int(widths[3][0]), int(widths[3][1])]

print(f'Full Width Half Maximum, k_B: {widths[0][0]}')
print(f'Full Width Half Maximum, k_a: {widths[0][1]}')

# Auflösungsvermögen
E_alpha_min = h * c/(2 * d_LiF * unp.sin(theta[ind_one[0]] * np.pi / 180))   * 6.242 * 10**(18)
E_alpha_max = h * c /(2 * d_LiF * unp.sin(theta[ind_two[0]] * np.pi / 180)) * 6.242 * 10**(18)
E_beta_min = h * c /(2 * d_LiF * unp.sin(theta[ind_one[1]] * np.pi / 180)) * 6.242 * 10**(18)
E_beta_max = h * c /(2 * d_LiF * unp.sin(theta[ind_two[1]] * np.pi / 180)) * 6.242 * 10**(18)

delta_E_alpha = E_alpha_min - E_alpha_max
delta_E_beta = E_beta_min - E_beta_max
print(delta_E_beta, delta_E_alpha)
A_alpha = E_alpha/delta_E_alpha
A_beta = E_beta/delta_E_beta

print(f'A_a = {A_alpha:.3f}, A_b = {A_beta:.3f}')

# Abschirmkonstanten
l = 3
m = 2
n = 1

R = 13.6 * 6.242 * 10**(18)
E_abs = ufloat(8987.9, 15)
E_abs *= 6.242 * 10**(18)

sigma_1 = - unp.sqrt(E_abs/R) + Z 
sigma_2 = - unp.sqrt((- E_alpha + R * (Z - sigma_1) / n**2)/R ) * m + Z
sigma_3 = - unp.sqrt((- E_beta + R * (Z - sigma_1) / n**2)/R ) * l + Z

print(f's_1 = {sigma_1:.3f}, s_2 = {sigma_2:.3f}, s_3 = {sigma_3:.3f}')

E_alpha_t = 8048.11
E_beta_t = 8906.9

sigma_1_t = - unp.sqrt(E_abs/R) + Z 
sigma_2_t = - unp.sqrt((- E_alpha_t + R * (Z - sigma_1_t) / n**2)/R ) * m + Z
sigma_3_t = - unp.sqrt((- E_beta_t + R * (Z - sigma_1_t) / n**2)/R ) * l + Z

print(f'Theoriewerte s_1 = {sigma_1_t:.3f}, s_2 = {sigma_2_t:.3f}, s_3 = {sigma_3_t:.3f}')

#### Montones de Plots. En realidad es solo uno pero tiene miles de partes
# plot Emission spectrum
plt.plot(theta, N, 'r.', label='Messwerte')
plt.plot(theta, N, 'b-', linewidth=0.5)

# plot full width at half maximum
plt.hlines(widths[1], [theta[ind_one[0]], theta[ind_one[1]]], [theta[ind_two[0]], theta[ind_two[1]]], label='Full Width at Half Maximum')

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
            xy = (kalpha, 5050.0), xycoords='data', xytext=(+10, -2),
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

# diferencia procentual a los valores de literatura
E_t_a = unp.uarray(8048.11, 45)
E_t_b = unp.uarray(8906.9, 12)

abw_alpha = 100 * (E_alpha - E_t_a)/E_t_a
abw_beta = 100 * (E_beta - E_t_b)/E_t_b

print(f'% abw. alpha: {abw_alpha}')
print(f'% abw. beta: {abw_beta}')