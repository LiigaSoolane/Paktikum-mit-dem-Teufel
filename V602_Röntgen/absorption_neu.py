import numpy as np
import matplotlib.pyplot as plt
from uncertainties import ufloat
import uncertainties.unumpy as unp
from scipy.signal import find_peaks
from scipy.signal import peak_widths
from scipy.optimize import curve_fit

######################## constants
eV = 6.242e18
si = 1.602118e-19
h_si = 6.62607015e-34
h_ev = h_si * eV
c = 299792458
d_lif = 201.4e-12
a = 7.297e-3
Zs = np.array([30, 31, 35, 37, 38, 40])
R = 3.289e15
R_inf = 13.6 #eV

######################## read values
theta_zn, N_zn = np.genfromtxt('Zink.dat', unpack=True)
theta_ga, N_ga = np.genfromtxt('Gallium.dat', unpack=True)
theta_br, N_br = np.genfromtxt('Brom.dat', unpack=True)
theta_rb, N_rb = np.genfromtxt('Rubidium.dat', unpack=True)
theta_st, N_st = np.genfromtxt('Strontium.dat', unpack=True)
theta_zk, N_zk = np.genfromtxt('Zirkonium.dat', unpack=True)

######################## circulo de repeticion
names = ['Zink', 'Gallium', 'Brom', 'Rubidium', 'Strontium', 'Zirkonium']
thetas = [theta_zn, theta_ga, theta_br, theta_rb, theta_st, theta_zk]
Ns = [N_zn, N_ga, N_br, N_rb, N_st, N_zk]
theos = unp.uarray([9668.55, 10377.76, 13483.86, 15207.74, 16115.26, 18008.15], [15, 16, 19, 22, 23, 26])

E_K_array = np.empty(len(Zs))
E_K_errs = np.empty(len(Zs))
i = 0


for name, theta, N, Z, theo in zip(names, thetas, Ns, Zs, theos):
    print(f'-> {name}')

    # make plot
    plt.figure()
    plt.plot(theta, N, 'b.', label = f'{name}')
    plt.xlabel('Theta [°]')
    plt.ylabel('Anzahl Impulse')
    plt.legend()
    plt.savefig(f'build/plot_{name.lower()}.pdf')
    plt.clf()
    
    # Fehlerrechnung -> assume poisson distribution
    N_err = np.sqrt(N)
    N = unp.uarray(N, N_err)

    # Kante finden
    I_K_max = np.amax(unp.nominal_values(N))
    I_K_min = np.amin(unp.nominal_values(N))

    ind_max = np.where(unp.nominal_values(N) == I_K_max)
    ind_min = np.where(unp.nominal_values(N) == I_K_min)

    I_K_max = N[ind_max]
    I_K_min = N[ind_min]

    I_K = I_K_min + (I_K_max - I_K_min)/2

    I_K = I_K[0]
    print(f'Intensitätsmax. = {I_K}')

    test = np.abs(N - I_K)
    mask_test = test == np.amin(test)
    theta_K = theta[mask_test]

    # Sonderbehandlung für Gallium
    if Z == 31:
        theta_K = 17.35

    theta_K = ufloat(theta_K, 0.5)
    print(f'Winkel = {theta_K}') 

    # energias
    E_K = h_ev * c/(2 * d_lif * unp.sin(theta_K * np.pi / 180))
    E_K_array[i] = unp.nominal_values(E_K)
    E_K_errs[i] = unp.std_devs(E_K)

    E_dev = 100 * (theo - E_K)/theo
    print(f'Absorptionsenergie = {E_K}, Abweichung = {E_dev}')
    
    # Abschirmkonstanten
    sigma_K = Z - unp.sqrt(E_K / R_inf - a**2 * Z**4 / 4)
    sigma_theo = Z - unp.sqrt(theo / R_inf - a**2 * Z**4 / 4)
    sigma_dev = 100 * (sigma_K - sigma_theo)/sigma_theo
    print(f'Abschirmkonstante = {sigma_K}, Theoriewert = {sigma_theo}, Abweichung = {sigma_dev}')

    i += 1


######################## Rydberg

def y(x, m, n):
    return m*x+n


E_K_array_si = E_K_array * si

E = unp.sqrt(E_K_array)

# fit
#params, cov = np.polyfit(Zs, unp.nominal_values(E), deg = 1, full = False, cov = True)
params, cov = curve_fit(y, unp.nominal_values(E), Zs)
errs = np.sqrt(np.diag(cov))

a = ufloat(params[0], errs[0])
b = ufloat(params[1], errs[1])

print(f'Koeefizienten: ({a}, {b})')

# Rydberg durch Moseley
E_R_m = 1/(a**2)
R_m = E_R_m * si / h_si

#E_R_m *= eV

E_R_dev = 100 * (R_inf - E_R_m)/R_inf
R_dev = 100 * (R - R_m)/R

print(f'Rydbergfrequenz = {R_m}, Abweichung = {R_dev}')
print(f'Rydbergenergie = {E_R_m}, Abweichung = {E_R_dev}')

# make plot
x = np.linspace(8000, 17000, 1000)
x_ = np.sqrt(x)

plt.figure()
plt.plot(np.sqrt(E_K_array), Zs, color = 'salmon', marker = '.', linestyle = '', label = 'Moseley´sches Gesetz (diskret)')
plt.plot(x_, params[0] * x_ + params[1], color = 'steelblue', linestyle = '-', label = 'Lineare Regression')
plt.xlabel(r'$\sqrt{\text{ABsorptionsenergie} \mathbin{/} \text{eV}}$')
plt.ylabel(r'Z')
plt.legend()
plt.savefig('moseley.pdf')