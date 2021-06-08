import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from scipy.optimize import curve_fit

# set const
h = 6.62607015 * 10**(-34)
c = 299792458
d_LiF = 201.4 
d_LiF *= 10**(-12)
alpha = 7.297 * 10**(-3)
h_eV = 4.1357 * 10**(-15)


# read vals
theta_zn, N_zn = np.genfromtxt('Zink.dat', unpack=True)
theta_ga, N_ga = np.genfromtxt('Gallium.dat', unpack=True)
theta_br, N_br = np.genfromtxt('Brom.dat', unpack=True)
theta_rb, N_rb = np.genfromtxt('Rubidium.dat', unpack=True)
theta_st, N_st = np.genfromtxt('Strontium.dat', unpack=True)
theta_zk, N_zk = np.genfromtxt('Zirkonium.dat', unpack=True)

names = ['Zink', 'Gallium', 'Brom', 'Rubidium', 'Strontium', 'Zirkonium']
thetas = [theta_zn, theta_ga, theta_br, theta_rb, theta_st, theta_zk]
Ns = [N_zn, N_ga, N_br, N_rb, N_st, N_zk]
theos = unp.uarray([9668.55, 10377.76, 13483.86, 15207.74, 16115.26, 18008.15], [15, 16, 19, 22, 23, 26])
Zs = [30, 31, 35, 37, 38, 40]

E_K_array = np.empty(len(Zs))
i = 0

# make loop for all other stuffs
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

    # Kante
    I_K_max = np.amax(N)
    I_K_min = np.amin(N)

    I_K = I_K_min + (I_K_max - I_K_min)/2

    print(f'Intensitätsmax. = {I_K}')
    test = I_K - theta
    mask_test = test == np.amin(test)
    theta_K = theta[mask_test]
    print(f'Winkel = {theta_K}')

    E_K = h * c/(2 * d_LiF * np.sin(theta_K * np.pi / 180)) * 6.242 * 10**(18)
    E_K_array[i] = E_K
    print(f'Absorptionsenergie = {E_K}')

    R_inf = 13.6
    # Abschirmkonstanten
    sigma_K = Z - np.sqrt(E_K / R_inf - alpha**2 * Z**4 / 4)
    sigma_theo = Z - unp.sqrt(theo / R_inf - alpha**2 * Z**4 / 4)
    sigma_dev = 100 * (sigma_K - sigma_theo)/sigma_theo
    print(f'Abschirmungskonstante = {sigma_K}, Theoriewert = {sigma_theo}, abweichung = {sigma_dev}')

    # Deviations
    dev = 100 * (E_K - theo)/theo
    print(f'%Abweichung von Energie = {dev}')

    i += 1


def y(x, m, n):
    return m * x + n

E_n = E_K_array * 1.60218 * 10**(-19)
E_jasd = np.sqrt(E_K_array)
E_K_sq = np.sqrt(E_n)
params, cov = curve_fit(y, E_jasd, Zs)
err = np.sqrt(np.diag(cov))
a = unp.uarray(params[0], err[0])
b = unp.uarray(params[1], err[1])

print(a, b)

# Rydberg bestimmt durch Moseley
#E_R_m = 6.242 * 10**(18)/(a**2) 
E_R_m = 1/a**2
R_m = E_R_m * 1.60218 * 10**(-19)/h

dev_ryd = 100 * (13.6 - E_R_m)/13.6 

print(f'Rydbergfrequenz = {R_m}, Rydbergenergie = {E_R_m}, deviation = {dev_ryd}')

x = np.linspace(8000, 17000, 1000)

#params, cov = curve_fit(y, E_K_array, Zs)
# Wurzel e - z Plot
plt.plot(np.sqrt(E_K_array), Zs, 'b.', label = 'Moseley´sches Gesetz (diskret)')
plt.plot(np.sqrt(x), params[0]*np.sqrt(x) + params[1], 'b-', label = 'Lineare Regression')
plt.xlabel('Absorptionsenergie [eV]')
plt.ylabel('Z')
plt.legend()
plt.savefig('moseley.pdf')