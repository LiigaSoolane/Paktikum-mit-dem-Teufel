import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import scipy.constants as const 
import uncertainties.unumpy as unp
from uncertainties import ufloat

###################### solo a un lado

x_rund, D_rund = np.genfromtxt('einseitig_rund_data.txt', unpack = True)
x_eckig, D_eckig = np.genfromtxt('einseitig_eckig_data.txt', unpack = True)

x_rechts, D_rechts = np.genfromtxt('beidseitig_rechts_data.txt', unpack=True)
x_links, D_links = np.genfromtxt('beidseitig_links_data.txt', unpack=True)

x_rechts *= 1e-2
x_links *= 1e-2

D_rechts *= 1e-3
D_links *= 1e-3

x_rund *= 1e-2
x_eckig *= 1e-2

D_rund *= 1e-3
D_eckig *=1e-3

# tu sei pazzo, mica van gogh
L_rund = 58.05e-2
L_eckig = 60e-2

G_rund = 356e-3
G_eckig = 535.7e-3

Q_rund = np.pi * (1e-2 / 2)**2
Q_eckig = 1.025e-2 * 1.025e-2

rho_rund = G_rund/(Q_rund * L_rund)
rho_eckig = G_eckig/(L_eckig * Q_eckig)

rho_rund *= 1e3 / (1e2)**(3)
rho_eckig *= 1e3 / (1e2)**(3)

I_rund = np.pi * (1e-2)**4 / 64
I_eckig = (1.025e-2)**4 /12

x_pos_eckig = 54e-2
x_pos_rund = 52e-2

M_rund = 528e-3
M_eckig = 650e-3
M_beid = 751e-3

F_eckig = M_eckig * const.g
F_rund = M_rund * const.g
F_beid = M_beid * const.g

x = np.linspace(0, 0.1, 1000)
y = np.linspace(0, 0.2, 1000)

# tu no sei piu sano, tu sei pazzo
print(f'rho, rund = {rho_rund}, eckig = {rho_eckig}')

######################## make fit

def hans_günther(L, x):
    return L * x**2 - x**3 / 3

term_rund = hans_günther(x_pos_rund, x_rund)
term_eckig = hans_günther(x_pos_eckig, x_eckig)

######################## now really make the fit y date prisa porque no tienes tiempo para nada
params_rund, cov_rund = np.polyfit(term_rund, D_rund, deg = 1, cov = True)
errs_rund = np.sqrt(np.diag(cov_rund))

a_rund = ufloat(params_rund[0], errs_rund[0])
b_rund = ufloat(params_rund[1], errs_rund[1])

params_eckig, cov_eckig = np.polyfit(term_eckig, D_eckig, deg = 1, cov = True)
errs_eckig = np.sqrt(np.diag(cov_eckig))

a_eckig = ufloat(params_eckig[0], errs_eckig[0])
b_eckig = ufloat(params_eckig[1], errs_eckig[1])

print(f'für rund: a = {a_rund}, b = {b_rund}')
print(f'für eckig: a = {a_eckig}, b = {b_eckig}')

plt.figure()
plt.plot(term_rund, D_rund, color='salmon', linestyle='', marker='.', label='Messwerte runder Stab (Messing)')
plt.plot(term_eckig, D_eckig, color='steelblue', linestyle='', marker='.', label='Messwerte eckiger Stab (Kupfer)')

plt.plot(x, params_eckig[0]*x+params_eckig[1], color='steelblue', linestyle='-', label='Lineare Ausgleichsgerade eckiger Stab')
plt.plot(x, params_rund[0]*x+params_rund[1], color='salmon', linestyle='-', label='Lineare Ausgleichsgerade runder Stab')

plt.legend()
plt.xlabel(r'$L x^2 - \frac{x^3}{3}')
plt.ylabel(r'$D(x)$')
plt.tight_layout()

plt.savefig('plots_einseitig.pdf')
plt.clf()


######################## e tu, dove vai a ballare?
E_rund = F_rund / (2 * a_rund * I_rund)
E_eckig = F_eckig / (2 * a_eckig * I_eckig)

print(f'Elast rund = {E_rund:.3f}')
print(f'Elast eckig = {E_eckig:.3f}')

######################## vieni a ballare en Puglia, Puglia, Puglia
######################## tremulo come una foglia, foglia, foglia

def karl_heinz(L, x):
    return 3 * L**2 * x - 4 * x**3

def ernst_wilhelm(L, x):
    return 4 * x**3 - 12 * L * x**2 + 9 * L**2 * x - L**3

term_links = ernst_wilhelm(0.56, x_links)
term_rechts = karl_heinz(0.56, x_rechts)

##################### dicono che mi spediranno en cielo come una soyuz
params_rechts, cov_rechts = np.polyfit(term_rechts, D_rechts, deg = 1, cov = True)
errs_rechts = np.sqrt(np.diag(cov_rechts))

a_rechts = ufloat(params_rechts[0], errs_rechts[0])
b_rechts = ufloat(params_rechts[1], errs_rechts[1])

params_links, cov_links = np.polyfit(term_links, D_links, deg = 1, cov = True)
errs_links = np.sqrt(np.diag(cov_links))

a_links = ufloat(params_links[0], errs_links[0])
b_links = ufloat(params_links[1], errs_links[1])

print(f'für rechts: a = {a_rechts}, b = {b_rechts}')
print(f'für links: a = {a_links}, b = {b_links}')

##################### la fine di gaia non arrivera
plt.figure()
plt.plot(term_rechts, D_rechts, linestyle='', marker='.', color='salmon', label='Messwerte')
plt.plot(y, params_rechts[0]*y + params_rechts[1], linestyle='-', color='steelblue', label='Lineare Ausgleichsgerade')

plt.legend()
plt.xlabel(r'$4 x^3 - 12 L x^2 - 9 L^2 x - L^3$')
plt.ylabel(r'D(x)')
plt.tight_layout()

plt.savefig('plot_rechts.pdf')
plt.clf()

plt.figure()
plt.plot(term_links, D_links, linestyle='', marker='.', color='salmon', label='Messwerte')
plt.plot(y, params_links[0]*y + params_links[1], linestyle='-', color='steelblue', label='Lineare Auslgleichsgerade')

plt.legend()
plt.xlabel(r'$ 3 L^2 x - 4 x^3 $')
plt.ylabel(r'D(x)')
plt.tight_layout()

plt.savefig('plot_links.pdf')
plt.clf()

##################### ok va bene, avrai ragione tu
E_rechts = F_beid / (48 * a_rechts * I_rund)
E_links = F_beid / (48 * a_links * I_rund)

print(f'Elast rechs = {E_rechts}')
print(f'Elast links = {E_links}')

# que falta? los I_eckig y los F, pienso
# y por supuesto todos los textos de latex
# no lo aguanto mas

# sono dalla parte del toro

E_theo_rund = 196e9

abw_rund = 100 * (E_theo_rund - E_rund)/E_theo_rund
abw_rechts = 100 * (E_theo_rund - E_rechts)/E_theo_rund
abw_links = 100 * (E_theo_rund - E_links)/E_theo_rund

E_eckig_theo = ufloat(100, 22)
E_eckig_theo *= 1e9

abw_eckig = 100 * (E_eckig_theo - E_eckig)/E_eckig_theo

print(f'abw einseitig: rund = {abw_rund}, eckig = {abw_eckig}')
print(f'abw beidseitig: rechts = {abw_rechts}, links = {abw_links}')