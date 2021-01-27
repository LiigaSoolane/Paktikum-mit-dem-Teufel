import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
from scipy.optimize import curve_fit

tv, Nv_ = np.genfromtxt('Vanadium.dat', unpack=True)
tr, Nr_ = np.genfromtxt('Rhodium.dat', unpack=True)

Nv_err = np.sqrt(Nv_)
Nv = unp.uarray(Nv_, Nv_err)

Nu_ = np.array([129.0, 143.0, 144.0, 136.0, 126.0, 158.0])
tu = 300

Nu_ *= 10**(-1)
Nu_err = np.sqrt(Nu_)
Nu = unp.uarray(Nu_, Nu_err)

Nu_mean = np.mean(Nu)
# Vanadium
# -> Curve fit
Nvtrue = Nv - Nu_mean
#print(Nvtrue)
# --> linregress
params1, cov1 = np.polyfit(tv, np.log(unp.nominal_values(Nvtrue)), deg = 1, cov = True)
errors1 = np.sqrt(np.diag(cov1))
params1[0] *= (-1)
print(params1, errors1)

# -> Plot
x1 = np.linspace(0, 1230)

plt.plot(x1, params1[1] -params1[0]*x1, 'b-', label="Lineare Regression")
plt.errorbar(tv, unp.nominal_values(unp.log(Nvtrue)), yerr = unp.std_devs(unp.log(Nvtrue)), fmt='g.', label="Messwerte mit Fehlern")
plt.yscale('log')
plt.xlabel('t [s]')
plt.ylabel('ln(N)')
plt.legend()
plt.tight_layout()
plt.savefig('plot.pdf')

# -> Halbwertszeit

Tv = np.log(2)/params1[0]
Tv_err = np.sqrt(np.log(2)**2 * errors1[0]**2)

print(Tv, Tv_err)

# -> Halbwertszeit genauer
params2, cov2 = np.polyfit(tv[0:15], np.log(unp.nominal_values(Nvtrue[0:15])), deg = 1, cov = True)
errors2 = np.sqrt(np.diag(cov2))
params2[0] *= (-1)
print(params2, errors2)

Tv = np.log(2)/params2[0]
Tv_err = np.sqrt(np.log(2)**2 * errors2[0]**2)
print(Tv, Tv_err)

# Rhodium

Num = np.mean(Nu_)

Nrtrue = unp.uarray(Nr_ - 1/2 * Num, np.sqrt(Nr_-1/2 * Num))

# linregress
params3, cov3 = np.polyfit(tr, np.log(unp.nominal_values(Nrtrue)), deg = 1, cov = True)
errors3 = np.sqrt(np.diag(cov3))
params3[0] *= (-1)
print(params3, errors3)

# setze willkürlich t* = 450s, also slot 26
tstern = 420

paramsl, covl = np.polyfit(tr[28:], np.log(unp.nominal_values(Nrtrue[28:])), deg= 1, cov = True)
errorsl = np.sqrt(np.diag(covl))
paramsl[0] *= -1
print (paramsl, errorsl)

Nkurz = Nrtrue -paramsl[1] * np.exp(- tr * paramsl[0])

#print(Nkurz) # muestra que desde el principio hasta numero 16 domina la corta mierda
tmax = 250

# y otra vez porque no ha sido bastante
paramsk, covk = np.polyfit(tr[:17], np.log(unp.nominal_values(Nkurz[:17])), deg=1, cov=True)
errorsk = np.sqrt(np.diag(covk))
paramsk[0] *= -1
print(paramsk, errorsk)

# declarar funcion aditiva
def add(x):
    return paramsk[1] - paramsk[0] * x + paramsl[1] - paramsl[0] * x

# -> Plot
x2 = np.linspace(0, 630)
xl = np.linspace(tstern, 630)
xk = np.linspace(0, tmax)

plt.clf()
plt.plot(x2, params3[1] -params3[0]*x2, 'b-', label="Lineare Regression")
plt.errorbar(tr, unp.nominal_values(unp.log(Nrtrue)), yerr = unp.std_devs(unp.log(Nrtrue)), fmt='b.', label="Messwerte mit Fehlern")
plt.plot(xl, paramsl[1] -paramsl[0]*xl, 'r--', label="Fit für das Eine")
plt.plot(xk, paramsk[1] -paramsk[0]*xk, 'y--', label="Fit für das Andere")
#plt.plot(x2, add(x2), 'g-', label="fuck this/me/you/off/everything/my life")
plt.yscale('log')
plt.xlabel('t [s]')
plt.ylabel('ln(N)')
plt.legend()
plt.tight_layout()
plt.savefig('plot2.pdf')

# -> Halbwertszeiten
Tl = np.log(2)/paramsl[0]
Tl_err = np.sqrt(np.log(2)**2 * errorsl[0]**2)

Tk = np.log(2)/paramsk[0]
Tk_err = np.sqrt(np.log(2)**2 * errorsk[0]**2)

print(Tl, Tl_err)
print(Tk, Tk_err)

