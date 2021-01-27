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
print(Nvtrue)
# --> linregress
params1, cov1 = np.polyfit(tv, np.log(unp.nominal_values(Nvtrue)), deg = 1, cov = True)
errors1 = np.sqrt(np.diag(cov1))
params1[0] *= (-1)
print(params1, errors1)

# -> Plot
x1 = np.linspace(0, 1230)

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

