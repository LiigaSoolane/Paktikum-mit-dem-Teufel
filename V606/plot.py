import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from uncertainties import ufloat
import math
from tabulate import tabulate
from scipy.stats import sem

def Gauss(x,a, b):
    return np.exp(-b*(x+a)**2)

U, F = np.genfromtxt('data.txt', unpack=True)
table ={'Spannung [V]': U, 'Frequenz [Hz]': F}
print(tabulate (table, tablefmt="latex"))
U = U/1000
F = F*1000
Umax = max(U)
U = U /Umax

params, covariance_matrix = curve_fit(Gauss, F, U, p0=(-36000, 0.00001))
uncertainties = np.sqrt(np.diag(covariance_matrix))
for name, value, uncertainty in zip('abc', params, uncertainties): 
    print(f'{name} = {value:8.8f} ± {uncertainty:.8f}')

x_dings = np.linspace(F[0], F[len(F)-1], 10000)
y_dings = Gauss(x_dings, *params)
plt.plot(F, U, '*', label="Messwerte")
plt.plot(x_dings, y_dings,'c', label='Curvefit')
plt.xlabel(r'Frequenz [Hz]')
plt.ylabel('U/U_A')
plt.legend(loc='best')
plt.savefig('build/plot.pdf')

nu_p = 0
nu_m = 0
nu = 0
sigma = 1/np.sqrt(2)
ymax=max(y_dings)
for i in range(10000):
    if nu_m==0:
        if y_dings[i] > sigma:
            nu_p=i
            nu_m=1
    else :
        if y_dings[i]==ymax:
            nu = i
        if y_dings[i]>1/np.sqrt(2):
            nu_m=i
G=nu/(nu_m-nu_p)
print(f"nu_p = ", nu_p)
print(f"nu = ", nu)
print(f"nu_m = ", nu_m)
print(f"Güte = ", G)

#########################################################################################

F = 0.866
gdohmv,_, gdohmn, gdspn = np.genfromtxt('gd.txt', unpack=True)
dyohmv,_, dyohmn, dyspn = np.genfromtxt('gd.txt', unpack=True)
gdspv=3.35/1000
dyspv=3.35/1000
gdohmv=gdohmv*5/1000
gdohmn=gdohmn*5/1000
dyohmv=dyohmv*5/1000
dyohmn=dyohmn*5/1000
gdspn=gdspn/1000
dyspn=dyspn/1000

table ={'1': gdohmv, '3': gdohmn, '4':gdspn}
print(tabulate (table, tablefmt="latex"))
table ={'1': dyohmv, '3': dyohmn, '4':dyspn}
print(tabulate (table, tablefmt="latex"))

gdohmv=ufloat(np.mean(gdohmv), sem(gdohmv))
gdohmn=ufloat(np.mean(gdohmn), sem(gdohmn))
dyohmv=ufloat(np.mean(dyohmv), sem(dyohmv))
dyohmn=ufloat(np.mean(dyohmn), sem(dyohmn))
gdspn=ufloat(np.mean(gdspn), sem(gdspn))
dyspn=ufloat(np.mean(dyspn), sem(dyspn))

gdM=14.08
gdp=7.40
gdl=15.7
gdQ=gdM/(gdp*gdl)*10000

dyM=14.39
dyp=7.8
dyl=15.7
dyQ=dyM/(dyp*dyl)*10000

gdxU=4*F*gdspv/(gdQ*gdspn)
gdxR=4*F*(gdohmv-gdohmn)/(gdQ*gdohmv)
dyxU=4*F*dyspv/(dyQ*dyspn)
dyxR=4*F*(dyohmv-dyohmn)/(dyQ*dyohmv)

print(f"gdxU = ", gdxU)
print(f"gdxR = ", gdxR)
print(f"dyxU = ", dyxU)
print(f"dyxR = ", dyxR)

############################################################################

#Therorie

m0 = 1.257e-6
mb = (-1.602e-19 * 6.626e-34)/(4*np.pi * 9.109e-31) 
kb = 1.38e-23
T = 17 + 273.15 
A = 6.02214086e23

#Gd2O3
print('Gd2O3')
J = 3.5
S = 3.5
L = 0
M = 362.4982 
dichte = 7.4e6 
N = 2*(A * dichte)/(M)
gj = (3*J*(J+1)+ (S*(S+1)-L*(L+1)))/(2*J*(J+1))
print('g_j:', gj )
chi_t = (m0 * (mb)**2 * (gj)**2 * N * J * (J+1) )/(3 * kb * T)
print('chi_T:', chi_t)

deltau = (chi_t - gdxU) / chi_t
print('ABWEICHUNG WiDERSTAND', deltau)
deltar = (chi_t - gdxR) / chi_t
print('ABWEICHUNG SPANNUNG', deltar)

# Dy2O3
print('Dy2O3')
J = 7.5
S = 2.5
L = 5
M = 372.9982 
dichte = 7.8e6 
N = 2*(A * dichte)/(M) 
gj = (3*J*(J+1)+ (S*(S+1)-L*(L+1)))/(2*J*(J+1))
print('g_j:', gj )
chi_t = (m0 * (mb)**2 * (gj)**2 * N * J * (J+1) )/(3 * kb * T)
print('chi_T:', chi_t)

deltau = (chi_t - dyxU) / chi_t
print('ABWEICHUNG WiDERSTAND', deltau)
deltar = (chi_t - dyxR) / chi_t
print('ABWEICHUNG SPANNUNG', deltar)

