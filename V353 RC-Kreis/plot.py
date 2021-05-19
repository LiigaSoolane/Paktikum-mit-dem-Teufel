import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from uncertainties import ufloat
import math
from tabulate import tabulate
from scipy.stats import sem

Sa = np.array([0.5 ,1   ,1.5 ,2   ,2.5 ,3   ,3.5 ,4 ,4.5 ,5   ,5.5 ,6   ,6.15,6.25,6.25]) #Spannung a) in V
t1 = np.array([0.1 ,0.19,0.23,0.3 ,0.42,0.6 ,0.8 ,1 ,1.3 ,1.6 ,2.09,3.2 ,4   ,6   ,8 ]) #Zeit in ms
t = t1/1000 #Zeit in s
F = np.array([20 ,120,220,320,420,620,820,1002 ,5020 ,10020,20020,30020,40020]) #Frequenz in Hz
Sb =np.array([3.25 ,2.3  ,1.7  ,1.325,1.025,0.73 ,0.59 ,0.48 ,0.105,0.052  ,0.026  ,0.01775,0.0135 ]) #Spannung b) in V
a = np.array([800 ,550 ,450 ,375 ,280 ,220 ,190 ,45  ,24  ,11.5 ,7.5 ,6   ]) #in mü s
b = np.array([8000,4200,3000,2300,1570,1160,960 ,195 ,96  ,48   ,31  ,24  ]) #b in mü s


#a)
y = np.log(1-(Sa/20))
def A(t,RC, U0):
#    return ((-1)/a)*t + b
    return U0*(1-np.exp(-1*(t/RC)))


par, cov = curve_fit(A, t, Sa, p0 = (0.004,20))
uncertainties = np.sqrt(np.diag(cov))

U0 = par[1]

for name, value, uncertainty in zip('ab', par, uncertainties): 
    print(f'{name} = {value:8.8f} ± {uncertainty:.8f}')

x = np.linspace(0,0.008)

plt.errorbar(t, Sa, fmt='r.', label=r'Messwerte')
plt.plot(x, A(x, par[0], par[1]), "r")
plt.xlabel(r'Zeit [s]')
plt.ylabel(r'U [V]')
plt.xscale("log")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('build/a).pdf')
plt.close('all')

#b)
def B(F, RC):
    #return np.exp((-a)*F+b)
    return 1 / (np.sqrt(1+(((2*np.pi*F)**2)*(RC**2))))

par2, cov2 = curve_fit(B, F, Sb/U0, p0 = (0.004))
uncertainties = np.sqrt(np.diag(cov2))

for name, value, uncertainty in zip('a', par2, uncertainties): 
    print(f'{name} = {value:8.8f} ± {uncertainty:.8f}')

x2 = np.linspace(20,40000)

plt.errorbar(F, Sb/U0, fmt='r.', label=r'Messwerte')
plt.plot(x2, B(x2, par2[0]), "r")
plt.xlabel(r'Frequenz [Hz]')
plt.ylabel(r'$U/U_0$')
plt.legend(loc='best')
plt.xscale("log")
plt.tight_layout()
plt.savefig('build/b).pdf')
plt.close('all')

#c)
Ph = (a/b)*2*np.pi
Fc = np.array([120,220,320,420,620,820,1002 ,5020 ,10020,20020,30020,40020]) #F ohne 20Hz

def C(F,RC):
    return np.arctan(-2*np.pi*F*RC)

par3, cov3 = curve_fit(C, Fc, Ph, p0 = (0.004))
uncertainties = np.sqrt(np.diag(cov3))

for name, value, uncertainty in zip('a', par3, uncertainties): 
    print(f'{name} = {value:8.8f} ± {uncertainty:.8f}')

plt.errorbar(Fc, Ph, fmt='r.', label=r'Messwerte')
plt.plot(x2, C(x2, par3[0]), "r")
plt.xlabel(r'Frequenz [Hz]')
plt.ylabel(r'$\Phi$')
plt.legend(loc='best')
plt.xscale("log")
plt.tight_layout()
plt.savefig('build/c).pdf')
plt.close('all')

#d)

def AU(F, Ph, RC):
    -np.sin(Ph)/(2*np.pi*F*RC)

RC = -par3[0]

AK = 1/ (np.sqrt(1+((2*np.pi*Fc)**2)*(RC**2)))

plt.figure()
plt.polar(Ph, AK, '.', label = 'Messdaten')

x = np.linspace(0, 100000, 10000)
phi = np.arcsin(((x*RC)/(np.sqrt(1+x**2*(RC)**2)))) 

A = 1/(np.sqrt(1+x**2*(-RC)**2))
plt.polar(phi, A, label = 'Berechnete Amplitude')
plt.legend()
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/d).pdf')

