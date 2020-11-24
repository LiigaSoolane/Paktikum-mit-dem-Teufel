import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from uncertainties import ufloat

p,t=np.genfromtxt("messwerte.txt", unpack=True)
t=t+273.15
p=p*10**(2)

x_plot=1/t
y_plot=np.log(p)
xfit = np.linspace(x_plot[0],x_plot[10],80)
params, cov = np.polyfit(x_plot, y_plot, deg=1, cov=True)
print('a,b = ',params)
errors2 = np.sqrt(np.diag(cov))
print('aerr,berr=',errors2)
L = ufloat(-params[0]*8.31446261815324,errors2[0]*8.31446261815324)
error=errors2[0]*8.31446261815324
print('Lerr=',error)
print('L = ', L)

plt.figure(0)
plt.plot(x_plot, y_plot, 'r*', label='Messdaten des Reservoir 1')
plt.plot(xfit, params[0]*xfit+params[1], 'g-', label='fit')
plt.ylabel(r'ln($\dfrac{p}{p_0}$)')
plt.xlabel(r'1/Temperatur [1/K]')
plt.legend(loc='best')
plt.savefig('plot1.pdf')



