import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from uncertainties import ufloat

p,t=np.genfromtxt("messwerte.txt", unpack=True)
t=t+273.15
p=p*10**(2)
r = 8.31446261815324

np.savetxt("aufgabe4.txt", np.column_stack([p, t]), fmt = "%10.4f", delimiter = " & ", newline = " \\\ ", header = " Tka Tsa ")

x_plot=1/t
y_plot=np.log(p)
xfit = np.linspace(x_plot[0],x_plot[49],80)
params, cov = np.polyfit(x_plot, y_plot, deg=1, cov=True)
print('a,b = ',params)
errors2 = np.sqrt(np.diag(cov))
print('aerr,berr=',errors2)
L = ufloat(-params[0]*8.31446261815324,errors2[0]*8.31446261815324)
print('L = ', L)

plt.figure(0)
plt.plot(x_plot, y_plot, 'r*', label='Messdaten')
plt.plot(xfit, params[0]*xfit+params[1], 'g-', label='fit')
plt.ylabel(r'ln($\dfrac{p}{p_0}$)')
plt.xlabel(r'1/Temperatur [1/K]')
plt.legend(loc='best')
plt.savefig('build/plot.pdf')


p,t=np.genfromtxt("messwerte.txt", unpack=True)
t=t+273.15
p=p*10**(5)
S=0.9
params, cov = np.polyfit(t, p, deg=3, cov=True)
print('a,b = ',params)
errors2 = np.sqrt(np.diag(cov))
print('aerr,berr=',errors2)
a = ufloat(params[0], errors2[0])
b = ufloat(params[1], errors2[1])
c = ufloat(params[2], errors2[2])
d = ufloat(params[3], errors2[3])
pt= params[0]*t**3+params[1]*t**2+params[2]*t+params[3]
upt = a*t**3+b*t**2+c*t+d
dpt= 3*params[0]*t**2+2*params[1]*t+params[2]
udpt= 3*a*t**2+2*b*t+c

plt.figure(1)
plt.plot(t, p, 'r*', label='Messdaten')
plt.plot(t, pt, 'g-', label='fit')
plt.ylabel(r'$p$ [Pa]')
plt.xlabel(r'Temperatur [K]')
plt.legend(loc='best')
plt.savefig('build/plot1.pdf')

#print(lt[2])
#lt = dpt/pt*(r*t/2+(-(r*t/2)**2+S*pt)**(1/2))
#plt.figure(2)
#plt.plot(t, lt, 'r*', label='L')
#plt.ylabel(r'L(T)')
#plt.xlabel(r'Temperatur [K]')
#plt.legend(loc='best')
#plt.savefig('build/plot2.pdf')
#
#lt = dpt/pt*(r*t/2-((r*t/2)**2-S*pt)**(1/2))
#plt.figure(3)
#plt.plot(t, lt, 'r*', label='L')
#plt.ylabel(r'L(T)')
#plt.xlabel(r'Temperatur [K]')
#plt.legend(loc='best')
#plt.savefig('build/plot3.pdf')
#