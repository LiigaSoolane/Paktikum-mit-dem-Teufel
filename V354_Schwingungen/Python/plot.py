import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from uncertainties import ufloat
from uncertainties import unumpy
import math

#########################################################################################

def expo(x, a, b):
    return a*np.exp(-2*math.pi*b*x)

#väärtuseid
L = ufloat(0.01278, 0.00009)
C = ufloat(2.066*10**(-9), 0.006*10**(-9))
time1, voltage1 = np.genfromtxt('Messung1.txt', unpack=True)
time2, voltage2 = np.genfromtxt('Messung2.txt', unpack=True)
voltagesisse, voltagevälja, freq= np.genfromtxt('Messung3.txt', unpack=True)
#sisse ringisse, välja ringist
time1=time1*10**(-6)
time2=time2*10**(-6)
voltage1=voltage1*10**(-3)
voltage2=voltage2*10**(-3)
voltagesisse=voltagesisse*10**(-3)
voltagevälja=voltagevälja*10**(-3)

##########################################################################################

#ajaline käitumine

#################################################################
#joonis üks

mid=np.mean(voltage1[1:len(voltage1)-3])

for i in range(len(voltage1)):
    if i%2==0:
        voltage1[i]=mid+mid-voltage1[i]
    voltage1[i]=voltage1[i]-mid

params1, cov1 = curve_fit(expo, time1, voltage1)
uncertainties = np.sqrt(np.diag(cov1))

for name, value, uncertainty in zip('abc', params1, uncertainties): 
    print(f'{name} = {value:8.3f} ± {uncertainty:.3f}')
a, b= params1[0], params1[1]
x = np.linspace(time1[0], time1[len(time1)-1], 100)
mu = ufloat(b, uncertainties[1])
tex = 1/(2*math.pi*mu)
reff = 4*math.pi*mu*L
rtheo = 2*(L/C)**(1/2)
print(f"mid=", mid)
print(f"Tex=", tex)
print(f"Reff=", reff)
print(f"Rtheo=", rtheo)


plt.figure(0)
plt.plot(time1, voltage1, 'c*', label= 'Messwerte')
plt.plot(x, a*np.exp(-2*math.pi*b*x), 'g-', label= 'gefittete Funktion')
plt.xlabel("Zeit [sec]")
plt.ylabel("U [V]")
plt.savefig('plot1.pdf')
plt.close()

#joonis kaks
plt.figure(1)
plt.plot(freq, voltagevälja/voltagesisse, 'g*', label= 'Messwerte')
plt.xlabel("f [kHz]")
plt.ylabel("$U_C/U_0$")
plt.savefig('plot2.pdf')
plt.close()