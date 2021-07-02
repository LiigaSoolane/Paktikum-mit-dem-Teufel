import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
from scipy import stats
from scipy.optimize import curve_fit
import scipy.constants as const
from tabulate import tabulate

###############################################################################

#funktsioonid

def einzel(phi,b, a):
    lamda=633*10**(-9)
    uff=(a**2)*(b**2)*((lamda/(np.pi*b*np.sin(phi)))**2)
    teil2=(np.sin((np.pi*b*np.sin(phi))/lamda))**2
    return(uff*teil2)

def doppel(phi,b,A_0,s):
    lamda=633*10**(-9)
    teil1=(A_0**2)*np.cos((np.pi*s*np.sin(phi)/lamda))**2
    teil2=(lamda/(np.pi*b*np.sin(phi)))**2
    teil3=np.sin((np.pi*b*np.sin(phi))/lamda)**2
    return(teil1*teil2*teil3)


def a(ag,at):
    return(np.absolute(ag-at)/(at))

#######################################################################################

#Einzelspalt

#reading data
L = 1.047
I_d = 2.27*10**(-9)

d, I= np.genfromtxt('data.txt',unpack=True)
d = (d*10**(-3))-0.026
I = (I*10**(-3))-I_d
phi=(d)/L



#走吧
table ={'1': d, '3': I}
print(tabulate (table, tablefmt="latex"))
I1=I[1:50]
phi1=phi[1:50]

params, covariance =curve_fit(einzel,phi1,I1,p0=[0.00015,70])
errors = np.sqrt(np.diag(covariance))
print('breite des Spaltes:', params[0],'+-',errors[0])
print('Abweichung', a(params[0],0.00015))
print('A_0:', params[1],'+-',errors[1])
#def einzel(phi,b,lamda,A_0)

x=np.linspace(-0.025,0.025,10000)

plt.figure(1)
plt.plot(phi,I,'cx',label=r'$Messwerte$')
plt.plot(x,einzel(x,*params),'b-',label=r'$Ausgleichsfunktion$')
plt.legend(loc='best')
plt.xlabel('Abstand vom Hauptmaximum [m]')
plt.ylabel(r'$\mathrm{Intensität \ I [A]}$')
plt.savefig('plot.pdf')


###############################################################################

#Doppelspalt

#reading data
d, I_2= np.genfromtxt('ddata.txt',unpack=True)
d = (d*10**(-3))-0.026
I_2 = (I_2*10**(-3))-I_d
phi=d/L

#走吧
table ={'1': d, '3': I}
print(tabulate (table, tablefmt="latex"))
I2=I_2[1:50]
phi2=phi[1:50]

params2, covariance2 =curve_fit(doppel,phi2,I2,p0=[0.000015,1, 0.000065])
errors2 = np.sqrt(np.diag(covariance2))
print('breite des Spaltes:', params2[0],'+-',errors2[0])
print('Abweichung', a(params2[0],0.00015))
print('A_0:', params2[1],'+-',errors2[1])
print(params2[2], '+-', errors2[2])


x=np.linspace(-0.025,0.025,1000)

plt.figure(2)
plt.plot(phi,I_2,'bx',label=r'$Messwerte$')
plt.plot(x,doppel(x,0.000138, 0.0543, 0.000498),'g-',label=r'$Ausgleichsfunktion$')
plt.legend(loc='best')
plt.xlabel('Abstand vom Hauptmaximum [m]')
plt.ylabel(r'$\mathrm{Intensität \ I [A]}$')
plt.savefig('plot1.pdf')

########################################################################################

#overlay

plt.figure(3)
plt.plot(phi,I_2,'bx',label=r'Messwerte des Doppelspaltes')
plt.plot(x,doppel(x, 0.000138, 0.0543, 0.000498),'g-',label=r'Ausgleichsfunktion des Doppelspaltes')
plt.plot(phi,I,'cx',label=r'Messwerte des Einzelspaltes')
plt.plot(x,einzel(x,*params),'b-',label=r'Ausgleichsfunktion des Einzelspaltes')
plt.legend(loc='best')
plt.xlabel('Abstand vom Hauptmaximum [m]')
plt.ylabel('Intensität I [A]')
plt.savefig('plot2.pdf')


###############################################################################