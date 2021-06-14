import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
from scipy import stats
from scipy.optimize import curve_fit
import scipy.constants as const
from tabulate import tabulate

#####################################################################################

#a)
print('###########################################################')
print('a)') 
R3 = np.array([390, 490, 242])
R2 = np.array([500, 332, 1000])
R4 = 1000-R3

R_x=R2*R3/R4
table ={'1': R2, '3': R3, '134': R4, '10': R_x}
print(tabulate (table, tablefmt="latex"))
R_mid=np.mean(R_x)
R_err=np.std(R_x)
print(f"R_mid: ", R_mid, " R_err: ", R_err)

R_theo = 319.5
R_abw = (R_mid - R_theo)/R_theo

#######################################################################################

#b)
print('###########################################################')
print('b)') 
R3 = 632
R2 = 228.5
R4 = 1000-R3
C2 = 750*10**(-9)

R_x=R2*R3/R4
print(f"R_x: ", R_x)
C_x=C2*R4/R3
print(f"C_x: ", C_x)

R_theo = 464.9
R_abw = (R_x - R_theo)/R_theo
C_theo = 433.7
C_abw = (C_x - C_theo)/C_theo

#########################################################################################

#c)
print('###########################################################')
print('c)') 
R3 = 608
R2 = 224
R4 = 1000-R3
L2 = 0.0201

R_x=R2*R3/R4
print(f"R_x: ", R_x)
L_x=L2*R4/R3
print(f"L_x: ", L_x)

R_theo = 360.5
R_abw = (R_x - R_theo)/R_theo
L_theo = 433.7
L_abw = (L_x - L_theo)/L_theo

#########################################################################################

#d)
print('###########################################################')
print('d)') 
R3 = 191
R2 = 332
R4 = 179
C2 = 750 * 10**(-6)

R_x=R2*R3/R4
print(f"R_x: ", R_x)
L_x=L2*R2*R3
print(f"L_x: ", L_x)

R_theo = 360.5
R_abw = (R_x - R_theo)/R_theo
L_theo = 433.7
L_abw = (L_x - L_theo)/L_theo

###########################################################################################

#e)
print('###########################################################')
print('e)') 

v, U= np.genfromtxt('data.txt',unpack=True)
U=U/2
table ={'1': v, '2':U}
print(tabulate (table, tablefmt="latex"))
U=U/10
R=1000
C=295*10**(-9)
v_0 = 1/(R*C)
print(v_0)

def hui(omma, a):
    return np.sqrt(((omma/a)**2-1)**2/(9*((1-(omma/a)**2)**2+ 9* (omma/a)**2)))

v_0, cov = curve_fit(hui, v, U, p0=(v_0))
print(v_0)
x=np.linspace(v[0], v[20], 10000)
y=hui(x, v_0)

plt.figure()
plt.plot(v/v_0, U, 'b*', label='Messdaten')
plt.plot(x/v_0, y, 'c')
plt.xscale('log')
plt.xlabel(r'$v/v_0$')
plt.ylabel(r'$U_{Br}/U_{Sp}$')
plt.legend(loc='best')
plt.savefig('build/plot.pdf')

#Klirrfaktor
print(f"Klirrfaktor: ", U[8]/hui(2, 1))