import matplotlib.pyplot as plt
import math
import numpy as np
from uncertainties import ufloat
import scipy.constants as const
from scipy.stats import sem

c = 2.9979*10**8
d = 5.85*10**(-2)

# Aufgabe 1
##########################################################################################
alpha = np.array([20.0, 30.0, 40.0, 45.0, 50.5, 60.0, 70.0])
beta = np.array([20.5, 31.0, 41.0, 46, 51.0, 60.5, 71.0])
x=np.linspace(np.min(alpha), np.max(alpha))
params,covariance_matrix=np.polyfit(alpha, beta, deg=1, cov=True)
errors = np.sqrt(np.diag(covariance_matrix))
for name, value, error in zip('ab', params, errors):
    print(f'{name} = {value:.3f} ± {error:.3f}')

plt.figure(0)
plt.plot(x, params[1]+x*params[0], "c", label="lineare Regression")
plt.plot(alpha ,beta, '*', label="Messdaten")
plt.xlabel(r"$\alpha_1$ [rad]")
plt.ylabel(r"$\alpha_2$ [rad]")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('plot.pdf')


# Aufgabe 2
#############################################################################################
alpha = np.array([10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0])
beta = np.array([7.25, 13.5, 20, 26, 31.5, 36, 39.5])
n = np.sin(alpha*math.pi/360.0)/np.sin(beta*math.pi/360.0)
nmid = np.mean(n)
nerr = sem(n)
nfloat = ufloat(nmid, nerr)
cm=c/nfloat

print(f"Brechungsindizes:", n)
print(f"Mittelwert:", nfloat)
print(f"Lichtgeschwindigkeit:", cm)



# Aufgabe 3
##########################################################################################

s = d*np.sin((alpha-beta)*math.pi/360)/np.cos(beta*math.pi/360)
betath = np.sin(alpha*math.pi/360)/nmid
print(betath)
print(beta)
sth = d*np.sin((alpha-betath)*math.pi/360)/np.cos(betath*math.pi/360)

print(f"s-Wert:", s)
print(f"theoretischer S-Wert:", sth)



# Aufgabe 4
##############################################################################################

alpha = np.array([26, 30, 40, 50, 60])
alphag = np.array([90, 74.5, 57.5, 45, 35])
alphar = np.array([87.5, 73.5, 57, 44.75, 35])

beta1 = np.arcsin(np.sin(alpha)/nmid)
beta2g = np.arcsin(np.sin(alphag)/nmid)
beta2r = np.arcsin(np.sin(alphar)/nmid)

delg = (alpha + alphag) - (beta1 + beta2g)
delr = (alpha + alphar) - (beta1 + beta2r)

print(f"Delta-Funktion grün:", delg)
print(f"Delta-Funktion rot:", delr)



#Aufgabe 5
##################################################################################

lambd=np.zeros(14)


d=0.001/600
phi = np.array([23.2, 22.7])
for i in range(2):
    lambd[i]=d*np.sin(phi[i])

d=0.001/300
phi = np.array([[11, 11.5],[22.7, 23.2]])
for i in range(2):
    lambd[i+2]=d*np.sin(phi[i][0])/(i+1)
    lambd[i+4]=d*np.sin(phi[i][1])/(i+1)

d=0.001/100
phi = np.array([[4, 3.6],[7.7, 7],[11.5, 11],[15, 15]])
for i in range(4):
    lambd[i+6]=d*np.sin(phi[i][0])/(i+1)
    lambd[i+10]=d*np.sin(phi[i][1])/(i+1)
print(f"Lambda bei 100mm:", lambd)

mid = np.mean(lambd)
err = sem(lambd)
midl = ufloat(mid, err)
print("Mittelwert von lambda:", midl)